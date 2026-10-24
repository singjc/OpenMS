#!/usr/bin/env python3
"""Sample Linux host load during one benchmark, without adjusting measured times."""

import argparse
import json
import math
import os
import platform
import signal
import stat
import statistics
import subprocess
import time
from pathlib import Path


def system_info():
    """Record hardware, not serial numbers or host/user names."""
    def optional_text(path):
        try:
            return Path(path).read_text().strip() or None
        except OSError:
            return None

    memory = dict(line.split(':', 1) for line in Path('/proc/meminfo').read_text().splitlines())
    cpu_model = next((line.split(':', 1)[1].strip() for line in Path('/proc/cpuinfo').read_text().splitlines()
                      if line.startswith('model name')), None)
    return {'machine': ' '.join(filter(None, [optional_text('/sys/class/dmi/id/sys_vendor'),
                                             optional_text('/sys/class/dmi/id/product_name')])) or None,
            'cpu_model': cpu_model, 'os': platform.platform(),
            'logical_cpus': os.cpu_count(), 'affinity_logical_cpus': len(os.sched_getaffinity(0)),
            'mem_total_mib': int(memory['MemTotal'].split()[0]) / 1024,
            'mem_available_mib_at_start': int(memory['MemAvailable'].split()[0]) / 1024,
            'captured_at': time.time()}


def process_tree(pid):
    """Snapshot processes, including children created by non-main threads."""
    processes = {}
    pending = [pid]
    while pending:
        current = pending.pop()
        if current in processes:
            continue
        try:
            stat = Path(f"/proc/{current}/stat").read_text().rsplit(")", 1)[1].split()
            processes[current] = stat
            for task in Path(f"/proc/{current}/task").iterdir():
                try:
                    pending.extend(int(child) for child in (task / "children").read_text().split())
                except (FileNotFoundError, ProcessLookupError, PermissionError):
                    continue
        except (FileNotFoundError, ProcessLookupError, PermissionError):
            continue
    return processes


def process_ticks(pid):
    """Estimate cumulative tree CPU, retaining CPU of waited children."""
    processes = process_tree(pid)
    # cutime/cstime retain reaped child work. Snapshot races still make this an estimate.
    ticks = sum(sum(int(stat[index]) for index in (11, 12, 13, 14)) for stat in processes.values())
    descendants = [child for child in processes if child != pid]
    return ticks, descendants


def stop_descendants(pid, grace_seconds=2):
    """Stop only the measured tree, retaining pidfds if children become orphans."""
    tracked = {}
    deadline = time.monotonic() + grace_seconds
    initial = process_tree(pid)
    root_started = initial[pid][19] if pid in initial else None
    try:
        while True:
            # Pin all identities before signalling a parent that could orphan its children.
            processes = process_tree(pid)
            if pid not in processes or processes[pid][19] != root_started:
                processes = {}
            # Orphans may themselves create children during their TERM cleanup.
            for child, started in list(tracked):
                subtree = process_tree(child)
                if child in subtree and subtree[child][19] == started:
                    processes.update(subtree)
            new_descriptors = []
            for child, stat in processes.items():
                identity = (child, stat[19])  # PID plus kernel start time
                if child == pid or stat[0] == "Z" or identity in tracked:
                    continue
                try:
                    descriptor = os.pidfd_open(child)
                except ProcessLookupError:
                    continue
                try:
                    current = Path(f"/proc/{child}/stat").read_text().rsplit(")", 1)[1].split()
                    if current[19] != stat[19]:
                        os.close(descriptor)
                        continue
                except (FileNotFoundError, ProcessLookupError):
                    os.close(descriptor)
                    continue
                tracked[identity] = descriptor
                new_descriptors.append(descriptor)
            for descriptor in reversed(new_descriptors):
                try:
                    signal.pidfd_send_signal(descriptor, signal.SIGTERM)
                except ProcessLookupError:
                    pass

            alive = []
            for (child, started), descriptor in tracked.items():
                try:
                    stat = Path(f"/proc/{child}/stat").read_text().rsplit(")", 1)[1].split()
                    if stat[19] == started and stat[0] != "Z":
                        alive.append(descriptor)
                except (FileNotFoundError, ProcessLookupError):
                    continue
            if not alive and tracked:
                return True
            if not alive and not Path(f"/proc/{pid}").exists():
                return True
            now = time.monotonic()
            if now >= deadline:
                for descriptor in alive:
                    try:
                        signal.pidfd_send_signal(descriptor, signal.SIGKILL)
                    except ProcessLookupError:
                        pass
                if now >= deadline + 2:
                    return not alive
            time.sleep(0.05)
    finally:
        for descriptor in tracked.values():
            os.close(descriptor)


def pressure(kind):
    try:
        first = Path(f"/proc/pressure/{kind}").read_text().splitlines()[0]
        return float(dict(item.split("=", 1) for item in first.split()[1:])["avg10"])
    except (OSError, KeyError, IndexError):
        return None


def job_memory(pid):
    """Sample simultaneous RSS and largest observed individual VmHWM, not their sum."""
    rss, high_water = [], []
    for child in process_tree(pid):
        try:
            fields = dict(line.split(':', 1) for line in Path(f'/proc/{child}/status').read_text().splitlines())
            if 'VmRSS' in fields:
                rss.append(int(fields['VmRSS'].split()[0]) * 1024)
            if 'VmHWM' in fields:
                high_water.append(int(fields['VmHWM'].split()[0]) * 1024)
        except (OSError, ValueError):
            continue
    return {'job_tree_rss_bytes_sampled': sum(rss) if rss else None,
            'job_largest_process_hwm_bytes_observed': max(high_water) if high_water else None}


def disk_bytes(path):
    """Apparent regular-file bytes, excluding symlinks and deduplicating hardlinks."""
    total = 0
    seen = set()
    for directory, _, files in os.walk(path):
        for name in files:
            try:
                item = (Path(directory) / name).lstat()
                identity = (item.st_dev, item.st_ino)
                if stat.S_ISREG(item.st_mode) and identity not in seen:
                    total += item.st_size
                    seen.add(identity)
            except OSError:
                continue
    return total


def cgroup_memory(path, pid):
    """Read only a requested dedicated scope containing a measured descendant."""
    if path is None:
        return {}
    try:
        relative = '/' + str(path.resolve().relative_to('/sys/fs/cgroup'))
        if not path.name.startswith('opendia-bench-') or path.suffix != '.scope':
            return {}
        if not any(f'0::{relative}' in Path(f'/proc/{child}/cgroup').read_text().splitlines()
                   for child in process_tree(pid)):
            return {}
        values = {}
        for filename in ('memory.current', 'memory.peak', 'memory.swap.current'):
            key = 'cgroup_' + filename.replace('.', '_') + '_bytes'
            try:
                values[key] = int((path / filename).read_text())
            except (OSError, ValueError):
                values[key] = None
        return values
    except (OSError, ValueError):
        return {}


def host_sample(pid, output):
    cpu = [int(value) for value in Path("/proc/stat").read_text().splitlines()[0].split()[1:9]]
    memory = {}
    for line in Path("/proc/meminfo").read_text().splitlines():
        key, value = line.split(":", 1)
        memory[key] = int(value.split()[0])
    ticks, descendants = process_ticks(pid)
    sample = {
        "timestamp": time.time(),
        "load1": os.getloadavg()[0],
        "mem_available_mib": memory["MemAvailable"] / 1024,
        "swap_used_mib": (memory["SwapTotal"] - memory["SwapFree"]) / 1024,
        "cpu_pressure_some_percent": pressure("cpu"),
        "memory_pressure_some_percent": pressure("memory"),
        "io_pressure_some_percent": pressure("io"),
        "disk_free_gib": os.statvfs(output).f_bavail * os.statvfs(output).f_frsize / 1024**3,
    }
    return sample, cpu, ticks, descendants


def scope_workload_finished(path, done):
    """Only the keeper and its direct polling sleep may remain in the scope."""
    try:
        keeper = int(done.read_text().split()[1])
        for process in (path / 'cgroup.procs').read_text().split():
            child = int(process)
            if child == keeper:
                continue
            try:
                fields = Path(f'/proc/{child}/stat').read_text().rsplit(')', 1)[1].split()
                if fields[0] == 'Z':
                    continue
                if int(fields[1]) == keeper and Path(f'/proc/{child}/comm').read_text().strip() == 'sleep':
                    continue
                return False
            except FileNotFoundError:
                continue
        return True
    except (OSError, ValueError, IndexError):
        return False


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pid", type=int)
    parser.add_argument("output", type=Path)
    parser.add_argument("--interval", type=float, default=1)
    parser.add_argument("--cgroup", type=Path, help="Explicit dedicated opendia-bench-*.scope only")
    parser.add_argument("--cgroup-done", type=Path)
    parser.add_argument("--cgroup-release", type=Path)
    parser.add_argument("--min-free-disk-gib", type=float, default=3)
    args = parser.parse_args()
    if bool(args.cgroup) != bool(args.cgroup_done and args.cgroup_release):
        parser.error('cgroup measurement requires both done and release paths')
    if args.pid <= 0 or not math.isfinite(args.interval) or not math.isfinite(args.min_free_disk_gib) or args.interval <= 0 or args.min_free_disk_gib < 0:
        parser.error("interval must be positive and disk reserve must not be negative")
    samples = []
    previous_cpu = previous_ticks = previous_time = None
    cores = os.cpu_count() or 1
    hardware = system_info()
    abort_reason = None
    cgroup_final = None
    with (args.output / "host_workload.jsonl").open("w", encoding="utf-8") as handle:
        while Path(f"/proc/{args.pid}").exists():
            root = process_tree(args.pid).get(args.pid)
            if root is None or root[0] == 'Z':
                break
            sample, cpu, ticks, descendants = host_sample(args.pid, args.output)
            sample.update(job_memory(args.pid))
            sample.update(cgroup_memory(args.cgroup, args.pid))
            sample['output_apparent_bytes_sampled'] = disk_bytes(args.output.parent)
            sample['temp_apparent_bytes_sampled'] = disk_bytes(args.output.parent / 'tmp')
            now = time.monotonic()
            if previous_cpu is not None:
                elapsed = now - previous_time
                total = sum(cpu) - sum(previous_cpu)
                idle = cpu[3] + cpu[4] - previous_cpu[3] - previous_cpu[4]
                sample["host_cpu_busy_percent"] = 100 * (1 - idle / total) if total > 0 else None
                # Idle only, excluding I/O wait. Thread equivalents are not reserved capacity.
                idle_ticks = cpu[3] - previous_cpu[3]
                sample['host_cpu_idle_percent'] = 100 * idle_ticks / total if total > 0 else None
                sample['host_cpu_idle_thread_equivalents'] = cores * idle_ticks / total if total > 0 else None
                # This is a workload estimate, not a correction to elapsed time.
                job_cpu = 100 * (ticks - previous_ticks) / os.sysconf("SC_CLK_TCK") / elapsed / cores
                sample["job_cpu_host_percent"] = job_cpu if job_cpu >= 0 else None
                sample["other_cpu_busy_estimate_percent"] = (
                    max(0, sample["host_cpu_busy_percent"] - job_cpu)
                    if job_cpu >= 0 and total > 0 else None
                )
            samples.append(sample)
            if len(samples) % 60 == 0:
                print("JCODE_PROGRESS " + json.dumps({
                    "message": f"{args.output.parent.name}: {len(samples)} workload samples, {sample['disk_free_gib']:.1f} GiB disk free",
                }), flush=True)
            handle.write(json.dumps(sample, allow_nan=False) + "\n")
            handle.flush()
            if args.cgroup and args.cgroup_done.exists() and cgroup_final is None:
                verified = cgroup_memory(args.cgroup, args.pid)
                finished = scope_workload_finished(args.cgroup, args.cgroup_done)
                # Read the authoritative peak after confirming all workload children exited.
                cgroup_final = cgroup_memory(args.cgroup, args.pid) if finished else {}
                if cgroup_final.get('cgroup_memory_peak_bytes') is None:
                    abort_reason = 'Requested cgroup final memory peak could not be verified'
                if verified and not finished:
                    abort_reason = 'Workload processes survived in the verified dedicated scope'
                    try:
                        stopped = subprocess.run(['systemctl', '--user', 'kill', '--kill-whom=all',
                                                  args.cgroup.name], capture_output=True, timeout=5)
                        if stopped.returncode:
                            abort_reason += '; scope cleanup failed'
                    except (OSError, subprocess.TimeoutExpired):
                        abort_reason += '; scope cleanup failed'
                # Release even on failure, so the workload exit and GNU time survive.
                args.cgroup_release.write_text('captured\n' if not abort_reason else 'failed\n')
            if sample["disk_free_gib"] < args.min_free_disk_gib:
                abort_reason = f"Free disk space fell below {args.min_free_disk_gib:g} GiB"
                if not descendants:
                    time.sleep(0.1)
                    continue
                # Leave GNU time alive to record failure, but do not abandon resistant children.
                if not stop_descendants(args.pid):
                    abort_reason += "; some descendants did not stop after SIGKILL"
                break
            previous_cpu, previous_ticks, previous_time = cpu, ticks, now
            time.sleep(args.interval)
    if args.cgroup and not cgroup_final:
        abort_reason = abort_reason or 'Requested cgroup measurement never reached a verified final read'
    summary = {
        "sample_count": len(samples),
        "logical_cpus": cores,
        "system": hardware,
        "started_at": samples[0]['timestamp'] if samples else None,
        "finished_at": samples[-1]['timestamp'] if samples else None,
        "interval_seconds": args.interval,
        "abort_reason": abort_reason,
        "resource_schema_version": 1,
        "cgroup_path_requested": str(args.cgroup) if args.cgroup else None,
        "cgroup_final": cgroup_final,
        "measurement_notes": {
            "job_tree_rss_bytes_sampled": "Sum of resident process RSS at each sample, shared pages may be counted more than once. Short-lived children may be missed.",
            "job_largest_process_hwm_bytes_observed": "Largest individual VmHWM observed, never a sum of process peaks.",
            "cgroup_memory_peak_bytes": "Kernel memory.peak samples include file cache and keeper overhead, not process RSS. cgroup_final is read after workload completion while dedicated scope keeper remains alive, before release.",
            "disk": "Sampled apparent regular-file bytes, not allocated blocks. Temp is a subset of output. Open unlinked files and external temp paths are not counted.",
        },
        "note": "Host CPU includes the benchmark. Job CPU and other CPU are estimates from process snapshots. Elapsed times are not corrected for workload.",
    }
    for key in set().union(*(sample.keys() for sample in samples)):
        if key == "timestamp":
            continue
        values = [sample[key] for sample in samples if sample.get(key) is not None]
        summary[key] = {"mean": statistics.fmean(values), "min": min(values), "max": max(values)} if values else None
    # CPU fields are absent in the first sample, so include them explicitly.
    for key in ("host_cpu_busy_percent", "host_cpu_idle_percent", "host_cpu_idle_thread_equivalents",
                "job_cpu_host_percent", "other_cpu_busy_estimate_percent"):
        values = [sample[key] for sample in samples if sample.get(key) is not None]
        summary[key] = {"mean": statistics.fmean(values), "min": min(values), "max": max(values)} if values else None
    (args.output / "host_workload.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    return 2 if abort_reason else 0


if __name__ == "__main__":
    raise SystemExit(main())
