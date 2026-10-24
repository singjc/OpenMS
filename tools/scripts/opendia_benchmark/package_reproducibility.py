#!/usr/bin/env python3
"""Local, stdlib-only packaging of completed OpenDIA benchmark evidence."""

import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import re
import shutil
import tarfile

JOBS = ('timsim/target_only', 'timsim/entrapment', 'astral')
SCRIPTS = ('lib.sh', 'build.sh', 'run_timsim.sh', 'run_astral.sh',
           'monitor_host.py', 'benchmark_metrics.py', 'validate_timsim.py',
           'metrics_astral.py', 'make_report.py', 'render_report.py',
           'package_reproducibility.py', 'compare_repeats.py', 'measure_outputs.py')


def safe_path(path):
    """Reject links in every component, including the supplied root's parents."""
    path = Path(path).absolute()
    for component in (path, *path.parents):
        if component.is_symlink():
            raise ValueError(f'Symlink not allowed: {component}')
    return path


def relative_name(name):
    path = PurePosixPath(name)
    if (not name or name == '.' or path.is_absolute() or '..' in path.parts or
            str(path) != name or '\\' in name or any(c in name for c in '\n\r\t')):
        raise ValueError(f'Unsafe relative path: {name!r}')
    return path


def digest(path):
    path = safe_path(path)
    if not path.is_file():
        raise ValueError(f'Not a regular file: {path}')
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def files_under(root):
    root = safe_path(root)
    if not root.is_dir():
        raise ValueError(f'Missing directory: {root}')
    for path in sorted(root.rglob('*')):
        safe_path(path)
        if not path.is_dir():
            if not path.is_file():
                raise ValueError(f'Not a regular file: {path}')
            relative_name(path.relative_to(root).as_posix())
            yield path


def read_json(path):
    return json.loads(safe_path(path).read_text())


def checksum_records(path, portable=False):
    records = []
    for line in safe_path(path).read_text().splitlines():
        match = re.fullmatch(r'([0-9a-fA-F]{64}) [ *](.+)', line)
        if not match:
            raise ValueError(f'Invalid checksum line in {path}')
        sha, name = match.groups()
        if portable:
            name = name.removeprefix('./')
            relative_name(name)
        else:
            if '..' in PurePosixPath(name).parts:
                raise ValueError(f'Unsafe captured checksum path: {name!r}')
            name = Path(name).name
        relative_name(name)
        records.append({'path': name, 'sha256': sha.lower()})
    if not records or len({r['path'] for r in records}) != len(records):
        raise ValueError(f'Empty or duplicate checksums: {path}')
    return records


REPLAY = '''# Replay and provenance limitations

This is local benchmark evidence, not a hermetic build or a CI attestation.
CI provenance was not captured. A dirty tree cannot be reproduced from its commit
alone. Source patches, dependencies, model assets and binaries are not bundled.
The scripts/ snapshot is from packaging time, not proven to be the run-time version.
Original-host commands are preserved verbatim in each job's meta/command.sh.
They contain host-specific paths and are evidence only. Do not execute them.
Input hashes in inputs.json are dataset-relative identifiers, not bundled inputs.
Supplied Astral inputs are checked against captured hashes, not copied into this
bundle. TimSim libraries are matched to captured hashes. Other TimSim input hashes
come from the supplied dataset manifest, not an independently captured run-time copy.
TimSim verification matches the runner's raw/, libraries/, truth/ and metadata/
scope. The full original SHA256SUMS is retained, including unchecked non-inputs.

On a compatible host, review scripts/ and the captured parameters first. From a
complete OpenMS checkout, copy the reviewed scripts snapshot into
`tools/scripts/opendia_benchmark/` there and use the existing runners with new
empty output paths (run the following from that directory):

```bash
export BUILD_DIR=/path/to/compatible/OpenMS-build
export OPEN_DIA="$BUILD_DIR/bin/OpenDIA" THREADS=4
# Review meta/resource_limits.txt and set BENCH_MEMORY_MAX,
# BENCH_MEMORY_SWAP_MAX, PERCOLATOR_THREADS and MIN_FREE_DISK_GIB as appropriate.
bash run_timsim.sh "$TIMSIM_DATASET" "$REPLAY_OUT/timsim" all
FASTA="$ASTRAL_FASTA" bash run_astral.sh "$ASTRAL_DATASET" "$REPLAY_OUT/astral"
```

Set TIMSIM_DATASET, ASTRAL_DATASET, ASTRAL_FASTA and REPLAY_OUT locally.
Runner defaults may differ from historical commands. Review Astral OXIDATION,
MAX_PEPTIDE_LENGTH, PREDICTION_THREADS and any extra captured options manually.
No captured shell is parsed or executed by this tool. Replay is not promised to
be bit-identical. Reference results are the captured final TSVs, which can include
post-search re-export. Recovery provenance is available only if captured in meta/logs.
Verify this directory with `python3 scripts/package_reproducibility.py verify .`.
SHA256SUMS covers manifest.json and every payload file, excluding SHA256SUMS itself.
These checks detect corruption, not malicious replacement of the entire bundle.
'''


def package(root, output, timsim_dataset=None, astral_dataset=None,
            astral_fasta=None, include_report=False, archive=None, evidence=None):
    root, output = safe_path(root), safe_path(output)
    if output.exists():
        raise ValueError(f'Destination already exists: {output}')
    if output.is_relative_to(root):
        raise ValueError('Destination must be outside the artifact root')
    archive = safe_path(archive) if archive else None
    if archive and (archive.exists() or archive.is_relative_to(output)):
        raise ValueError('Archive must be new and outside the bundle')
    selected = []

    def add_file(path, name):
        path = safe_path(path)
        if not path.is_file():
            raise ValueError(f'Missing required file: {path}')
        relative_name(name)
        selected.append((path, name))

    def add_tree(path, prefix):
        for file in files_under(path):
            add_file(file, f'{prefix}/{file.relative_to(path).as_posix()}')

    add_file(root / 'build-manifest.txt', 'build-manifest.txt')
    provenance = {}
    inputs = {}
    for job in JOBS:
        base = root / job
        if safe_path(base / 'meta/exit_code.txt').read_text().strip() != '0':
            raise ValueError(f'Unsuccessful job: {job}')
        monitor = read_json(base / 'logs/host_workload.json')
        if 'abort_reason' not in monitor or monitor['abort_reason']:
            raise ValueError(f'Monitor incomplete or aborted: {job}')
        for required in ('meta/command.sh', 'meta/source_revision.txt',
                         'meta/source_status.txt', 'meta/OpenDIA.sha256',
                         'logs/resources.txt', 'results/OpenDIA.results.tsv'):
            if not safe_path(base / required).is_file():
                raise ValueError(f'Missing required file: {job}/{required}')
        resources = safe_path(base / 'logs/resources.txt').read_text()
        measured_exits = re.findall(r'Exit status:\s*(\d+)', resources)
        if measured_exits != ['0'] or 'Command terminated by signal' in resources:
            raise ValueError(f'Nonzero measured exit: {job}')
        if not (base / 'results/OpenDIA.results.tsv').stat().st_size:
            raise ValueError(f'Empty reference result: {job}')
        for folder in ('meta', 'logs'):
            add_tree(base / folder, f'{job}/{folder}')
        # Explicitly exclude database/raw/intermediate files, even inside results.
        results = [p for p in files_under(base / 'results')
                   if p.parent == base / 'results' and
                   p.name.startswith('OpenDIA') and p.suffix == '.tsv']
        for result in results:
            add_file(result, f'{job}/results/{result.name}')
        provenance[job] = {
            'revision': (base / 'meta/source_revision.txt').read_text().strip(),
            'dirty': bool((base / 'meta/source_status.txt').read_text().strip()),
            'command_kind': 'original-host command, evidence only',
        }
    add_file(root / 'timsim/validation.json', 'timsim/validation.json')
    add_tree(root / 'timsim/meta', 'timsim/meta')
    add_file(root / 'astral/metrics.json', 'astral/metrics.json')
    for name in ('build.log', 'timsim-runner.log', 'astral-runner.log'):
        if (root / name).exists():
            add_file(root / name, name)
    for name in ('raw_files', 'fasta'):
        inputs[f'astral/{name}'] = checksum_records(root / f'astral/meta/{name}.sha256')
    for job in JOBS[:2]:
        inputs[f'{job}/library'] = checksum_records(root / job / 'meta/library.sha256')
    if timsim_dataset:
        dataset = safe_path(timsim_dataset)
        records = checksum_records(dataset / 'SHA256SUMS', portable=True)
        inputs['timsim/dataset'] = records
        add_file(dataset / 'SHA256SUMS', 'inputs/timsim/SHA256SUMS')
        add_file(dataset / 'metadata/fixture_manifest.json',
                 'inputs/timsim/metadata/fixture_manifest.json')
        for record in records:
            # Match run_timsim.sh: documentation and expected reports are not inputs.
            if PurePosixPath(record['path']).parts[0] not in {
                    'raw', 'libraries', 'truth', 'metadata'}:
                continue
            if digest(dataset / record['path']) != record['sha256']:
                raise ValueError(f'TimSim checksum mismatch: {record["path"]}')
        indexed = {r['path']: r['sha256'] for r in records}
        for job in JOBS[:2]:
            for record in inputs[f'{job}/library']:
                if indexed.get('libraries/' + record['path']) != record['sha256']:
                    raise ValueError(f'TimSim library differs from captured input: {job}')
    if astral_dataset:
        dataset = safe_path(astral_dataset)
        for record in inputs['astral/raw_files']:
            if digest(dataset / 'raw' / record['path']) != record['sha256']:
                raise ValueError(f'Astral RAW checksum mismatch: {record["path"]}')
    if astral_fasta:
        records = inputs['astral/fasta']
        if len(records) != 1 or digest(astral_fasta) != records[0]['sha256']:
            raise ValueError('Astral FASTA checksum mismatch')
    if include_report:
        for file in files_under(root / 'report'):
            if file.suffix in ('.md', '.pdf', '.svg', '.png'):
                add_file(file, f'report/{file.relative_to(root / "report").as_posix()}')
    if evidence:
        add_tree(safe_path(evidence), 'evidence')
    for name in SCRIPTS:
        add_file(Path(__file__).absolute().parent / name, f'scripts/{name}')
    output.mkdir(parents=True, exist_ok=False)
    try:
        for source, name in selected:
            destination = output / name
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(source, destination)
        (output / 'REPLAY.md').write_text(REPLAY)
        (output / 'inputs.json').write_text(json.dumps(inputs, indent=2) + '\n')
        inventory = {p.relative_to(output).as_posix():
                     {'size': p.stat().st_size, 'sha256': digest(p)}
                     for p in files_under(output)}
        manifest = {'schema_version': 1, 'jobs': provenance, 'files': inventory,
                    'ci_provenance': None, 'scripts_snapshot': 'packaging-time',
                    'limitations': ['Dirty source is not reconstructable from commit alone',
                                    'Inputs, binaries and intermediates are excluded',
                                    'Historical script and recovery provenance may be missing'],
                    'timsim_verification_scope': ['raw/', 'libraries/', 'truth/', 'metadata/'],
                    'inputs_verified_now': {'timsim': bool(timsim_dataset),
                                            'astral_raw': bool(astral_dataset),
                                            'astral_fasta': bool(astral_fasta)}}
        (output / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
        (output / 'SHA256SUMS').write_text(''.join(
            f'{digest(p)}  {p.relative_to(output).as_posix()}\n'
            for p in files_under(output) if p != output / 'SHA256SUMS'))
        verify(output)
        if archive:
            with archive.open('xb') as stream:
                with tarfile.open(fileobj=stream, mode='w:gz') as tar:
                    tar.add(output, arcname=output.name)
    except Exception:
        # Never remove an existing destination. This directory was created by us.
        shutil.rmtree(output)
        raise
    return output


def verify(root):
    root = safe_path(root)
    actual = {p.relative_to(root).as_posix() for p in files_under(root)}
    manifest = read_json(root / 'manifest.json')
    if manifest.get('schema_version') != 1:
        raise ValueError('Unsupported manifest schema')
    inventory = manifest['files']
    for name, record in inventory.items():
        relative_name(name)
        file = safe_path(root / name)
        if file.stat().st_size != record['size'] or digest(file) != record['sha256']:
            raise ValueError(f'Size/hash mismatch: {name}')
    expected = set(inventory) | {'manifest.json', 'SHA256SUMS'}
    if actual != expected:
        raise ValueError('Bundle has missing or unlisted files')
    checksums = checksum_records(root / 'SHA256SUMS', portable=True)
    if {r['path'] for r in checksums} != expected - {'SHA256SUMS'}:
        raise ValueError('SHA256SUMS inventory mismatch')
    for record in checksums:
        if digest(root / record['path']) != record['sha256']:
            raise ValueError(f'Checksum mismatch: {record["path"]}')


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='action', required=True)
    pack = sub.add_parser('package')
    pack.add_argument('--artifact-root', type=Path, required=True)
    pack.add_argument('--output', type=Path, required=True)
    pack.add_argument('--timsim-dataset', type=Path)
    pack.add_argument('--astral-dataset', type=Path, help='Dataset containing raw/*.raw')
    pack.add_argument('--astral-fasta', type=Path)
    pack.add_argument('--evidence', type=Path,
                      help='Optional local comparison/provenance supplement directory')
    pack.add_argument('--include-report', action='store_true')
    pack.add_argument('--archive', type=Path, help='Optional new .tar.gz path')
    check = sub.add_parser('verify')
    check.add_argument('bundle', type=Path)
    args = vars(parser.parse_args(argv))
    action = args.pop('action')
    try:
        if action == 'verify':
            verify(args['bundle'])
        else:
            args['root'] = args.pop('artifact_root')
            package(**args)
    except (OSError, ValueError, KeyError, TypeError) as error:
        parser.exit(1, f'Error: {error}\n')
    print('Bundle verified' if action == 'verify' else 'Bundle packaged and verified')


if __name__ == '__main__':
    main()
