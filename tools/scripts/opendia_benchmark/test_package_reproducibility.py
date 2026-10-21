#!/usr/bin/env python3
"""Focused stdlib tests for local reproducibility bundles."""

import json
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import tempfile
import unittest

import package_reproducibility as bundle


class PackageTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.base = Path(self.temp.name)
        self.root = self.base / 'artifacts'
        self.output = self.base / 'bundle'
        self.dataset = self.base / 'timsim-input'
        self.astral = self.base / 'astral-input'
        self.write(self.root / 'build-manifest.txt', 'source_dirty=true\n')
        self.write(self.root / 'timsim/meta/raw_files.txt', '/original/raw.d\n')
        for job in bundle.JOBS:
            for name, text in {
                'meta/exit_code.txt': '0\n',
                'meta/command.sh': 'env ORIGINAL=1 /original/OpenDIA "a b"\n',
                'meta/source_revision.txt': 'abc123\n',
                'meta/source_status.txt': ' M tracked.cpp\n',
                'meta/OpenDIA.sha256': 'a' * 64 + '  /original/OpenDIA\n',
                'logs/host_workload.json': '{"abort_reason": null}',
                'logs/resources.txt': 'Exit status: 0\n',
                'logs/stdout.log': 'done\n',
                'results/OpenDIA.results.tsv': 'header\nresult\n',
                'results/OpenDIA.protein.matrix.tsv': 'matrix\n',
                'intermediates/large.osw': 'exclude',
                'results/large.raw': 'exclude',
            }.items():
                self.write(self.root / job / name, text)
        self.write(self.root / 'timsim/validation.json', '{}')
        self.write(self.root / 'astral/metrics.json', '{}')
        self.write(self.root / 'report/FINAL_REPORT.md', '# Report\n')
        records = []
        for job, name in zip(bundle.JOBS[:2], ('target.tsv', 'entrapment.tsv')):
            path = self.dataset / 'libraries' / name
            self.write(path, job)
            sha = bundle.digest(path)
            records.append(f'{sha}  ./libraries/{name}\n')
            self.write(self.root / job / 'meta/library.sha256',
                       f'{sha}  /original/libraries/{name}\n')
        self.write(self.dataset / 'metadata/fixture_manifest.json', '{}')
        self.write(self.dataset / 'SHA256SUMS', ''.join(records) +
                   '0' * 64 + '  ./README.md\n')
        for name, location in (('raw_files', 'raw/run.raw'), ('fasta', 'db.fasta')):
            path = self.astral / location
            self.write(path, location)
            self.write(self.root / f'astral/meta/{name}.sha256',
                       f'{bundle.digest(path)}  /original/{location}\n')

    @staticmethod
    def write(path, text):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text)

    def package(self, **kwargs):
        return bundle.package(self.root, self.output, **kwargs)

    def test_roundtrip_portable_inputs_archive_and_snapshot(self):
        archive = self.base / 'bundle.tar.gz'
        self.package(timsim_dataset=self.dataset, astral_dataset=self.astral,
                     astral_fasta=self.astral / 'db.fasta', include_report=True,
                     archive=archive)
        bundle.verify(self.output)
        manifest = bundle.read_json(self.output / 'manifest.json')
        self.assertIsNone(manifest['ci_provenance'])
        self.assertTrue(all(j['dirty'] for j in manifest['jobs'].values()))
        self.assertTrue(all(manifest['inputs_verified_now'].values()))
        self.assertEqual(manifest['scripts_snapshot'], 'packaging-time')
        self.assertTrue((self.output / 'inputs/timsim/SHA256SUMS').is_file())
        self.assertTrue((self.output / 'report/FINAL_REPORT.md').is_file())
        for job in bundle.JOBS:
            self.assertEqual((self.output / job / 'meta/command.sh').read_bytes(),
                             (self.root / job / 'meta/command.sh').read_bytes())
        self.assertFalse(any('intermediates/' in p for p in manifest['files']))
        self.assertFalse(any(p.endswith('.raw') for p in manifest['files']))
        self.assertNotIn('/original/', (self.output / 'inputs.json').read_text())
        for script in bundle.SCRIPTS:
            self.assertTrue((self.output / 'scripts' / script).is_file())
        with tarfile.open(archive) as tar:
            self.assertTrue(all(m.isfile() or m.isdir() for m in tar.getmembers()))
            self.assertIn('bundle/manifest.json', tar.getnames())
        moved = self.base / 'moved'
        shutil.move(self.output, moved)
        result = subprocess.run([sys.executable,
                                 str(moved / 'scripts/package_reproducibility.py'),
                                 'verify', str(moved)], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_nonzero_missing_and_aborted_jobs(self):
        for path, text in (
            ('meta/exit_code.txt', '1'),
            ('logs/host_workload.json', '{"abort_reason":"disk"}'),
            ('logs/host_workload.json', '{}'),
            ('logs/resources.txt', 'Exit status: 9\n'),
            ('logs/resources.txt', ''),
            ('logs/resources.txt', 'Command terminated by signal 9\nExit status: 0\n'),
            ('logs/resources.txt', 'Exit status: 0\nExit status: 0\n'),
            ('results/OpenDIA.results.tsv', ''),
        ):
            with self.subTest(path=path, text=text):
                target = self.root / bundle.JOBS[0] / path
                old = target.read_text()
                target.write_text(text)
                with self.assertRaises(ValueError):
                    self.package()
                self.assertFalse(self.output.exists())
                target.write_text(old)
        (self.root / 'astral/meta/exit_code.txt').unlink()
        with self.assertRaises(OSError):
            self.package()

    def test_never_overwrites(self):
        self.output.mkdir()
        self.write(self.output / 'keep', 'untouched')
        with self.assertRaises(ValueError):
            self.package()
        self.assertEqual((self.output / 'keep').read_text(), 'untouched')
        with self.assertRaises(ValueError):
            bundle.package(self.root, self.root / 'nested')

    def test_symlink_in_payload_or_parent(self):
        link = self.root / 'astral/logs/link'
        link.symlink_to(self.root / 'build-manifest.txt')
        with self.assertRaises(ValueError):
            self.package()
        link.unlink()
        parent = self.base / 'linked'
        parent.symlink_to(self.root, target_is_directory=True)
        with self.assertRaises(ValueError):
            bundle.package(parent, self.output)

    def test_input_mismatch_and_traversal(self):
        self.write(self.astral / 'raw/run.raw', 'changed')
        with self.assertRaisesRegex(ValueError, 'RAW checksum mismatch'):
            self.package(astral_dataset=self.astral)
        for path in ('../outside', '/absolute', './a/../b', 'a\\b', '.'):
            with self.subTest(path=path), self.assertRaises(ValueError):
                bundle.relative_name(path)
        self.write(self.dataset / 'SHA256SUMS', 'a' * 64 + '  ../outside\n')
        with self.assertRaises(ValueError):
            self.package(timsim_dataset=self.dataset)
        self.assertFalse(self.output.exists())

    def test_verify_rejects_corruption_extra_missing_and_manifest_traversal(self):
        self.package()
        file = self.output / 'astral/results/OpenDIA.results.tsv'
        old = file.read_bytes()
        file.write_bytes(b'x' * len(old))
        with self.assertRaisesRegex(ValueError, 'Size/hash mismatch'):
            bundle.verify(self.output)
        file.write_bytes(old + b'x')
        with self.assertRaises(ValueError):
            bundle.verify(self.output)
        file.write_bytes(old)
        extra = self.output / 'extra'
        extra.write_text('extra')
        with self.assertRaises(ValueError):
            bundle.verify(self.output)
        extra.unlink()
        file.unlink()
        with self.assertRaises((OSError, ValueError)):
            bundle.verify(self.output)
        file.write_bytes(old)
        manifest = bundle.read_json(self.output / 'manifest.json')
        manifest['files']['../outside'] = {'size': 1, 'sha256': 'a' * 64}
        (self.output / 'manifest.json').write_text(json.dumps(manifest))
        with self.assertRaises(ValueError):
            bundle.verify(self.output)

    def test_verify_symlink_and_checksum_tampering(self):
        self.package()
        file = self.output / 'astral/results/OpenDIA.results.tsv'
        file.unlink()
        file.symlink_to(self.root / 'astral/results/OpenDIA.results.tsv')
        with self.assertRaises(ValueError):
            bundle.verify(self.output)
        file.unlink()
        shutil.copyfile(self.root / 'astral/results/OpenDIA.results.tsv', file)
        sums = self.output / 'SHA256SUMS'
        sums.write_text(sums.read_text().replace('manifest.json', '../manifest.json'))
        with self.assertRaises(ValueError):
            bundle.verify(self.output)

    def test_input_fasta_library_and_fixture_mismatches(self):
        self.write(self.astral / 'db.fasta', 'changed')
        with self.assertRaisesRegex(ValueError, 'FASTA checksum mismatch'):
            self.package(astral_fasta=self.astral / 'db.fasta')
        library = self.dataset / 'libraries/target.tsv'
        library.write_text('changed')
        with self.assertRaisesRegex(ValueError, 'TimSim checksum mismatch'):
            self.package(timsim_dataset=self.dataset)
        library.write_text(bundle.JOBS[0])
        captured = self.root / bundle.JOBS[0] / 'meta/library.sha256'
        captured.write_text('0' * 64 + '  /original/libraries/target.tsv\n')
        with self.assertRaisesRegex(ValueError, 'differs from captured input'):
            self.package(timsim_dataset=self.dataset)

    def test_existing_archive_is_preserved(self):
        archive = self.base / 'existing.tar.gz'
        archive.write_bytes(b'keep')
        with self.assertRaises(ValueError):
            self.package(archive=archive)
        self.assertEqual(archive.read_bytes(), b'keep')
        self.assertFalse(self.output.exists())

    def test_optional_evidence_and_symlink_rejection(self):
        evidence = self.base / 'supplement'
        self.write(evidence / 'nested/tolerances.json', '{"predeclared": true}')
        self.package(evidence=evidence)
        self.assertEqual((self.output / 'evidence/nested/tolerances.json').read_bytes(),
                         (evidence / 'nested/tolerances.json').read_bytes())
        bundle.verify(self.output)
        shutil.rmtree(self.output)
        (evidence / 'bad-link').symlink_to(self.root / 'build-manifest.txt')
        with self.assertRaises(ValueError):
            self.package(evidence=evidence)
        self.assertFalse(self.output.exists())

    def test_captured_shell_is_never_executed(self):
        sentinel = self.base / 'must-not-exist'
        command = self.root / 'astral/meta/command.sh'
        command.write_text(f'touch {sentinel}\n')
        self.package()
        bundle.verify(self.output)
        self.assertFalse(sentinel.exists())
        self.assertFalse((self.output / 'report').exists())


if __name__ == '__main__':
    unittest.main()
