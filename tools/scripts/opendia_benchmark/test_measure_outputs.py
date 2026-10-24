#!/usr/bin/env python3
"""Small synthetic output measurement tests, no benchmark execution."""

import hashlib
import io
import json
from pathlib import Path
import sqlite3
import tempfile
import unittest
from unittest.mock import patch
import zipfile

from measure_outputs import main, measure_artifact


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


class MeasurementTests(unittest.TestCase):
    def database(self, path):
        db = sqlite3.connect(path)
        db.executescript('CREATE TABLE FEATURE(ID INTEGER PRIMARY KEY);'
                         'INSERT INTO FEATURE VALUES (1), (2), (3);'
                         'CREATE TABLE FEATURE_TRANSITION(FEATURE_ID INTEGER);'
                         'INSERT INTO FEATURE_TRANSITION VALUES (1), (1), (2);'
                         'CREATE TABLE SCORE_MS2(FEATURE_ID INTEGER);')
        db.commit()
        return db

    def test_sqlite_counts_fanout_and_no_mutation(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'test.osw'
            self.database(path).close()
            before = digest(path)
            result = measure_artifact(path)
            self.assertEqual(result['tables']['FEATURE_TRANSITION'], {'status': 'ok', 'rows': 3})
            self.assertEqual(result['tables']['SCORE_MS2'], {'status': 'ok', 'rows': 0})
            self.assertEqual(result['tables']['PROTEIN'], {'status': 'absent', 'rows': None})
            self.assertEqual(result['feature_transition_fanout'], {
                'status': 'ok', 'feature_count': 3, 'zero_transition_features': 1,
                'min': 0, 'max': 2, 'mean': 1})
            self.assertEqual(result['final_bytes'], path.stat().st_size)
            self.assertTrue(all(item['status'] == 'absent' for item in result['sidecars']))
            self.assertEqual(before, digest(path))
            self.assertEqual(list(Path(directory).iterdir()), [path])

    def test_wal_visible_without_input_sidecar_writes(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'test.osw'
            db = self.database(path)
            db.execute('PRAGMA journal_mode=WAL')
            db.execute('INSERT INTO FEATURE_TRANSITION VALUES (3)')
            db.commit()
            files = list(Path(directory).iterdir())
            before = {file.name: digest(file) for file in files}
            result = measure_artifact(path)
            self.assertEqual(result['tables']['FEATURE_TRANSITION']['rows'], 4)
            self.assertEqual(before, {file.name: digest(file) for file in files})
            self.assertEqual(set(files), set(Path(directory).iterdir()))
            self.assertEqual(result['sidecars'][0]['status'], 'present')
            db.close()

    def test_empty_features_and_absent_input(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'test.osw'
            db = self.database(path)
            db.execute('DELETE FROM FEATURE')
            db.commit()
            db.close()
            result = measure_artifact(path)
            self.assertEqual(result['feature_transition_fanout']['feature_count'], 0)
            self.assertIsNone(result['feature_transition_fanout']['mean'])
            absent = Path(directory) / 'absent.osw'
            self.assertEqual(measure_artifact(absent)['status'], 'absent')
            self.assertFalse(absent.exists())

    def test_archive_missing_pyarrow_explicit_unavailable(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'test.oswpq'
            with zipfile.ZipFile(path, 'w') as archive:
                archive.writestr('../../features.parquet', b'not read without pyarrow')
            before = digest(path)
            with patch.dict('sys.modules', {'pyarrow.parquet': None}):
                result = measure_artifact(path)
            self.assertEqual(result['tables']['FEATURE'], {'status': 'unavailable', 'rows': None})
            self.assertEqual(result['tables']['FEATURE_TRANSITION']['status'], 'absent')
            self.assertEqual(result['parquet_members'][0]['status'], 'unavailable')
            self.assertEqual(before, digest(path))

    def test_native_parquet_metadata_and_partition_fanout(self):
        try:
            import pyarrow as pa
            import pyarrow.parquet as pq
        except ImportError:
            self.skipTest('optional pyarrow not installed')
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'test.oswpq'
            with zipfile.ZipFile(path, 'w') as archive:
                for run, feature_ids, transition_ids in ((1, [1, 2], [1, 1]), (2, [1], [1])):
                    for name, ids in (('features', feature_ids), ('feature_transition', transition_ids)):
                        data = io.BytesIO()
                        pq.write_table(pa.table({'feature_id': pa.array(ids, type=pa.int64())}), data)
                        archive.writestr(f'runs/run_id={run}/{name}.parquet', data.getvalue())
                data = io.BytesIO()
                pq.write_table(pa.table({'value': pa.array([], type=pa.int64())}), data)
                # Hostile path is counted by metadata but never extracted.
                archive.writestr('../../outside.parquet', data.getvalue())
            before = digest(path)
            result = measure_artifact(path)
            self.assertEqual(result['tables']['FEATURE']['rows'], 3)
            self.assertEqual(result['tables']['FEATURE_TRANSITION']['rows'], 3)
            self.assertEqual(result['parquet_members'][-1]['rows'], 0)
            self.assertEqual(result['feature_transition_fanout']['zero_transition_features'], 1)
            self.assertEqual(result['feature_transition_fanout']['mean'], 1)
            self.assertEqual(list(Path(directory).iterdir()), [path])
            self.assertEqual(before, digest(path))

    def test_cli_schema_provenance_and_input_protection(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            path = root / 'test.osw'
            self.database(path).close()
            resources = root / 'resources.json'
            resources.write_text(json.dumps({'resource_schema_version': 1, 'cgroup_final': None}))
            rescore = root / 'rescore.json'
            rescore.write_text(json.dumps({'elapsed_seconds': 0}))
            output = root / 'measure_outputs.json'
            scratch = root / 'scratch'
            scratch.mkdir()
            command = [str(path), '--dataset', 'TimSim', '--job', str(root),
                       '--temp-dir', str(scratch),
                       '--resources-json', str(resources), '--rescore-json', str(rescore),
                       '--output-json', str(output)]
            with patch('measure_outputs.tempfile.TemporaryDirectory', wraps=tempfile.TemporaryDirectory) as temporary:
                self.assertEqual(main(command), 0)
                self.assertEqual(temporary.call_args.kwargs['dir'], scratch)
            result = json.loads(output.read_text())
            self.assertEqual(result['schema_version'], 1)
            self.assertEqual(result['provenance'], {'dataset': 'TimSim', 'job': str(root)})
            self.assertEqual(result['temp_directory'], str(scratch))
            self.assertEqual(list(scratch.iterdir()), [])
            self.assertEqual(result['standalone_rescore']['elapsed_seconds'], 0)
            self.assertEqual(result['stages']['status'], 'unavailable')
            self.assertIsNone(result['resources']['cgroup_final'])
            with self.assertRaises(SystemExit):
                main(command[:-1] + [str(path)])


if __name__ == '__main__':
    unittest.main()
