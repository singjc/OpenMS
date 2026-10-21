#!/usr/bin/env python3
"""Small regression checks for strict, order-independent benchmark repeat validation."""

import copy
import csv
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from compare_repeats import POLICY, compare_jobs, compare_numeric


def write_table(path, rows, reverse=False):
    fields = list(rows[0])
    if reverse:
        fields.reverse()
        rows = list(reversed(rows))
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fields, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def job(root):
    (root / 'meta').mkdir(parents=True)
    (root / 'logs').mkdir()
    (root / 'results').mkdir()
    (root / 'meta/exit_code.txt').write_text('0\n')
    (root / 'meta/OpenDIA.sha256').write_text('a' * 64 + '  /binary\n')
    (root / 'logs/host_workload.json').write_text('{"abort_reason": null}')
    rows = [dict(run_name='run1', FullPeptideName=peptide, Sequence=peptide, Charge='2',
                 decoy='0', peak_group_rank='1', feature_id=str(i), filename='/old/input',
                 ProteinName='P1;P2', m_score='.005', m_score_peptide_global='.002',
                 m_score_protein_global='.003', d_score='2.5', Intensity='100', RT='123.45',
                 aggr_Peak_Area='10;20') for i, peptide in enumerate(('PEPTIDE', 'ANOTHER'))]
    write_table(root / 'results/OpenDIA.results.tsv', rows)
    for level in ('precursor', 'peptide', 'protein'):
        matrix = []
        for row in rows:
            record = {'ProteinName': row['Sequence']} if level == 'protein' else {'FullPeptideName': row['FullPeptideName']}
            if level == 'precursor':
                record.update(Charge='2', transition_group_id=row['feature_id'])
            record.update(run1='100', run2='')
            matrix.append(record)
        write_table(root / f'results/OpenDIA.{level}.matrix.tsv', matrix)
    return rows


class RepeatComparisonTest(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.reference = Path(self.temporary.name) / 'reference'
        self.candidate = Path(self.temporary.name) / 'candidate'
        self.rows = job(self.reference)
        job(self.candidate)

    def compare(self):
        return compare_jobs(self.reference, self.candidate, POLICY)

    def update(self, rows):
        write_table(self.candidate / 'results/OpenDIA.results.tsv', rows)

    def test_row_column_order_and_volatile_ids(self):
        for row in self.rows:
            row['feature_id'] = 'different'
            row['filename'] = '/different/input'
            row['ProteinName'] = 'P2;P1'
        write_table(self.candidate / 'results/OpenDIA.results.tsv', self.rows, reverse=True)
        path = self.candidate / 'results/OpenDIA.precursor.matrix.tsv'
        with path.open() as stream:
            matrix = list(csv.DictReader(stream, delimiter='\t'))
        write_table(path, matrix, reverse=True)
        self.assertTrue(self.compare()['passed'])

    def test_numeric_tolerance_and_large_changes(self):
        self.rows[0]['Intensity'] = '100.0005'
        self.rows[0]['d_score'] = '2.500001'
        self.update(self.rows)
        self.assertTrue(self.compare()['passed'])
        self.rows[0]['Intensity'] = '100.1'
        self.update(self.rows)
        self.assertFalse(self.compare()['passed'])

    def test_identity_changes_even_when_counts_match(self):
        self.rows[0]['FullPeptideName'] = 'DIFFERENT'
        self.update(self.rows)
        result = self.compare()['files'][0]
        self.assertEqual(result['missing_identities'], 1)
        self.assertEqual(result['added_identities'], 1)
        self.assertFalse(result['passed'])

    def test_threshold_crossing_cannot_hide_inside_numeric_tolerance(self):
        reference_rows = copy.deepcopy(self.rows)
        reference_rows[0]['m_score'] = '.01'
        write_table(self.reference / 'results/OpenDIA.results.tsv', reference_rows)
        self.rows[0]['m_score'] = '.0100000001'
        self.update(self.rows)
        result = self.compare()['files'][0]
        self.assertEqual(result['changed_cells'], 0)
        self.assertEqual(result['accepted_identity_changes']['precursor_only_1pct']['missing'], 1)
        self.assertFalse(result['passed'])

    def test_duplicate_keys_and_missing_globals_fail(self):
        self.update(self.rows + [self.rows[0]])
        self.assertFalse(self.compare()['passed'])
        self.rows[0]['m_score_peptide_global'] = ''
        self.update(self.rows)
        self.assertIn('global scores', self.compare()['files'][0]['error'])

    def test_missing_matrix_and_changed_missingness_fail(self):
        path = self.candidate / 'results/OpenDIA.peptide.matrix.tsv'
        path.unlink()
        self.assertFalse(self.compare()['passed'])
        write_table(path, [dict(FullPeptideName='PEPTIDE', run1='100', run2='0'),
                           dict(FullPeptideName='ANOTHER', run1='100', run2='')])
        self.assertFalse(self.compare()['files'][2]['passed'])

    def test_nonfinite_and_vector_changes_fail(self):
        for value in ('nan', 'inf', '-inf'):
            self.assertFalse(compare_numeric(value, value, POLICY['quantities'])[0])
        self.assertFalse(compare_numeric('1;2', '1', POLICY['quantities'])[0])
        self.assertFalse(compare_numeric('', '0', POLICY['quantities'])[0])
        self.assertTrue(compare_numeric('', '', POLICY['quantities'])[0])

    def test_binary_mismatch_and_interrupted_job_fail(self):
        (self.candidate / 'meta/OpenDIA.sha256').write_text('b' * 64 + '  /binary')
        with self.assertRaisesRegex(ValueError, 'Binary SHA'):
            self.compare()
        (self.candidate / 'meta/exit_code.txt').write_text('137')
        with self.assertRaisesRegex(ValueError, 'successfully'):
            self.compare()

    def test_fragment_ids_and_order_are_not_scientific_changes(self):
        for row in self.rows:
            row.update(aggr_Fragment_Annotation='12_b2_1;13_y3_1', aggr_Peak_Area='10;20', aggr_Peak_Apex='1;2')
        write_table(self.reference / 'results/OpenDIA.results.tsv', self.rows)
        changed = copy.deepcopy(self.rows)
        for row in changed:
            row.update(aggr_Fragment_Annotation='99_y3_1;98_b2_1', aggr_Peak_Area='20;10', aggr_Peak_Apex='2;1')
        self.update(changed)
        self.assertTrue(self.compare()['passed'])
        changed[0]['aggr_Fragment_Annotation'] = '99_y3_1;98_y3_1'
        self.update(changed)
        self.assertFalse(self.compare()['passed'])

    def test_nonfinite_exact_numeric_fields_and_missing_monitor_fail(self):
        for root in (self.reference, self.candidate):
            rows = copy.deepcopy(self.rows)
            rows[0]['mz'] = 'nan'
            write_table(root / 'results/OpenDIA.results.tsv', rows)
        self.assertFalse(self.compare()['passed'])
        (self.candidate / 'logs/host_workload.json').write_text('{}')
        with self.assertRaisesRegex(ValueError, 'monitor'):
            self.compare()

    def test_cli_failure_writes_evidence_and_returns_nonzero(self):
        self.rows[0]['d_score'] = '3.5'
        self.update(self.rows)
        output = Path(self.temporary.name) / 'comparison.json'
        process = subprocess.run([sys.executable, str(Path(__file__).with_name('compare_repeats.py')),
                                  str(self.reference), str(self.candidate), '--output', str(output)], capture_output=True)
        self.assertEqual(process.returncode, 1)
        self.assertFalse(json.loads(output.read_text())['passed'])
        original = output.read_bytes()
        again = subprocess.run([sys.executable, str(Path(__file__).with_name('compare_repeats.py')),
                                str(self.reference), str(self.candidate), '--output', str(output)], capture_output=True)
        self.assertNotEqual(again.returncode, 0)
        self.assertEqual(output.read_bytes(), original)


if __name__ == '__main__':
    unittest.main()
