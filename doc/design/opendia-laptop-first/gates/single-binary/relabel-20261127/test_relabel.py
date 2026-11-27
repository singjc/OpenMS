"""Synthetic check of the relabel evaluators: LABELS=old equals the frozen evaluator exactly; LABELS=corrected applies
instr (a 'P_target' human accession is not entrapment) and the RELABEL list (an entrapment peptide on it counts as human)."""
import os, sqlite3, subprocess, sys, tempfile, unittest
from pathlib import Path
import pandas as pd
HERE = Path(__file__).parent
PY = sys.executable

def library(path):
    c = sqlite3.connect(path)
    c.executescript("""create table PRECURSOR(ID int, TRAML_ID text, DECOY int, CHARGE int);
      create table PRECURSOR_PEPTIDE_MAPPING(PRECURSOR_ID int, PEPTIDE_ID int);
      create table PEPTIDE(ID int, UNMODIFIED_SEQUENCE text, MODIFIED_SEQUENCE text, DECOY int);
      create table PEPTIDE_PROTEIN_MAPPING(PEPTIDE_ID int, PROTEIN_ID int);
      create table PROTEIN(ID int, PROTEIN_ACCESSION text);""")
    # 0-59 human targets (accession HUMAN_target), 60-79 entrapment (X_p_target), 80-89 human ending in P ("AP_target"),
    # 90-99 entrapment that are exact human substrings (on the RELABEL list). Each target i has decoy 1000+i.
    rows = []
    for i in range(100):
        acc = 'H%d_target' % i if i < 60 else 'E%d_p_target' % i if i < 80 else 'H%dP_target' % i if i < 90 else 'S%d_p_target' % i
        seq = 'PEP%03dK' % i
        for d, pid in ((0, i), (1, 1000 + i)):
            c.execute('insert into PRECURSOR values(?,?,?,2)', (pid, ('DECOY_' if d else '') + 't%d' % i, d))
            c.execute('insert into PRECURSOR_PEPTIDE_MAPPING values(?,?)', (pid, pid))
            c.execute('insert into PEPTIDE values(?,?,?,?)', (pid, seq[::-1] if d else seq, seq[::-1] if d else seq, d))
            c.execute('insert into PEPTIDE_PROTEIN_MAPPING values(?,?)', (pid, pid))
            c.execute('insert into PROTEIN values(?,?)', (pid, ('DECOY_' if d else '') + acc))
    c.commit(); c.close()

def runs(root):
    for arm in ('cur',):
        for s in (0, 1):
            d = root / f'D_{arm}_s{s}'; d.mkdir(parents=True)
            rows = [dict(precursor_id=i, is_decoy=0, score=200 - i, q_value=0.001 * i) for i in range(100)]
            rows += [dict(precursor_id=1000 + i, is_decoy=1, score=50 - i, q_value=0.5) for i in range(100)]
            pd.DataFrame(rows).to_csv(d / 'run.precursors.tsv', sep='\t', index=False)

class T(unittest.TestCase):
    def setUp(self):
        self.t = Path(tempfile.mkdtemp()); library(self.t / 'lib.pqp'); runs(self.t / 'runs')
        pd.DataFrame({'seq': ['PEP%03dK' % i for i in range(90, 100)]}).to_csv(self.t / 'rl.tsv', sep='\t', index=False)

    def run_eval(self, script, labels=None):
        env = dict(os.environ); env.pop('LABELS', None)
        if labels: env.update(LABELS=labels, RELABEL=str(self.t / 'rl.tsv'))
        out = subprocess.run([PY, str(HERE / script), str(self.t / 'lib.pqp'), str(self.t / 'runs'), 'D'], env=env,
                             capture_output=True, text=True, check=True).stdout
        return [l for l in out.splitlines() if not l.startswith('labels:')]

    def test_old_equals_frozen(self):
        self.assertEqual(self.run_eval('eval_pdec.py'), self.run_eval('eval_pdec_relabel.py', 'old'))

    def test_corrected_counts(self):
        old = '\n'.join(self.run_eval('eval_pdec_relabel.py', 'old'))
        new = '\n'.join(self.run_eval('eval_pdec_relabel.py', 'corrected'))
        # old: entrapment = 20 (_p_target) + 10 (P_target via LIKE wildcard) + 10 (substring) = 40 of 100, r = 40/60
        self.assertIn('r = 0.6667', old)
        # corrected: entrapment = 20 only (P_target excluded by instr, substrings relabelled as human), r = 20/80
        self.assertIn('r = 0.2500', new)

    def test_diann_old_equals_frozen_and_corrected(self):
        rep = self.t / 'report.parquet'
        pd.DataFrame({'Modified.Sequence': ['PEP%03dK' % i for i in range(100)], 'Precursor.Charge': [2] * 100,
                      'Q.Value': [0.0001 * i for i in range(100)], 'PEP': [0.001 * i for i in range(100)]}).to_parquet(rep)
        def ev(script, labels=None):
            env = dict(os.environ); env.pop('LABELS', None)
            if labels: env.update(LABELS=labels, RELABEL=str(self.t / 'rl.tsv'))
            out = subprocess.run([PY, str(HERE / script), str(self.t / 'lib.pqp'), str(rep)], env=env,
                                 capture_output=True, text=True, check=True).stdout
            return [l for l in out.splitlines() if not l.startswith('labels:')]
        self.assertEqual(ev('eval_diann_r5.py'), ev('eval_diann_relabel.py', 'old'))
        self.assertIn('r = 0.6667', '\n'.join(ev('eval_diann_relabel.py', 'old')))
        new = '\n'.join(ev('eval_diann_relabel.py', 'corrected'))
        self.assertIn('r = 0.2500', new)
        self.assertIn('E 20,', new)

if __name__ == '__main__':
    unittest.main()
