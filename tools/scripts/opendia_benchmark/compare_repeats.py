#!/usr/bin/env python3
"""Compare completed same-binary benchmark jobs without relying on row order or local IDs."""

import argparse
import csv
import json
import math
import re
from pathlib import Path

from benchmark_metrics import deduplicate, valid_q

POLICY = {
    'scores': {'relative_tolerance': 1e-5, 'absolute_tolerance': 1e-8},
    'quantities': {'relative_tolerance': 1e-5, 'absolute_tolerance': 1e-6},
    'retention_time': {'relative_tolerance': 1e-6, 'absolute_tolerance': 1e-4},
}
# Database-generated IDs and absolute input paths are not scientific identities.
IGNORED = {'run_id', 'filename', 'feature_id', 'peptide_id', 'transition_group_id',
           'precursor_id', 'alignment_group_id', 'alignment_reference_feature_id'}
RAW_KEY = ('run_name', 'FullPeptideName', 'Charge', 'decoy', 'peak_group_rank')
QUANTITIES = {'Intensity', 'aggr_prec_Peak_Area', 'aggr_prec_Peak_Apex', 'aggr_Peak_Area', 'aggr_Peak_Apex'}
RT = {'RT', 'assay_rt', 'delta_rt', 'iRT', 'assay_iRT', 'delta_iRT', 'leftWidth', 'rightWidth', 'alignment_reference_rt'}
METADATA = {'transition_group_id', 'Sequence', 'FullPeptideName', 'Charge', 'ProteinName', 'GeneName'}
EXACT_NUMERIC = {'mz', 'EXP_IM', 'IM_leftWidth', 'IM_rightWidth'}


def normalize_fragments(row):
    """Align fragment quantities by ion annotation, not generated transition IDs/order."""
    if 'aggr_Fragment_Annotation' not in row:
        return row
    columns = ('aggr_Fragment_Annotation', 'aggr_Peak_Area', 'aggr_Peak_Apex')
    parts = [row[column].split(';') if row[column] else [] for column in columns]
    if len({len(values) for values in parts}) != 1:
        raise ValueError('Fragment annotation/quantity lengths differ')
    annotations = [re.sub(r'^\d+_', '', value) for value in parts[0]]
    if any(not value for value in annotations):
        raise ValueError('Blank semantic fragment annotation')
    # ponytail: exports lack fragment m/z, so duplicate labels use paired multisets,
    # not a claim of exact fragment-transition identity. Preserve multiplicity.
    aligned = sorted(zip(annotations, parts[1], parts[2]), key=lambda item: (item[0], float(item[1]), float(item[2])))
    result = dict(row)
    for index, column in enumerate(columns):
        result[column] = ';'.join(values[index] for values in aligned)
    return result



def read_table(path):
    with path.open(newline='', encoding='utf-8') as stream:
        reader = csv.DictReader(stream, delimiter='\t')
        fields = reader.fieldnames or []
        if not fields or len(fields) != len(set(fields)):
            raise ValueError(f'Empty or duplicate header: {path.name}')
        rows = list(reader)
    if not rows or any(None in row or any(value is None for value in row.values()) for row in rows):
        raise ValueError(f'Empty or malformed table: {path.name}')
    return fields, rows


def index_rows(rows, keys):
    result = {}
    for row in rows:
        identity = tuple(row[column] for column in keys)
        if any(not value for value in identity) or identity in result:
            raise ValueError(f'Blank or duplicate scientific identity: {identity}')
        result[identity] = row
    return result


def compare_numeric(left, right, tolerance):
    if left == '' or right == '':
        return left == right, None
    xs, ys = left.split(';'), right.split(';')
    if len(xs) != len(ys):
        return False, None
    largest = 0.0
    for x, y in zip(xs, ys):
        try:
            x, y = float(x), float(y)
        except ValueError:
            return False, None
        if not math.isfinite(x) or not math.isfinite(y):
            return False, None
        difference = abs(y - x)
        if not math.isfinite(difference):
            return False, None
        largest = max(largest, difference)
        if difference > tolerance['absolute_tolerance'] + tolerance['relative_tolerance'] * abs(x):
            return False, largest
    return True, largest


def accepted(rows, standard):
    chosen = deduplicate(rows)
    if standard:
        globals_ = ('m_score_peptide_global', 'm_score_protein_global')
        if any(not all(valid_q(row.get(c), 1) for c in globals_) for row in chosen.values()):
            raise ValueError('Standard view unavailable: missing or invalid global scores')
        chosen = {key: row for key, row in chosen.items() if all(valid_q(row[c], .01) for c in globals_)}
    return set(chosen)


def compare_table(reference, candidate, keys, policy, matrix=False):
    fields, rows = read_table(reference)
    other_fields, other_rows = read_table(candidate)
    if set(fields) != set(other_fields) or not set(keys) <= set(fields):
        raise ValueError(f'Column schema mismatch: {reference.name}')
    if not matrix:
        rows = [normalize_fragments(row) for row in rows]
        other_rows = [normalize_fragments(row) for row in other_rows]
    before, after = index_rows(rows, keys), index_rows(other_rows, keys)
    missing, added = before.keys() - after.keys(), after.keys() - before.keys()
    result = {'file': reference.name, 'reference_rows': len(rows), 'candidate_rows': len(other_rows),
              'missing_identities': len(missing), 'added_identities': len(added),
              'changed_cells': 0, 'numeric_cells_checked': 0, 'max_absolute_difference': {}, 'examples': []}
    for identity in sorted(before.keys() & after.keys()):
        for column in fields:
            if column in IGNORED:
                continue
            kind = None
            if matrix and column not in METADATA:
                kind = 'quantities'
            elif column in QUANTITIES:
                kind = 'quantities'
            elif column in RT:
                kind = 'retention_time'
            elif column.startswith('m_score') or column.endswith(('_pep', '_m_score', '_qvalue')) or column in {'pep', 'd_score'}:
                kind = 'scores'
            left, right = before[identity][column], after[identity][column]
            if kind:
                equal, difference = compare_numeric(left, right, policy[kind])
                result['numeric_cells_checked'] += 1
                if difference is not None:
                    result['max_absolute_difference'][kind] = max(result['max_absolute_difference'].get(kind, 0), difference)
            elif column in EXACT_NUMERIC:
                equal, _ = compare_numeric(left, right, {'relative_tolerance': 0, 'absolute_tolerance': 0})
            elif column in {'ProteinName', 'GeneName'}:
                equal = sorted(left.split(';')) == sorted(right.split(';'))
            else:
                equal = left == right
            if not equal:
                result['changed_cells'] += 1
                if len(result['examples']) < 10:
                    result['examples'].append({'identity': identity, 'column': column, 'reference': left, 'candidate': right})
    if not matrix:
        result['accepted_identity_changes'] = {}
        for standard, name in ((False, 'precursor_only_1pct'), (True, 'standard_1pct')):
            left, right = accepted(rows, standard), accepted(other_rows, standard)
            result['accepted_identity_changes'][name] = {'reference': len(left), 'candidate': len(right),
                                                        'missing': len(left - right), 'added': len(right - left)}
    result['passed'] = not (missing or added or result['changed_cells']) and all(
        not (view['missing'] or view['added']) for view in result.get('accepted_identity_changes', {}).values())
    return result


def completed_job(root):
    if (root / 'meta/exit_code.txt').read_text().strip() != '0':
        raise ValueError(f'Job did not complete successfully: {root}')
    workload = json.loads((root / 'logs/host_workload.json').read_text())
    if not isinstance(workload, dict) or 'abort_reason' not in workload or workload['abort_reason'] is not None:
        raise ValueError(f'Missing or aborted workload monitor: {root}')
    checksum = (root / 'meta/OpenDIA.sha256').read_text().split()
    if not checksum or not re.fullmatch(r'[0-9a-fA-F]{64}', checksum[0]):
        raise ValueError(f'Missing or invalid binary SHA-256: {root}')
    return checksum[0].lower()


def compare_jobs(reference, candidate, policy):
    if completed_job(reference) != completed_job(candidate):
        raise ValueError('Binary SHA-256 differs: this is not a same-binary repeat')
    files = [('OpenDIA.results.tsv', RAW_KEY, False),
             ('OpenDIA.precursor.matrix.tsv', ('FullPeptideName', 'Charge'), True),
             ('OpenDIA.peptide.matrix.tsv', ('FullPeptideName',), True),
             ('OpenDIA.protein.matrix.tsv', ('ProteinName',), True)]
    checks = []
    for name, keys, matrix in files:
        try:
            checks.append(compare_table(reference / 'results' / name, candidate / 'results' / name, keys, policy, matrix))
        except (OSError, ValueError, KeyError) as error:
            checks.append({'file': name, 'passed': False, 'error': str(error)})
    return {'passed': all(check['passed'] for check in checks), 'policy': policy,
            'scope': 'Output comparison only. Input/settings/model/dependency equivalence must be verified separately from the package provenance.',
            'provenance_checked': ['OpenDIA binary SHA-256', 'successful exits', 'explicit non-aborted monitor summaries'],
            'identity_rule': 'Exact stable scientific keys and accepted identity sets, independent of row/column order',
            'ignored_columns': sorted(IGNORED), 'files': checks}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('reference', type=Path)
    parser.add_argument('candidate', type=Path)
    parser.add_argument('--policy', type=Path, help='Predeclared JSON containing scores, quantities and retention_time tolerances')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Comparison output already exists, choose a new path')
    try:
        policy = json.loads(args.policy.read_text()) if args.policy else POLICY
        for kind in POLICY:
            for key in ('relative_tolerance', 'absolute_tolerance'):
                value = policy[kind][key]
                if not isinstance(value, (int, float)) or isinstance(value, bool) or not math.isfinite(value) or value < 0:
                    raise ValueError(f'Invalid {kind} {key}')
    except (OSError, ValueError, KeyError, TypeError) as error:
        parser.error(f'Invalid tolerance policy: {error}')
    try:
        result = compare_jobs(args.reference, args.candidate, policy)
    except (OSError, ValueError, KeyError) as error:
        result = {'passed': False, 'error': str(error), 'policy': policy}
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    print(f"{'PASS' if result['passed'] else 'FAIL'}: {args.output}")
    raise SystemExit(0 if result['passed'] else 1)


if __name__ == '__main__':
    main()
