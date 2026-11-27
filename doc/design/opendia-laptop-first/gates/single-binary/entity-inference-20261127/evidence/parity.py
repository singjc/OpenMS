#!/usr/bin/env python3
"""Independent SQL + pinned PyProphet parity gate. Never edits the upstream checkout.

Usage: /path/pyprophet/.venv/bin/python parity.py FIXTURE_DIR
Frozen absolute tolerances: p=1e-12, q=1e-10, PEP=1e-2. PEP reflects
OpenMS MultipleTesting_test's existing scipy/R-density regression tolerance,
not a tolerance fitted to this fixture. Entity selection and scores are exact.
"""
import argparse
import json
from pathlib import Path
import sqlite3
import subprocess
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pyprophet
from pyprophet.io.levels_context.osw import OSWReader
from pyprophet.stats import error_statistics, lookup_values_from_error_table
from scipy_compat import enable

COMPATIBILITY = enable()

TOLERANCES = {"p": 1e-12, "q": 1e-10, "pep": 1e-2}
PIN = "4344e2b"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fixture", type=Path)
    args = parser.parse_args()
    upstream = Path(pyprophet.__file__).resolve().parents[1]
    revision = subprocess.check_output(["git", "-C", str(upstream), "rev-parse", "HEAD"], text=True).strip()
    assert revision.startswith(PIN), (revision, PIN)
    for source in ("pyprophet/stats.py", "pyprophet/io/levels_context/osw.py"):
        pinned = subprocess.check_output(["git", "-C", str(upstream), "show", f"{PIN}:{source}"])
        assert (upstream / source).read_bytes() == pinned, f"Oracle source differs from pinned commit: {source}"
    try:
        error_statistics(np.ones(240), np.ones(240), False, False, np.arange(.1, .5, .05), pi0_method="bootstrap", compute_lfdr=True)
    except Exception as exc:
        degenerate_error = f"{type(exc).__name__}: {exc}"
    else:
        raise AssertionError("Pinned upstream unexpectedly estimated all-tie fixture")
    lib = pd.read_csv(args.fixture / "precursors.tsv", sep="\t", keep_default_na=False)
    windows = pd.read_csv(args.fixture / "windows.tsv", sep="\t")
    native = pd.read_csv(args.fixture / "native.tsv", sep="\t")
    # OSW IDs must encode both modified sequence and target/decoy class. SQL
    # operates on those independently derived IDs, never native entity mappings.
    peptide_keys = sorted(set(zip(lib.sequence, lib.decoy)))
    peptide_ids = {key: i for i, key in enumerate(peptide_keys)}
    protein_keys = sorted({(p, row.decoy) for row in lib.itertuples() for p in row.proteins.split(";")})
    protein_ids = {key: i for i, key in enumerate(protein_keys)}
    db = args.fixture / "oracle.osw"
    if db.exists():
        db.unlink()
    con = sqlite3.connect(db)
    con.executescript("""
        CREATE TABLE PEPTIDE(ID INTEGER PRIMARY KEY, MODIFIED_SEQUENCE TEXT);
        CREATE TABLE PROTEIN(ID INTEGER PRIMARY KEY, PROTEIN_ACCESSION TEXT);
        CREATE TABLE PRECURSOR(ID INTEGER PRIMARY KEY, DECOY INTEGER);
        CREATE TABLE PRECURSOR_PEPTIDE_MAPPING(PRECURSOR_ID INTEGER, PEPTIDE_ID INTEGER);
        CREATE TABLE PEPTIDE_PROTEIN_MAPPING(PEPTIDE_ID INTEGER, PROTEIN_ID INTEGER);
        CREATE TABLE FEATURE(ID INTEGER PRIMARY KEY, RUN_ID INTEGER, PRECURSOR_ID INTEGER);
        CREATE TABLE SCORE_MS2(FEATURE_ID INTEGER, SCORE REAL);
    """)
    con.executemany("INSERT INTO PEPTIDE VALUES (?,?)", [(i, key[0]) for key, i in peptide_ids.items()])
    con.executemany("INSERT INTO PROTEIN VALUES (?,?)", [(i, key[0]) for key, i in protein_ids.items()])
    mappings = set()
    for row in lib.itertuples():
        peptide_id = peptide_ids[(row.sequence, row.decoy)]
        con.execute("INSERT INTO PRECURSOR VALUES (?,?)", (row.precursor, row.decoy))
        con.execute("INSERT INTO PRECURSOR_PEPTIDE_MAPPING VALUES (?,?)", (row.precursor, peptide_id))
        for protein in row.proteins.split(";"):
            mappings.add((peptide_id, protein_ids[(protein, row.decoy)]))
    con.executemany("INSERT INTO PEPTIDE_PROTEIN_MAPPING VALUES (?,?)", sorted(mappings))
    for feature, row in enumerate(windows.itertuples()):
        con.execute("INSERT INTO FEATURE VALUES (?,?,?)", (feature, row.run, row.precursor))
        con.execute("INSERT INTO SCORE_MS2 VALUES (?,?)", (feature, row.score))
    con.commit()
    reports, errors = [], []
    for level in ("peptide", "protein"):
        for context in ("run-specific", "experiment-wide", "global"):
            # Bypass only connection management to avoid network INSTALL of DuckDB
            # sqlite_scanner. Execute actual upstream SQLite reader SQL unchanged.
            reader = object.__new__(OSWReader)
            reader.config = SimpleNamespace(context_fdr=context, level=level)
            oracle = reader._read_using_sqlite(con)
            keys = peptide_keys if level == "peptide" else protein_keys
            oracle["entity"] = oracle[level + "_id"].map(lambda i: keys[int(i)][0])
            oracle["run"] = oracle.run_id.fillna(-1).astype(int)
            batches = oracle.groupby("run", sort=True) if context == "run-specific" else [(None, oracle)]
            for run, batch in batches:
                table, pi0 = error_statistics(
                    batch.loc[batch.decoy == 0, "score"], batch.loc[batch.decoy == 1, "score"],
                    parametric=False, pfdr=False, pi0_lambda=np.arange(.1, .5, .05),
                    pi0_method="bootstrap", pi0_smooth_df=3, pi0_smooth_log_pi0=False,
                    compute_lfdr=True, lfdr_trunc=True, lfdr_monotone=True,
                    lfdr_transf="probit", lfdr_adj=1.5, lfdr_eps=1e-8)
                p, _, pep, q = lookup_values_from_error_table(batch.score.to_numpy(), table)
                expected = batch.assign(p=p, q=q, pep=pep)
                actual = native[(native.level == level) & (native.context == context)]
                if run is not None:
                    actual = actual[actual.run == run]
                index = ["entity", "decoy", "run"]
                expected = expected.set_index(index).sort_index()
                actual = actual.set_index(index).sort_index()
                assert actual.index.is_unique and expected.index.is_unique
                assert actual.index.equals(expected.index), (level, context, run, "entity/run selection differs", actual.index.difference(expected.index), expected.index.difference(actual.index))
                np.testing.assert_array_equal(actual.score, expected.score)
                report = {"level": level, "context": context, "run": None if run is None else int(run), "rows": len(batch), "pi0": float(pi0["pi0"])}
                for metric, tolerance in TOLERANCES.items():
                    delta = np.abs(actual[metric].to_numpy() - expected[metric].to_numpy())
                    report[metric + "_max_abs"] = float(np.max(delta))
                    bad = ~np.isfinite(delta) | (delta > tolerance)
                    if bad.any():
                        worst = int(np.nanargmax(delta)) if np.isfinite(delta).any() else 0
                        errors.append({**report, "metric": metric, "tolerance": tolerance, "key": str(actual.index[worst]), "native": float(actual[metric].iloc[worst]), "pyprophet": float(expected[metric].iloc[worst]), "failures": int(bad.sum())})
                reports.append(report)
    payload = {"pyprophet_revision": revision, "compatibility": COMPATIBILITY, "all_tie_upstream_error": degenerate_error, "tolerances": TOLERANCES, "reports": reports, "errors": errors, "passed": not errors}
    (args.fixture / "parity_results.json").write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))
    if errors:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
