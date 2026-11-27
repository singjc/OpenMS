"""Descriptive target entity counts at 1% (selected context + global), plus inference warnings."""
import sys, pyarrow.parquet as pq, pandas as pd
for out in sys.argv[1:]:
    print("==", out)
    for name in ("peptides", "protein_groups"):
        try: t = pd.read_csv(f"{out}/{name}.tsv", sep="\t", na_values=["NA"])
        except FileNotFoundError: print(name, "missing"); continue
        tgt = t[~t.decoy.astype(str).str.lower().isin(["true","1"])]
        idc = tgt.columns[0]
        print(f"{name}: target entities passing={tgt[tgt['pass'].astype(str).str.lower().isin(['true','1'])][idc].nunique()} of {tgt[idc].nunique()}")
    try:
        md = pq.read_schema(f"{out}/qpx/opendia.feature.parquet").metadata
        w = md.get(b"opendia_entity_inference_warnings", b"").decode()
        print("warnings:", w if w else "(none)")
    except Exception as e: print("no qpx", e)
