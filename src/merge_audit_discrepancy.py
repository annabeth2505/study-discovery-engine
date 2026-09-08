"""
merge_audit_discrepancy.py
Merge the audit layer (catalog_audit.tsv) and the discrepancy layer
(discrepancies.tsv) into yes_catalog.tsv.

Both passes ran long ago and wrote their own TSVs, but neither was ever merged into the
catalog -- CLAUDE.md documents them as layers 4 and 5 and quotes their results, yet the
columns were not in the file. This closes that gap.

Collisions (the catalog's own copy always wins):
  body_site       -> merged as audit_body_site. The catalog already has the
                     body_site / body_sites / llm_body_site_signal trio; a distinct
                     name keeps the audit's verdict from silently overwriting any of them.
  host            -> merged as audit_host. Does not collide, but 'host' beside
                     host_species / host_ena_value / host_resolved would be ambiguous.
  library_source  -> SKIPPED, the catalog already has it.
  n_wgs_samples   -> SKIPPED. The audit carries the OLD buggy value; fix_wgs_counts.py
                     has since recomputed the catalog's copy, which must not be undone.

NOT_CHECKED, not blank:
  The discrepancy pass only ran on the 560 paper-linked studies -- it compares paper
  against deposit, so a deposit-only study has nothing to compare. Those 484 studies get
  discrepancy_flags = "NOT_CHECKED" rather than an empty cell, so "never checked" stays
  distinguishable from "checked, nothing found" (which reads "none"). A silent blank
  would let unchecked studies pass as clean.

Backs up first. Idempotent: re-running refreshes these columns instead of duplicating.
"""

import shutil
import pandas as pd

CATALOG = "../results/yes_catalog.tsv"
BACKUP = "../results/yes_catalog.tsv.bak4"
AUDIT = "../results/catalog_audit.tsv"
DISCREP = "../results/discrepancies.tsv"

AUDIT_RENAME = {"body_site": "audit_body_site", "host": "audit_host"}
AUDIT_SKIP = {"library_source", "n_wgs_samples"}
DISCREP_COLS = ["discrepancy_flags", "discrepancy_host", "discrepancy_host_missing",
                "discrepancy_gut", "discrepancy_multihost", "discrepancy_notes"]


def run():
    cat = pd.read_csv(CATALOG, sep="\t", low_memory=False)
    aud = pd.read_csv(AUDIT, sep="\t")
    dis = pd.read_csv(DISCREP, sep="\t")
    n_before, cols_before = len(cat), len(cat.columns)

    shutil.copy(CATALOG, BACKUP)
    print(f"backup -> {BACKUP}\n")

    # ---------- audit ----------
    keep_cols = [c for c in aud.columns
                 if c != "study_accession" and c not in AUDIT_SKIP]
    aud = aud[["study_accession"] + keep_cols].rename(columns=AUDIT_RENAME)
    print(f"audit: merging {len(keep_cols)} columns, skipping {sorted(AUDIT_SKIP)}")
    print(f"       renamed {AUDIT_RENAME}")

    # ---------- discrepancy ----------
    dis = dis[["study_accession"] + [c for c in DISCREP_COLS if c in dis.columns]]
    print(f"discrepancy: merging {len(dis.columns)-1} columns for {len(dis)} checked studies")

    incoming = list(aud.columns) + list(dis.columns)
    already = [c for c in incoming if c != "study_accession" and c in cat.columns]
    if already:
        print(f"refreshing {len(already)} column(s) already present: {already}")
        cat = cat.drop(columns=already)

    merged = (cat.merge(aud, on="study_accession", how="left")
                 .merge(dis, on="study_accession", how="left"))

    # never-checked must not look clean
    checked = set(pd.read_csv(DISCREP, sep="\t").study_accession)
    merged["discrepancy_flags"] = [
        (f if isinstance(f, str) and f.strip() else "none") if acc in checked
        else "NOT_CHECKED"
        for acc, f in zip(merged.study_accession, merged.discrepancy_flags)]

    bad = [c for c in merged.columns if c.endswith("_x") or c.endswith("_y")]
    assert not bad, f"merge created suffixed columns: {bad}"
    assert len(merged) == n_before, f"row count changed: {n_before} -> {len(merged)}"

    merged.to_csv(CATALOG, sep="\t", index=False)

    print(f"\nrows    {n_before} -> {len(merged)}  (unchanged)")
    print(f"columns {cols_before} -> {len(merged.columns)}  (+{len(merged.columns)-cols_before})")
    print("\naudit results now in the catalog:")
    print("  keep      :", merged.keep.value_counts().to_dict())
    print("  is_human  :", merged.is_human.value_counts().to_dict())
    print("  gut       :", merged.gut.value_counts().to_dict())
    print("\ndiscrepancy_flags:")
    vc = merged.discrepancy_flags.value_counts()
    print(f"  NOT_CHECKED (no paper to compare): {vc.get('NOT_CHECKED', 0)}")
    print(f"  none (checked, nothing found)    : {vc.get('none', 0)}")
    print(f"  flagged                          : {int(len(merged) - vc.get('NOT_CHECKED',0) - vc.get('none',0))}")
    print(f"\nllm_notes populated: {merged.llm_notes.notna().sum()} / {len(merged)}")
    print(f"\nSaved -> {CATALOG}")
    return merged


if __name__ == "__main__":
    run()
