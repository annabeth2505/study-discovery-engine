"""
merge_host_sources.py
Merge the host-source columns from host_sources.tsv into yes_catalog.tsv.

The catalog has grown many layers, and merge collisions are a known failure mode here
(a duplicated name becomes _x/_y and silently breaks column selection downstream). So
this script:
  1. backs the catalog up first (yes_catalog.tsv.bak)
  2. DROPS any incoming column that already exists in the catalog, rather than letting
     pandas suffix it -- biome_labels is the known overlap, and it came FROM the catalog
     in the first place, so the catalog's copy is authoritative
  3. asserts afterwards that no _x/_y column was created and the row count is unchanged

Adds (14 columns):
  has_paper, host_ena_referenced, host_ena_value, host_title_referenced,
  host_title_value, host_abstract_referenced, host_abstract_value, is_model_organism,
  model_organism_disease, human_by_biome_label, host_source_summary, host_source_note,
  sample_host_pattern, n_host_populated

Idempotent: re-running refreshes the same columns instead of duplicating them.
"""

import shutil
import pandas as pd

CATALOG = "../results/yes_catalog.tsv"
BACKUP = "../results/yes_catalog.tsv.bak"
SOURCES = "../results/host_sources.tsv"

ADD = ["has_paper", "host_ena_referenced", "host_ena_value",
       "host_title_referenced", "host_title_value",
       "host_abstract_referenced", "host_abstract_value",
       "is_model_organism", "model_organism_disease", "human_by_biome_label",
       "host_source_summary", "host_source_note",
       "sample_host_pattern", "n_host_populated"]


def run():
    cat = pd.read_csv(CATALOG, sep="\t", low_memory=False)
    src = pd.read_csv(SOURCES, sep="\t")
    n_before, cols_before = len(cat), len(cat.columns)

    shutil.copy(CATALOG, BACKUP)
    print(f"backup -> {BACKUP}")

    # re-running should refresh, not duplicate
    already = [c for c in ADD if c in cat.columns]
    if already:
        print(f"refreshing {len(already)} column(s) already present")
        cat = cat.drop(columns=already)

    dropped = sorted(set(src.columns) & set(cat.columns) - {"study_accession"})
    if dropped:
        print(f"dropping from incoming (catalog copy wins): {dropped}")

    merged = cat.merge(src[["study_accession"] + ADD], on="study_accession", how="left")

    # ---- guards ----
    bad = [c for c in merged.columns if c.endswith("_x") or c.endswith("_y")]
    assert not bad, f"merge created suffixed columns: {bad}"
    assert len(merged) == n_before, f"row count changed: {n_before} -> {len(merged)}"
    assert len(merged.columns) == cols_before - len(already) + len(ADD)

    merged.to_csv(CATALOG, sep="\t", index=False)

    print(f"\nrows   {n_before} -> {len(merged)}  (unchanged)")
    print(f"columns {cols_before} -> {len(merged.columns)}  (+{len(merged.columns)-cols_before})")
    print(f"\npopulated in the new columns:")
    for c in ADD:
        n = merged[c].notna().sum()
        print(f"  {c:26s} {n:5d}/{len(merged)}")
    print(f"\nSaved -> {CATALOG}")
    return merged


if __name__ == "__main__":
    run()
