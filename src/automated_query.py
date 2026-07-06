"""
sweep_underrepresented.py
Run in the study-triage env (from notebooks/ or repo root).

Pipeline:
1. Load search_vocabulary.json
2. Resolve each clade's tax_anchor + example_species to NCBI taxids (Entrez, free)
3. For each host taxid, query ENA: host_tax_id=<taxid> AND library_source=METAGENOMIC
4. Pool all study_accessions, dedupe
5. Diff against existing yes_catalog -> NEW vs DUPLICATE
6. Save results/sweep_hits.tsv and results/sweep_new_studies.tsv

Why host_tax_id (not scientific_name): microbiome deposits carry the biome
(e.g. "fish gut metagenome") as scientific_name; the host lives in host_tax_id.
"""

import json
import time
import io
import requests
import pandas as pd
from Bio import Entrez
from concurrent.futures import ThreadPoolExecutor

# --- config: Entrez must already be configured in your session ---
# from src.fetcher import configure_entrez; configure_entrez()

ENA_URL = "https://www.ebi.ac.uk/ena/portal/api/search"
VOCAB_PATH = "../results/search_vocabulary.json"
CATALOG_PATH = "../results/yes_catalog.tsv"


# ----------------------------------------------------------------------
# STEP 1 — resolve a scientific name to an NCBI taxid (real, not guessed)
# ----------------------------------------------------------------------
def name_to_taxid(name):
    """Resolve a scientific/common name to an NCBI taxid via Entrez taxonomy search."""
    try:
        handle = Entrez.esearch(db="taxonomy", term=name, retmode="xml")
        rec = Entrez.read(handle)
        handle.close()
        time.sleep(0.34)  # NCBI: <=3/sec without API key, safe with
        ids = rec.get("IdList", [])
        return int(ids[0]) if ids else None
    except Exception as e:
        print(f"  name_to_taxid failed for '{name}': {e}")
        return None


# ----------------------------------------------------------------------
# STEP 2 — query ENA for metagenomic deposits with a given host taxid
# ----------------------------------------------------------------------
def ena_studies_for_host(host_taxid):
    """Return set of study_accessions with metagenomic reads for this host taxid."""
    query = f'host_tax_id={host_taxid} AND library_source="METAGENOMIC"'
    params = {
        "result": "read_run",
        "query": query,
        "fields": "study_accession,scientific_name,host_scientific_name,host_tax_id,library_source,library_strategy",
        "format": "tsv",
        "limit": 0,  # 0 = no limit
    }
    try:
        r = requests.get(ENA_URL, params=params, timeout=60)
        if r.status_code != 200 or not r.text.strip():
            return pd.DataFrame()
        df = pd.read_csv(io.StringIO(r.text), sep="\t")
        return df
    except Exception as e:
        print(f"  ENA query failed for host_tax_id={host_taxid}: {e}")
        return pd.DataFrame()


# ----------------------------------------------------------------------
# MAIN
# ----------------------------------------------------------------------
def run_sweep():
    with open(VOCAB_PATH) as f:
        vocab = json.load(f)
    clades = vocab["clades"]

    # ---- collect all target names (anchor + example species) per clade ----
    targets = []  # (clade, gap_type, name)
    for c in clades:
        names = [c["tax_anchor"]] + c.get("example_species", [])
        for n in names:
            targets.append((c["clade"], c["gap_type"], n))
    print(f"Resolving {len(targets)} target names to taxids...")

    # ---- resolve names to taxids (dedupe names first) ----
    unique_names = sorted({t[2] for t in targets})
    name_taxid = {}
    for i, name in enumerate(unique_names):
        if i % 20 == 0:
            print(f"  resolving {i}/{len(unique_names)}")
        name_taxid[name] = name_to_taxid(name)
    resolved = {n: tid for n, tid in name_taxid.items() if tid}
    print(f"Resolved {len(resolved)}/{len(unique_names)} names to taxids")

    # ---- query ENA per unique taxid ----
    unique_taxids = sorted(set(resolved.values()))
    print(f"\nQuerying ENA for {len(unique_taxids)} unique host taxids...")

    taxid_to_studies = {}
    all_rows = []
    for i, tid in enumerate(unique_taxids):
        if i % 10 == 0:
            print(f"  ENA query {i}/{len(unique_taxids)}")
        df = ena_studies_for_host(tid)
        if len(df):
            taxid_to_studies[tid] = set(df["study_accession"].dropna().unique())
            all_rows.append(df)
        time.sleep(0.2)

    if not all_rows:
        print("No ENA hits at all. Check connectivity / query.")
        return

    hits = pd.concat(all_rows, ignore_index=True).drop_duplicates("study_accession")
    hits.to_csv("../results/sweep_hits.tsv", sep="\t", index=False)
    print(f"\nTotal unique studies found by sweep: {len(hits)}")

    # ---- map studies back to which clade(s) found them ----
    taxid_to_names = {}
    for name, tid in resolved.items():
        taxid_to_names.setdefault(tid, []).append(name)

    # ---- diff against existing catalog ----
    catalog = pd.read_csv(CATALOG_PATH, sep="\t")
    existing = set(catalog["study_accession"].dropna().unique())

    sweep_accessions = set(hits["study_accession"].dropna().unique())
    new_studies = sweep_accessions - existing
    duplicates = sweep_accessions & existing

    print(f"\n{'='*50}")
    print(f"SWEEP RESULTS vs CATALOG")
    print(f"{'='*50}")
    print(f"Studies found by sweep:  {len(sweep_accessions)}")
    print(f"Already in catalog:      {len(duplicates)}")
    print(f"NEW (missed by pipeline): {len(new_studies)}")

    new_df = hits[hits["study_accession"].isin(new_studies)].copy()
    new_df.to_csv("../results/sweep_new_studies.tsv", sep="\t", index=False)
    print(f"\nSaved {len(new_df)} new studies to sweep_new_studies.tsv")

    # ---- per-clade breakdown: how many NEW each clade surfaced ----
    print(f"\nNEW studies by host organism:")
    if "host_scientific_name" in new_df.columns:
        print(new_df["host_scientific_name"].value_counts().head(25).to_string())

    return hits, new_df


if __name__ == "__main__":
    run_sweep()