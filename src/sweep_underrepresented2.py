"""
sweep_underrepresented_v2.py
Two-tier ENA sweep for underrepresented animal gut metagenomics studies.

Tier 1 (anchor): query ENA host_tax_id=<clade anchor taxid> for every clade.
Tier 2 (species): for high-potential clades, enumerate descendant species
                  taxids and query each. "High-potential" = gap_type in
                  {search, mixed} OR the clade returned >=1 hit at anchor level
                  (proof deposits exist in that lineage).

Query per host taxid: host_tax_id=<taxid> AND library_source="METAGENOMIC"
(microbiome deposits store the host in host_tax_id; scientific_name is the biome).

Outputs:
  results/sweep_hits_v2.tsv        - all unique studies found
  results/sweep_new_studies_v2.tsv - studies NOT already in yes_catalog
"""

import json
import time
import io
import requests
import pandas as pd
from Bio import Entrez

ENA_URL = "https://www.ebi.ac.uk/ena/portal/api/search"
VOCAB_PATH = "../results/search_vocabulary.json"
CATALOG_PATH = "../results/yes_catalog.tsv"

MAX_SPECIES_PER_CLADE = 30   # cap Tier-2 breadth
NCBI_SLEEP = 0.34            # ~3 req/sec ceiling
ENA_SLEEP = 0.2


# ----------------------------------------------------------------------
def name_to_taxid(name):
    try:
        h = Entrez.esearch(db="taxonomy", term=name, retmode="xml")
        rec = Entrez.read(h); h.close()
        time.sleep(NCBI_SLEEP)
        ids = rec.get("IdList", [])
        return int(ids[0]) if ids else None
    except Exception as e:
        print(f"  name_to_taxid failed for '{name}': {e}")
        return None


def taxid_to_name(taxid):
    try:
        h = Entrez.efetch(db="taxonomy", id=str(taxid), retmode="xml")
        rec = Entrez.read(h); h.close()
        time.sleep(NCBI_SLEEP)
        return rec[0]["ScientificName"]
    except Exception:
        return None


def get_descendant_species(clade_name, max_species=MAX_SPECIES_PER_CLADE):
    """Species-rank taxids under a clade subtree."""
    try:
        term = f"{clade_name}[subtree] AND species[rank]"
        h = Entrez.esearch(db="taxonomy", term=term, retmax=max_species)
        rec = Entrez.read(h); h.close()
        time.sleep(NCBI_SLEEP)
        return [int(x) for x in rec.get("IdList", [])]
    except Exception as e:
        print(f"  descendant lookup failed for '{clade_name}': {e}")
        return []


def ena_studies_for_host(host_taxid):
    query = f'host_tax_id={host_taxid} AND library_source="METAGENOMIC"'
    params = {
        "result": "read_run",
        "query": query,
        "fields": "study_accession,scientific_name,host_scientific_name,host_tax_id,library_source,library_strategy",
        "format": "tsv",
        "limit": 0,
    }
    try:
        r = requests.get(ENA_URL, params=params, timeout=60)
        if r.status_code != 200 or not r.text.strip():
            return pd.DataFrame()
        return pd.read_csv(io.StringIO(r.text), sep="\t")
    except Exception as e:
        print(f"  ENA query failed host_tax_id={host_taxid}: {e}")
        return pd.DataFrame()


# ----------------------------------------------------------------------
def run_sweep_v2():
    with open(VOCAB_PATH) as f:
        clades = json.load(f)["clades"]

    # ---- Tier 1: anchor taxids for every clade ----
    print("TIER 1 — resolving anchors and querying ENA...")
    anchor_rows, clade_anchor_hits = [], {}
    clade_info = []  # (clade, gap_type, anchor_name, anchor_taxid, tier1_hitcount)

    for c in clades:
        anchor_name = c["tax_anchor"]
        taxid = name_to_taxid(anchor_name)
        hitcount = 0
        if taxid:
            df = ena_studies_for_host(taxid)
            time.sleep(ENA_SLEEP)
            if len(df):
                df["_clade"] = c["clade"]
                df["_tier"] = "anchor"
                anchor_rows.append(df)
                hitcount = df["study_accession"].nunique()
        clade_info.append((c["clade"], c["gap_type"], anchor_name, taxid, hitcount))
        print(f"  {c['clade']:32s} anchor={anchor_name:18s} hits={hitcount}")

    # ---- decide which clades get Tier 2 ----
    expand = []
    for (clade, gap_type, aname, ataxid, hits) in clade_info:
        if ataxid is None:
            continue
        if gap_type in ("search", "mixed") or hits >= 1:
            expand.append((clade, aname, ataxid))
    print(f"\nTIER 2 — expanding {len(expand)} high-potential clades to species level...")

    species_rows = []
    for clade, aname, ataxid in expand:
        cname = taxid_to_name(ataxid) or aname
        species = get_descendant_species(cname)
        print(f"  {clade:32s} -> {len(species)} species")
        for stid in species:
            df = ena_studies_for_host(stid)
            time.sleep(ENA_SLEEP)
            if len(df):
                df["_clade"] = clade
                df["_tier"] = "species"
                species_rows.append(df)

    # ---- pool + dedupe ----
    all_rows = anchor_rows + species_rows
    if not all_rows:
        print("No hits at all.")
        return
    hits = pd.concat(all_rows, ignore_index=True)
    hits = hits.drop_duplicates("study_accession")
    hits.to_csv("../results/sweep_hits_v2.tsv", sep="\t", index=False)

    # ---- diff vs catalog ----
    catalog = pd.read_csv(CATALOG_PATH, sep="\t")
    existing = set(catalog["study_accession"].dropna().unique())
    found = set(hits["study_accession"].dropna().unique())
    new = found - existing
    dup = found & existing

    print(f"\n{'='*52}")
    print("TWO-TIER SWEEP RESULTS vs CATALOG")
    print(f"{'='*52}")
    print(f"Studies found:            {len(found)}")
    print(f"Already in catalog:       {len(dup)}")
    print(f"NEW (missed by pipeline): {len(new)}")

    new_df = hits[hits["study_accession"].isin(new)].copy()
    new_df.to_csv("../results/sweep_new_studies_v2.tsv", sep="\t", index=False)
    print(f"Saved {len(new_df)} new studies to sweep_new_studies_v2.tsv")

    print(f"\nNEW studies by clade:")
    print(new_df["_clade"].value_counts().to_string())
    print(f"\nNEW studies by tier:")
    print(new_df["_tier"].value_counts().to_string())
    print(f"\nTop hosts among NEW:")
    print(new_df["host_scientific_name"].value_counts().head(20).to_string())

    return hits, new_df


if __name__ == "__main__":
    run_sweep_v2()