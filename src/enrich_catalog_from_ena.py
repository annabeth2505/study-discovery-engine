"""
enrich_catalog_from_ena.py
One ENA fetch per study to populate deposit-derived columns into yes_catalog.tsv.

Adds (matching Sam's GMTOL schema where possible):
  - library_strategy_mix   e.g. "WGS:26; AMPLICON:151"  (per-sample strategy counts)
  - n_host_species         distinct host_scientific_name count
  - host_species_list      "; "-joined distinct hosts
  - body_sites             distinct host_body_site values
  - n_wgs_samples          count of WGS/METAGENOMIC-strategy samples (the "real shotgun" number)
Also caches audit-supporting evidence so the audit step needn't re-fetch:
  - biome_labels           distinct sample scientific_name
  - isolation_sources      distinct isolation_source
  - sample_titles_cache    a few sample titles

Checkpointed + resumable. Writes back into yes_catalog.tsv.
"""

import os
import io
import json
import time
import requests
import pandas as pd
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

ENA_URL = "https://www.ebi.ac.uk/ena/portal/api/search"
CATALOG = "../results/yes_catalog.tsv"
CHECKPOINT = "../results/ena_enrichment.json"


def make_session():
    s = requests.Session()
    retry = Retry(total=5, backoff_factor=0.5,
                  status_forcelist=[429, 500, 502, 503, 504], allowed_methods=["GET"])
    s.mount("https://", HTTPAdapter(max_retries=retry))
    return s


SESSION = make_session()


def fetch_runs(acc):
    params = {
        "result": "read_run",
        "query": f'study_accession="{acc}"',
        "fields": ("run_accession,sample_accession,library_strategy,library_source,"
                   "host_scientific_name,host_body_site,scientific_name,"
                   "isolation_source,sample_title,country"),
        "format": "tsv", "limit": 0,
    }
    try:
        r = SESSION.get(ENA_URL, params=params, timeout=90)
        if r.status_code == 200 and r.text.strip():
            return pd.read_csv(io.StringIO(r.text), sep="\t")
    except Exception:
        pass
    return pd.DataFrame()


def summarize(df):
    if df.empty:
        return None

    # library strategy mix, per SAMPLE (dedupe to sample level first)
    samp = df.drop_duplicates("sample_accession") if "sample_accession" in df.columns else df
    strat = samp["library_strategy"].dropna().value_counts() if "library_strategy" in samp else pd.Series(dtype=int)
    mix = "; ".join(f"{k}:{v}" for k, v in strat.items()) if len(strat) else ""
    n_wgs = int(sum(v for k, v in strat.items() if k in ("WGS", "METAGENOMIC")))

    def distinct(col, n=50):
        if col not in df.columns:
            return []
        return sorted(set(str(x) for x in df[col].dropna() if str(x).strip()))

    hosts = distinct("host_scientific_name")
    biomes = distinct("scientific_name")
    sites = distinct("host_body_site")
    iso = distinct("isolation_source")
    stitles = distinct("sample_title")[:5]

    return {
        "library_strategy_mix": mix,
        "n_wgs_samples": n_wgs,
        "n_host_species": len(hosts),
        "host_species_list": "; ".join(hosts[:40]),
        "body_sites": "; ".join(sites) if sites else "",
        "biome_labels": "; ".join(biomes[:20]),
        "isolation_sources": "; ".join(iso[:20]),
        "sample_titles_cache": "; ".join(stitles),
        "n_total_samples": int(samp["sample_accession"].nunique()) if "sample_accession" in samp else len(samp),
    }


def run():
    cat = pd.read_csv(CATALOG, sep="\t")
    accs = cat["study_accession"].dropna().unique().tolist()

    done = {}
    if os.path.exists(CHECKPOINT):
        with open(CHECKPOINT) as f:
            done = json.load(f)
        print(f"Resuming — {len(done)} already enriched")

    todo = [a for a in accs if a not in done]
    print(f"Enriching {len(todo)} studies from ENA...")

    for i, acc in enumerate(todo):
        if i % 25 == 0:
            print(f"  {i}/{len(todo)}")
            with open(CHECKPOINT, "w") as f:
                json.dump(done, f)
        summ = summarize(fetch_runs(acc))
        done[acc] = summ if summ else {"library_strategy_mix": "", "n_host_species": 0,
                                       "host_species_list": "", "body_sites": "",
                                       "n_wgs_samples": 0}
        time.sleep(0.25)

    with open(CHECKPOINT, "w") as f:
        json.dump(done, f)

    # ---- write columns back into catalog ----
    new_cols = ["library_strategy_mix", "n_wgs_samples", "n_host_species",
                "host_species_list", "body_sites", "biome_labels",
                "isolation_sources", "sample_titles_cache", "n_total_samples"]
    for col in new_cols:
        cat[col] = cat["study_accession"].map(lambda a: (done.get(a) or {}).get(col))

    cat.to_csv(CATALOG, sep="\t", index=False)

    print(f"\n{'='*52}\nENRICHMENT COMPLETE\n{'='*52}")
    print(f"Studies enriched: {len(done)}")
    print(f"With library_strategy_mix: {cat['library_strategy_mix'].astype(str).str.len().gt(0).sum()}")
    print(f"With host_species_list:    {cat['host_species_list'].astype(str).str.len().gt(0).sum()}")
    print(f"With body_sites:           {cat['body_sites'].astype(str).str.len().gt(0).sum()}")
    print(f"\nMulti-host studies (n_host_species > 1): {(cat['n_host_species'] > 1).sum()}")
    print(f"  ...of which >5 hosts (likely surveys):  {(cat['n_host_species'] > 5).sum()}")
    print(f"\nStudies with 0 WGS samples (all-amplicon?): {(cat['n_wgs_samples'] == 0).sum()}")
    print(f"Saved -> {CATALOG}")
    return cat


if __name__ == "__main__":
    run()