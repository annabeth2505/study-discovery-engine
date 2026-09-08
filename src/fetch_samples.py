"""
fetch_samples.py
Phase 1 of the sample-level catalog: pull ALL individual samples for every study
in yes_catalog.tsv from ENA. No LLM anywhere in here -- every row comes from a
real ENA response.

Fetch strategy: BULK, not a per-study loop. The ENA portal search endpoint takes
an OR'd query, so one POST covers a whole batch of studies (default 100).
Measured: ~100 studies / ~17k runs in a few seconds, vs ~1.6s per study looping.
If a batch fails as a whole, it falls back to per-study filereport calls for that
batch only, so one bad accession can't sink 99 good ones.

Output:
  ../results/samples.tsv             -- one row per SAMPLE (see grain note below)
  ../results/samples_fetch_status.tsv-- per-study OK / EMPTY / ERROR, so studies
                                        that returned nothing stay VISIBLE
  ../results/samples_checkpoint.json -- resume state

Grain note: ENA returns one row per RUN, and ~12% of samples have >1 run (up to 39).
Rows are deduped to one per (study_accession, sample_accession). EVERY run accession
is kept -- first_run_accession is one representative, all_run_accessions is the full
";"-joined list, n_runs is the count -- so the run->sample collapse loses nothing and
the column names can't be mistaken for "the sample's run".

Host resolution: host_scientific_name -> host -> scientific_name, first non-generic
value wins (same precedence as ena_fetcher.resolve_host_species, applied per SAMPLE
instead of per study). Raw fields are kept untouched alongside host_resolved, and
host_resolved_from records which field it came from -- UNRESOLVED when none qualify,
so a metadata gap stays visible instead of looking like a real answer.

Checkpointed + resumable: samples.tsv is appended batch by batch and the checkpoint
records which studies are done, so a rerun skips them.
"""

import os
import io
import json
import time
import argparse
import requests
import pandas as pd
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

ENA_SEARCH = "https://www.ebi.ac.uk/ena/portal/api/search"
ENA_FILEREPORT = "https://www.ebi.ac.uk/ena/portal/api/filereport"

CATALOG = "../results/yes_catalog.tsv"
OUT = "../results/samples.tsv"
STATUS = "../results/samples_fetch_status.tsv"
CHECKPOINT = "../results/samples_checkpoint.json"

BATCH = 100          # studies per bulk query
PAUSE = 0.25         # politeness pause between batches

FIELDS = ("run_accession,sample_accession,study_accession,host,host_scientific_name,"
          "host_tax_id,host_body_site,isolation_source,scientific_name,"
          "library_strategy,library_source,sample_title,country")

# raw ENA columns carried through untouched
RAW = ["host", "host_scientific_name", "host_tax_id", "host_body_site",
       "isolation_source", "scientific_name", "library_strategy",
       "library_source", "sample_title", "country"]

# final column order
COLS = (["sample_accession", "study_accession"]
        + RAW
        + ["host_resolved", "host_resolved_from",
           "n_runs", "first_run_accession", "all_run_accessions"])

# values that are present but mean nothing -- ENA returns these as literal strings
SENTINELS = {"missing", "na", "n/a", "not applicable", "none", "null",
             "unknown", "not collected", "not provided", "uncalculated"}

# a metagenome label is a biome, not a host binomial
GENERIC_HOST = {"metagenome", "gut metagenome", "human gut metagenome",
                "mouse gut metagenome", "bovine gut metagenome", "feces metagenome",
                "soil metagenome", "food metagenome", "environmental metagenome",
                "marine metagenome", "freshwater metagenome",
                "microbial mat metagenome"}


def _clean(v):
    """Return a usable string, or None for blank/sentinel values."""
    if v is None or pd.isna(v):
        return None
    s = str(v).strip()
    return None if not s or s.lower() in SENTINELS else s


def _is_generic(v):
    s = v.lower()
    # catch the long tail ("chicken gut metagenome") as well as the known list
    return s in GENERIC_HOST or "metagenome" in s or "microbiome" in s


def resolve_host(row):
    """Per-sample host: first non-generic of host_scientific_name -> host -> scientific_name."""
    for field in ("host_scientific_name", "host", "scientific_name"):
        v = _clean(row.get(field))
        if v and not _is_generic(v):
            return v, field
    return None, "UNRESOLVED"


def make_session():
    s = requests.Session()
    retry = Retry(total=5, backoff_factor=0.5,
                  status_forcelist=[429, 500, 502, 503, 504],
                  allowed_methods=["GET", "POST"])
    s.mount("https://", HTTPAdapter(max_retries=retry))
    return s


SESSION = make_session()


def _parse(text):
    if not text or not text.strip():
        return pd.DataFrame()
    return pd.read_csv(io.StringIO(text), sep="\t", dtype=str)


def fetch_batch(accs, timeout=300):
    """One bulk OR query for a batch of studies. Returns run-level DataFrame."""
    query = " OR ".join(f'study_accession="{a}"' for a in accs)
    r = SESSION.post(ENA_SEARCH,
                     data={"result": "read_run", "query": query,
                           "fields": FIELDS, "format": "tsv", "limit": 0},
                     timeout=timeout)
    if r.status_code != 200:
        raise RuntimeError(f"ENA {r.status_code}: {r.text[:200]}")
    return _parse(r.text)


def fetch_one(acc, timeout=120):
    """Per-study fallback via filereport, used only when a bulk batch fails."""
    try:
        r = SESSION.get(ENA_FILEREPORT,
                        params={"accession": acc, "result": "read_run",
                                "fields": FIELDS, "format": "tsv", "limit": 0},
                        timeout=timeout)
        if r.status_code == 200:
            return _parse(r.text)
    except Exception as e:
        print(f"    {acc}: {e}")
    return pd.DataFrame()


def to_sample_grain(df):
    """Collapse run rows -> one row per (study, sample), keeping every run accession."""
    if df.empty:
        return df
    for c in RAW + ["run_accession"]:
        if c not in df.columns:
            df[c] = pd.NA

    key = ["study_accession", "sample_accession"]
    # sort first so first_run_accession is reproducible, not ENA response order
    runs = (df.groupby(key, dropna=False)["run_accession"]
              .agg(n_runs="size",
                   all_run_accessions=lambda s: ";".join(sorted(x for x in s if pd.notna(x))))
              .reset_index())
    runs["first_run_accession"] = runs["all_run_accessions"].str.split(";").str[0]

    out = df.drop_duplicates(subset=key).drop(columns=["run_accession"]).merge(runs, on=key, how="left")

    resolved = out.apply(resolve_host, axis=1, result_type="expand")
    out["host_resolved"], out["host_resolved_from"] = resolved[0], resolved[1]

    return out[COLS]


def run(limit=None, batch_size=BATCH, out=OUT, status=STATUS, checkpoint=CHECKPOINT):
    cat = pd.read_csv(CATALOG, sep="\t", low_memory=False)
    accs = cat["study_accession"].dropna().astype(str).unique().tolist()
    if limit:
        accs = accs[:limit]

    done = {}
    if os.path.exists(checkpoint):
        with open(checkpoint) as f:
            done = json.load(f)
        print(f"Resuming -- {len(done)} studies already fetched")

    todo = [a for a in accs if a not in done]
    print(f"Fetching samples for {len(todo)} studies "
          f"({len(accs) - len(todo)} skipped) in batches of {batch_size}...")

    header_written = os.path.exists(out) and os.path.getsize(out) > 0
    t0 = time.time()

    for i in range(0, len(todo), batch_size):
        batch = todo[i:i + batch_size]
        t = time.time()
        try:
            runs = fetch_batch(batch)
        except Exception as e:
            print(f"  batch {i//batch_size + 1} bulk failed ({e}) -- per-study fallback")
            parts = [fetch_one(a) for a in batch]
            parts = [p for p in parts if not p.empty]
            runs = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()

        samples = to_sample_grain(runs)

        if not samples.empty:
            samples.to_csv(out, sep="\t", index=False,
                           mode="a" if header_written else "w",
                           header=not header_written)
            header_written = True

        # record status per study -- studies that came back empty stay visible.
        # NOTE: ENA resolves umbrella BioProjects to their CHILD study accessions, so
        # the response can carry accessions we never asked for. Intersect with the
        # batch before counting, or "empty" goes negative and the extras look requested.
        all_returned = set(samples["study_accession"]) if not samples.empty else set()
        returned = all_returned & set(batch)
        extras = all_returned - set(batch)
        counts = samples.groupby("study_accession").size().to_dict() if not samples.empty else {}
        for a in batch:
            done[a] = {"status": "OK" if a in returned else "EMPTY",
                       "n_samples": int(counts.get(a, 0))}
        for a in extras:                      # child studies pulled in by an umbrella
            done.setdefault(a, {"status": "CHILD_OF_UMBRELLA",
                                "n_samples": int(counts.get(a, 0))})

        with open(checkpoint, "w") as f:
            json.dump(done, f)

        print(f"  batch {i//batch_size + 1}: {len(batch)} studies, "
              f"{len(samples)} samples, {len(batch) - len(returned)} empty "
              f"({time.time() - t:.1f}s)")
        time.sleep(PAUSE)

    elapsed = time.time() - t0

    # ---- resolve EMPTY studies: umbrella project, or genuinely no data? ----
    empties = [a for a in accs if done.get(a, {}).get("status") == "EMPTY"]
    umbrella_map = {}
    for a in empties:
        kids = fetch_one(a)
        if not kids.empty and "study_accession" in kids.columns:
            children = sorted(set(kids["study_accession"].dropna()) - {a})
            if children:
                umbrella_map[a] = children
                done[a] = {"status": "UMBRELLA", "n_samples": 0,
                           "children": ";".join(children)}
    if umbrella_map:
        print(f"\n{len(umbrella_map)} umbrella project(s) resolved to child studies:")
        for a, kids in umbrella_map.items():
            print(f"  {a} -> {len(kids)} children: {';'.join(kids)}")
    with open(checkpoint, "w") as f:
        json.dump(done, f)

    # ---- dedupe: umbrella children can arrive twice (direct + via the parent) ----
    full = pd.read_csv(out, sep="\t", dtype=str)
    n_before = len(full)
    full = full.drop_duplicates(["study_accession", "sample_accession"], keep="first")
    if len(full) < n_before:
        full.to_csv(out, sep="\t", index=False)
        print(f"\nDeduped {n_before - len(full)} repeated (study, sample) rows.")

    # ---- status table: every requested study appears, whatever its outcome ----
    st = pd.DataFrame([{"study_accession": a,
                        "fetch_status": done[a]["status"],
                        "n_samples": done[a]["n_samples"],
                        "children": done[a].get("children", "")}
                       for a in accs if a in done])
    st.to_csv(status, sep="\t", index=False)

    n_empty = int((st["fetch_status"] == "EMPTY").sum())
    n_umb = int((st["fetch_status"] == "UMBRELLA").sum())
    print(f"\n{'='*52}\nSAMPLE FETCH COMPLETE\n{'='*52}")
    print(f"Studies requested: {len(accs)}")
    print(f"  returned samples: {len(st) - n_empty - n_umb}")
    print(f"  UMBRELLA (data sits under child accessions): {n_umb}")
    print(f"  EMPTY (genuinely no ENA runs): {n_empty}")
    print(f"Rows in {os.path.basename(out)}: {len(full)}")
    print(f"Total samples (requested studies only): {int(st['n_samples'].sum())}")
    print(f"Elapsed: {elapsed:.1f}s")
    print(f"Saved -> {out}")
    print(f"Status -> {status}")
    return st


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--limit", type=int, default=None,
                   help="only fetch the first N studies")
    p.add_argument("--batch", type=int, default=BATCH)
    p.add_argument("--test", action="store_true",
                   help="dry run on --limit studies, writing to *_test files")
    a = p.parse_args()

    if a.test:
        run(limit=a.limit or 5, batch_size=a.batch,
            out="../results/samples_test.tsv",
            status="../results/samples_fetch_status_test.tsv",
            checkpoint="../results/samples_checkpoint_test.json")
    else:
        run(limit=a.limit, batch_size=a.batch)
