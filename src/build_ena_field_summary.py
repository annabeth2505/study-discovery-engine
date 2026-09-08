"""
build_ena_field_summary.py
One column per study listing every ENA field the study actually uses, with its value(s).

Output: ../results/ena_field_summary.tsv   (study_accession, ena_field_summary, n_fields,
                                            n_samples, n_runs)
Deliberately NOT written into yes_catalog.tsv -- that file has 79 columns across many
layers and merge collisions are a known failure mode. Merge later, on purpose.

WHY THE FULL FIELD LIST IS REQUESTED EXPLICITLY
-----------------------------------------------
ENA's returnFields endpoint lists every field available for the read_run result (195 of
them). If you let the portal default, it returns only the fields it feels like, and
sentinel values ("missing", "not provided") never surface. Asking for all 175
non-plumbing fields by name is what makes a study's real annotation footprint visible.

EXCLUDED FIELDS -- two groups, named so every exclusion is auditable
--------------------------------------------------------------------
1. Checksum / URL / byte plumbing (20):
     bam_aspera, bam_bytes, bam_ftp, bam_galaxy, bam_md5,
     fastq_aspera, fastq_bytes, fastq_ftp, fastq_galaxy, fastq_md5,
     sra_aspera, sra_bytes, sra_ftp, sra_galaxy, sra_md5,
     submitted_aspera, submitted_bytes, submitted_ftp, submitted_galaxy, submitted_md5
2. File-role / format plumbing (4) -- these describe the FILES, not the samples, and
   match no checksum pattern, so they must be named:
     fastq_file_role (GENERATED_FILE), bam_file_role (ARCHIVAL_FILE),
     submitted_file_role (SUBMISSION_FILE), submitted_format (FASTQ)
   NOT excluded: submitted_read_type. It looks like plumbing but carries PAIRED/SINGLE,
   which is real read structure, so it stays.
Plus study_accession: the row key, not a finding.
Everything else -- biological, environmental, sequencing, geographic -- is kept.

GRAMMAR OF THE COLUMN (documented so it can be parsed unambiguously)
-------------------------------------------------------------------
    field = v1; v2 | field2 = v1 | field3 = [N distinct]

  |            separates FIELDS
  =            separates a field name from its value(s)   (split on the FIRST " = ")
  ;            separates DISTINCT VALUES within one field
  [N distinct] replaces the values of a high-cardinality field (> 8 distinct values, and
               always for sample_title / sample_accession / run_accession) so one study
               cannot produce a 9,000-item string
  sentinels    "missing", "not provided", "not collected" etc. are shown LITERALLY, as
               ENA returned them -- never converted to blank. A field that is real for
               some samples and sentinel for others shows both:
                   host = Gallus gallus; missing
  absent       a field NO sample uses is omitted entirely, never shown as empty
  ...[truncated]  a SINGLE value longer than MAX_VALUE_LEN (200) chars, shortened for
               display. This is separate from [N distinct]: a field can hit either rule,
               both, or neither -- long values are a length problem, not a cardinality
               one (library_construction_protocol has 2 distinct values of ~1,000 chars).

  CAVEAT: long values are truncated at 200 chars FOR DISPLAY ONLY; values differing only
  past the cutoff may appear identical, but the distinct-value count reflects the FULL
  untruncated values. Two lab protocols differing only in the sequencer named at the end
  therefore render as two identical-looking truncated strings -- that is correct, not a
  duplicate.

DELIMITER SAFETY: free text can itself contain | ; or =. Inside a VALUE those three
characters are replaced with / , and - respectively (see SAFE). This is the one place a
value is altered, and it is done so the grammar above is unambiguous. sample_title, the
worst offender, is high-cardinality and never has its content emitted at all.

ROLL-UP: distinct is case-insensitive and whitespace-trimmed, but the value is displayed
as ENA returned it (first spelling wins). A field appears if ANY sample uses it.

Checkpointed + resumable; the raw per-batch response is rolled up and discarded, so
memory stays flat across the full catalog.
"""

import io
import os
import re
import json
import time
import argparse
import requests
import pandas as pd
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

R = "../results/"
CATALOG = R + "yes_catalog.tsv"
OUT = R + "ena_field_summary.tsv"
CHECKPOINT = R + "ena_field_summary_checkpoint.json"

RETURN_FIELDS = "https://www.ebi.ac.uk/ena/portal/api/returnFields"
ENA_SEARCH = "https://www.ebi.ac.uk/ena/portal/api/search"

BATCH = 25                       # studies per request; keeps payloads ~5MB
HIGH_CARD = 8                    # > this many distinct values -> [N distinct]
ALWAYS_COUNT = {"sample_title", "sample_accession", "run_accession",
                "experiment_accession", "run_alias", "sample_alias",
                "experiment_alias", "secondary_sample_accession"}
MAX_VALUE_LEN = 200              # a SINGLE value longer than this is truncated

# --- exclusions, group 1: checksum / URL / byte plumbing (20 fields) ---
# bam|fastq|sra|submitted  x  _ftp _md5 _bytes _aspera _galaxy
CHECKSUM_PLUMBING = re.compile(r"(ftp|md5|bytes|_aspera|_galaxy)$")

# --- exclusions, group 2: file-role / format plumbing ---
# These describe the FILES, not the samples, and match no checksum pattern, so they
# have to be named. Confirmed from validation data:
#   fastq_file_role      GENERATED_FILE,GENERATED_FILE
#   bam_file_role        ARCHIVAL_FILE
#   submitted_file_role  SUBMISSION_FILE,SUBMISSION_FILE
#   submitted_format     FASTQ,FASTQ
# submitted_read_type is NOT excluded: it carries PAIRED / SINGLE, which is real read
# structure, not plumbing.
FILE_ROLE_PLUMBING = {"fastq_file_role", "bam_file_role", "submitted_file_role",
                      "submitted_format"}

SAFE = {"|": "/", ";": ",", "=": "-"}


def session():
    s = requests.Session()
    s.mount("https://", HTTPAdapter(max_retries=Retry(
        total=5, backoff_factor=0.5, status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=["GET", "POST"])))
    return s


SESSION = session()


def get_fields():
    """Every read_run field ENA offers, minus file plumbing and the row key."""
    r = SESSION.get(RETURN_FIELDS, params={"result": "read_run", "format": "json"},
                    timeout=60)
    r.raise_for_status()
    names = [x["columnId"] for x in r.json()]
    checksum = [n for n in names if CHECKSUM_PLUMBING.search(n)]
    roles = [n for n in names if n in FILE_ROLE_PLUMBING]
    keep = [n for n in names
            if n not in set(checksum) | set(roles) and n != "study_accession"]
    # print the exclusions so every dropped field is a visible decision, not a silent filter
    print(f"  excluded, checksum/URL/byte plumbing ({len(checksum)}): "
          f"{', '.join(sorted(checksum))}")
    print(f"  excluded, file-role/format plumbing ({len(roles)}): "
          f"{', '.join(sorted(roles))}")
    print("  KEPT despite looking like plumbing: submitted_read_type (carries PAIRED/SINGLE)")
    return keep


def sanitize(v):
    """Make a value safe to embed in the grammar above, and cap its length.

    Truncation is for DISPLAY ONLY. Distinctness is computed on the full untruncated
    value by the caller, so two protocols differing only past the cutoff stay two
    distinct values even though they render identically.
    """
    s = str(v).strip()
    for bad, good in SAFE.items():
        s = s.replace(bad, good)
    s = re.sub(r"\s+", " ", s)
    if len(s) > MAX_VALUE_LEN:
        s = s[:MAX_VALUE_LEN] + "\u2026[truncated]"
    return s


def fetch_batch(accs, fields):
    q = " OR ".join(f'study_accession="{a}"' for a in accs)
    r = SESSION.post(ENA_SEARCH, data={
        "result": "read_run", "query": q,
        "fields": "study_accession," + ",".join(fields),
        "format": "tsv", "limit": 0}, timeout=600)
    if r.status_code != 200:
        raise RuntimeError(f"ENA {r.status_code}: {r.text[:200]}")
    if not r.text.strip():
        return pd.DataFrame()
    return pd.read_csv(io.StringIO(r.text), sep="\t", dtype=str, keep_default_na=False)


def summarize_study(df, fields):
    """Roll a study's run rows into the one-line field summary."""
    parts = []
    for f in fields:
        if f not in df.columns:
            continue
        seen = {}                                   # lowercased -> as-returned
        for v in df[f]:
            if v is None:
                continue
            s = str(v).strip()
            if not s or s.lower() == "nan":
                continue                            # truly absent, not a sentinel
            seen.setdefault(s.lower(), s)
        if not seen:
            continue                                # field unused -> omit entirely
        if f in ALWAYS_COUNT or len(seen) > HIGH_CARD:
            parts.append(f"{f} = [{len(seen)} distinct]")
        else:
            vals = "; ".join(sanitize(v) for v in seen.values())
            parts.append(f"{f} = {vals}")
    return " | ".join(parts)


def run(limit=None, studies=None, out=OUT, checkpoint=CHECKPOINT, batch=BATCH):
    fields = get_fields()
    print(f"ENA read_run fields: {len(fields)} requested "
          f"(plumbing + study_accession excluded)")

    cat = pd.read_csv(CATALOG, sep="\t", low_memory=False)
    accs = cat.study_accession.dropna().astype(str).unique().tolist()
    if studies:
        accs = [a for a in accs if a in set(studies)]
    if limit:
        accs = accs[:limit]

    done = {}
    if os.path.exists(checkpoint):
        with open(checkpoint) as f:
            done = json.load(f)
        print(f"Resuming -- {len(done)} studies already summarized")

    todo = [a for a in accs if a not in done]
    print(f"Summarizing {len(todo)} studies in batches of {batch}...\n")
    t0 = time.time()

    for i in range(0, len(todo), batch):
        chunk = todo[i:i + batch]
        t = time.time()
        try:
            df = fetch_batch(chunk, fields)
        except Exception as e:
            print(f"  batch {i//batch+1} failed ({e}) -- retrying one by one")
            df = pd.concat([fetch_batch([a], fields) for a in chunk], ignore_index=True)

        for acc in chunk:
            sub = df[df.study_accession == acc] if not df.empty else pd.DataFrame()
            if sub.empty:
                done[acc] = {"ena_field_summary": "", "n_fields": 0,
                             "n_samples": 0, "n_runs": 0, "status": "NO_ENA_ROWS"}
                continue
            summary = summarize_study(sub, fields)
            done[acc] = {
                "ena_field_summary": summary,
                "n_fields": summary.count(" | ") + 1 if summary else 0,
                "n_samples": int(sub.sample_accession.nunique())
                             if "sample_accession" in sub else 0,
                "n_runs": len(sub),
                "status": "OK"}

        with open(checkpoint, "w") as f:
            json.dump(done, f)
        print(f"  batch {i//batch+1}: {len(chunk)} studies, {len(df)} runs "
              f"({time.time()-t:.1f}s)")

    res = pd.DataFrame([{"study_accession": a, **done[a]} for a in accs if a in done])
    res = res[["study_accession", "ena_field_summary", "n_fields",
               "n_samples", "n_runs", "status"]]
    res.to_csv(out, sep="\t", index=False)

    print(f"\n{'='*60}\nENA FIELD SUMMARY COMPLETE\n{'='*60}")
    print(f"Studies: {len(res)}   elapsed {time.time()-t0:.0f}s")
    print(f"  no ENA rows: {int((res.status == 'NO_ENA_ROWS').sum())}")
    print(f"  fields used per study: min {res.n_fields.min()}  "
          f"median {int(res.n_fields.median())}  max {res.n_fields.max()}")
    print(f"\nSaved -> {out}")
    return res


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--limit", type=int, default=None)
    p.add_argument("--studies", type=str, default=None)
    p.add_argument("--out", type=str, default=OUT)
    p.add_argument("--checkpoint", type=str, default=CHECKPOINT)
    p.add_argument("--batch", type=int, default=BATCH)
    a = p.parse_args()
    run(limit=a.limit,
        studies=[s.strip() for s in a.studies.split(",")] if a.studies else None,
        out=a.out, checkpoint=a.checkpoint, batch=a.batch)
