"""
build_ena_field_summary_xml.py
Rebuild ena_field_summary from SAMPLE XML, so the submitter's OWN vocabulary and their
CUSTOM attributes survive.

WHY NOT THE PORTAL API
----------------------
The Portal API (returnFields) normalizes tag names and drops genuinely custom ones.
Verified on SAMD00891581 (minke whale):
    XML tag                Portal API
    env_broad_scale   ->   broad_scale_environmental_context   (renamed)
    ref_biomaterial   ->   (dropped)
    source_material_id ->  (dropped)
Sample attributes are where submitters customize, so the XML is the only faithful source.

SCOPE (Option B): SAMPLE-level attributes only. Run/experiment/study fields
(library_strategy, instrument_model, read_count, study_title...) are standardized, not
customizable, and already live in other catalog columns.

SOURCE: https://www.ebi.ac.uk/ena/browser/api/xml/<accession>  -- comma-separated
accessions are accepted, up to ~500 per request (1,000 returns HTTP 414). Measured
~2.4s per 500 samples, so the full 144,665-sample catalog runs in roughly 12 minutes.
Fetches CANNOT be deduped: samples within a study share tag NAMES but carry different
VALUES (sample_name, collection_date, isolation_source), and the roll-up needs every
distinct value.

EXCLUDED TAGS -- ENA/NCBI bookkeeping that submitters do not write:
    ENA-FIRST-PUBLIC, ENA-LAST-UPDATE, INSDC secondary accession,
    BioSampleModel, NCBI submission package, ENA-CHECKLIST
Everything else is kept, including custom tags.

GRAMMAR
    field = v1; v2 | field2 = [N distinct] | field3 = value...[truncated]
  |             separates FIELDS
  =             separates field from value(s)   (split on the FIRST " = ")
  ;             separates DISTINCT VALUES within one field
  [N distinct]  a field with more than HIGH_CARD (8) distinct values
  ...[truncated] a single value longer than MAX_VALUE_LEN (200) chars
  sentinels     "not applicable", "missing" etc. shown LITERALLY
  absent        a tag no sample uses is omitted entirely

DELIMITER DECISION (chosen rule, applied consistently)
------------------------------------------------------
ENA already uses "|" INSIDE a value as its own multi-value separator, e.g.
    env_broad_scale = fecal environment [ENVO:01001029]|digestive tract environment [...]
That is semantically identical to our ";" between distinct values, so an internal "|" is
SPLIT into separate distinct values rather than escaped or replaced. Nothing is lost and
the result reads naturally:
    env_broad_scale = fecal environment [ENVO:01001029]; digestive tract environment [...]
Any residual delimiter left inside a single value after that split is replaced:
"|"->"/", ";"->",", "="->"-" (see SAFE). This is the only place a value is altered.

CAVEAT: truncation is for DISPLAY only. Distinctness is computed on FULL values, so two
values differing only past the 200-char cutoff render identically but still count as two.

Checkpointed + resumable.
"""

import os
import re
import json
import time
import argparse
import requests
import pandas as pd
import xml.etree.ElementTree as ET
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

R = "../results/"
SAMPLES = R + "samples.tsv"
OUT = R + "ena_field_summary.tsv"
CHECKPOINT = R + "ena_field_summary_xml_checkpoint.json"

XML_URL = "https://www.ebi.ac.uk/ena/browser/api/xml/"

BATCH = 400              # 500 works; 400 leaves headroom for long accessions (URI limit)
HIGH_CARD = 8
MAX_VALUE_LEN = 200

EXCLUDE_TAGS = {
    "ena-first-public", "ena-last-update", "insdc secondary accession",
    "biosamplemodel", "ncbi submission package", "ena-checklist",
}
SAFE = {"|": "/", ";": ",", "=": "-"}


def session():
    s = requests.Session()
    s.mount("https://", HTTPAdapter(max_retries=Retry(
        total=5, backoff_factor=0.5, status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=["GET"])))
    return s


SESSION = session()


def sanitize(v):
    """Safe to embed in the grammar, and length-capped. Display only."""
    s = re.sub(r"\s+", " ", str(v).strip())
    for bad, good in SAFE.items():
        s = s.replace(bad, good)
    if len(s) > MAX_VALUE_LEN:
        s = s[:MAX_VALUE_LEN] + "…[truncated]"
    return s


def fetch_xml(chunk, depth=0):
    """Fetch a batch, BISECTING on failure to isolate the unavailable accessions.

    ENA returns 404 for the WHOLE request if any single accession in it is missing
    (typically a recent NCBI BioSample not yet mirrored into ENA). Treating that as
    "no attributes" silently discarded every other sample in the batch, so a failed
    batch is split until the bad accessions are isolated and reported by name.

    Returns (parsed, missing).
    """
    if not chunk:
        return {}, []
    try:
        r = SESSION.get(XML_URL + ",".join(chunk), timeout=300)
        if r.status_code == 200:
            parsed = parse_batch(r.text)
            # ENA also returns 200 with a PARTIAL SAMPLE_SET, silently omitting
            # accessions it does not hold. Anything requested but absent from the
            # response is unavailable -- detect it here rather than losing it.
            absent = [a for a in chunk if a not in parsed]
            if not absent:
                return parsed, []
            if len(chunk) == 1:
                return parsed, absent
            # retry just the absent ones, so a partial response cannot hide them
            retry, gone = fetch_xml(absent, depth + 1)
            parsed.update(retry)
            return parsed, gone
    except Exception:
        pass
    if len(chunk) == 1:
        return {}, list(chunk)                      # this one accession is unavailable
    mid = len(chunk) // 2
    left, ml = fetch_xml(chunk[:mid], depth + 1)
    right, mr = fetch_xml(chunk[mid:], depth + 1)
    left.update(right)
    return left, ml + mr


def parse_batch(xml_text):
    """-> {sample_accession: [(tag, value), ...]}"""
    out = {}
    root = ET.fromstring(xml_text)
    for samp in root.iter("SAMPLE"):
        acc = samp.get("accession")
        if not acc:
            pid = samp.find("IDENTIFIERS/PRIMARY_ID")
            acc = pid.text if pid is not None else None
        if not acc:
            continue
        pairs = []
        for sa in samp.iter("SAMPLE_ATTRIBUTE"):
            tag = (sa.findtext("TAG") or "").strip()
            val = (sa.findtext("VALUE") or "").strip()
            if not tag or not val:
                continue
            if tag.lower() in EXCLUDE_TAGS:
                continue
            # ENA's own multi-value separator -> our distinct values
            for piece in val.split("|"):
                piece = piece.strip()
                if piece:
                    pairs.append((tag, piece))
        out[acc] = pairs
    return out


def summarize(tags):
    """{tag: {lc_value: as_returned}} -> the one-line summary string."""
    parts = []
    for tag in sorted(tags, key=str.lower):
        seen = tags[tag]
        if not seen:
            continue
        if len(seen) > HIGH_CARD:
            parts.append(f"{tag} = [{len(seen)} distinct]")
        else:
            parts.append(f"{tag} = " + "; ".join(sanitize(v) for v in seen.values()))
    return " | ".join(parts)


def run(limit=None, studies=None, out=OUT, checkpoint=CHECKPOINT, batch=BATCH):
    s = pd.read_csv(SAMPLES, sep="\t", dtype=str,
                    usecols=["study_accession", "sample_accession"])
    if studies:
        s = s[s.study_accession.isin(studies)]
    if limit:
        keep = s.study_accession.drop_duplicates().head(limit)
        s = s[s.study_accession.isin(set(keep))]

    samp2study = dict(zip(s.sample_accession, s.study_accession))
    accs = sorted(samp2study)
    print(f"{s.study_accession.nunique()} studies, {len(accs)} samples")

    acc_tags, start, covered, missing = {}, 0, {}, []
    if os.path.exists(checkpoint):
        with open(checkpoint) as f:
            ck = json.load(f)
        acc_tags, start = ck["acc_tags"], ck["next_index"]
        covered, missing = ck.get("covered", {}), ck.get("missing", [])
        print(f"Resuming at sample {start}/{len(accs)}")

    t0 = time.time()
    for i in range(start, len(accs), batch):
        chunk = accs[i:i + batch]
        parsed, gone = fetch_xml(chunk)
        missing.extend(gone)

        for acc, pairs in parsed.items():
            study = samp2study.get(acc)
            if not study:
                continue
            covered[study] = covered.get(study, 0) + 1
            d = acc_tags.setdefault(study, {})
            for tag, val in pairs:
                d.setdefault(tag, {}).setdefault(val.lower(), val)

        if (i // batch) % 20 == 0 or i + batch >= len(accs):
            with open(checkpoint, "w") as f:
                json.dump({"acc_tags": acc_tags, "next_index": i + batch,
                           "covered": covered, "missing": missing}, f)
            done = min(i + batch, len(accs))
            rate = done - start
            print(f"  {done}/{len(accs)} samples  ({time.time()-t0:.0f}s"
                  + (f", {rate/(time.time()-t0):.0f}/s)" if time.time() > t0 else ")"))

    with open(checkpoint, "w") as f:
        json.dump({"acc_tags": acc_tags, "next_index": len(accs),
                   "covered": covered, "missing": missing}, f)

    rows = []
    for study in sorted(s.study_accession.unique()):
        tags = acc_tags.get(study, {})
        summary = summarize(tags)
        n_tot = int((s.study_accession == study).sum())
        n_cov = int(covered.get(study, 0))
        if tags:
            st = "OK" if n_cov >= n_tot else "PARTIAL"
        else:
            st = "NOT_IN_ENA_XML" if n_cov == 0 else "NO_SAMPLE_ATTRIBUTES"
        rows.append({
            "study_accession": study,
            "ena_field_summary": summary,
            "n_fields": len(tags),
            "n_samples": n_tot,
            "n_samples_with_xml": n_cov,
            "status": st})

    res = pd.DataFrame(rows)
    res.to_csv(out, sep="\t", index=False)

    ok = res[res.status == "OK"]
    print(f"\n{'='*60}\nENA FIELD SUMMARY (from sample XML)\n{'='*60}")
    print(f"Studies: {len(res)}   elapsed {time.time()-t0:.0f}s")
    for k, v in res.status.value_counts().items():
        print(f"  {k:22s} {v}")
    print(f"  samples unavailable in ENA XML (404): {len(set(missing))}")
    if len(ok):
        L = ok.ena_field_summary.str.len()
        print(f"  tags per study : min {ok.n_fields.min()} median "
              f"{int(ok.n_fields.median())} max {ok.n_fields.max()}")
        print(f"  summary length : max {L.max()} mean {L.mean():.0f} "
              f"median {L.median():.0f}")
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
        studies=[x.strip() for x in a.studies.split(",")] if a.studies else None,
        out=a.out, checkpoint=a.checkpoint, batch=a.batch)
