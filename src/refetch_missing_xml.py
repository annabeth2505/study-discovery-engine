"""
refetch_missing_xml.py
Bounded recovery pass for the studies whose sample XML came back missing or partial.

WHY A SEPARATE, BOUNDED PASS
----------------------------
The main XML build lost data for 40 studies because ENA sometimes returns HTTP 200 with
a PARTIAL SAMPLE_SET -- silently omitting accessions it does not hold instead of 404ing.
The general fix (bisect a failed batch to isolate the absent accessions) is efficient
when failures are rare, but these 40 studies are exactly where failures are dense: ~64%
of their samples are genuinely absent from ENA (recent NCBI BioSamples not yet mirrored).
Bisection then costs ~2 requests per accession purely to PROVE absence -- roughly 10,000
requests for 5,049 samples, about 8 hours.

The bound: probe at most PROBE_LIMIT (40) samples per study, individually. If none of
those 40 are available, record the study as NOT_IN_ENA_XML and move on rather than
proving all 493 of its samples are missing. Recovery of real data is unaffected -- a
study with any available samples will almost certainly show it within the first 40 --
while the cost drops to at most 40 studies x 40 requests = 1,600 requests.

The trade-off is recorded honestly per study:
    n_probed          how many samples were actually tested
    n_samples         how many the study has in total
    probe_exhaustive  True when every sample was tested, False when the cap stopped it
So a study marked NOT_IN_ENA_XML off a 40-sample probe is never presented as if all of
its samples had been checked.

Output: ../results/ena_refetch.tsv, and the recovered tags are merged into the main XML
checkpoint so the two-column build picks them up.
"""

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
SUMMARY = R + "ena_field_summary.tsv"
XML_CHECKPOINT = R + "ena_field_summary_xml_checkpoint.json"
OUT = R + "ena_refetch.tsv"

XML_URL = "https://www.ebi.ac.uk/ena/browser/api/xml/"
PROBE_LIMIT = 40         # max samples tested per study before giving up
CHUNK = 40               # try this many at once first; fall back to individual probes

EXCLUDE_TAGS = {
    "ena-first-public", "ena-last-update", "insdc secondary accession",
    "biosamplemodel", "ncbi submission package", "ena-checklist",
}


def session():
    s = requests.Session()
    s.mount("https://", HTTPAdapter(max_retries=Retry(
        total=3, backoff_factor=0.4, status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=["GET"])))
    return s


SESSION = session()


def parse(xml_text):
    out = {}
    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError:
        return out
    for samp in root.iter("SAMPLE"):
        acc = samp.get("accession")
        if not acc:
            continue
        pairs = []
        for sa in samp.iter("SAMPLE_ATTRIBUTE"):
            tag = (sa.findtext("TAG") or "").strip()
            val = (sa.findtext("VALUE") or "").strip()
            if not tag or not val or tag.lower() in EXCLUDE_TAGS:
                continue
            for piece in val.split("|"):
                piece = piece.strip()
                if piece:
                    pairs.append((tag, piece))
        out[acc] = pairs
    return out


def get(accs):
    try:
        r = SESSION.get(XML_URL + ",".join(accs), timeout=180)
        return parse(r.text) if r.status_code == 200 else {}
    except Exception:
        return {}


def probe_study(accs):
    """Bounded probe. -> (parsed, n_probed, exhaustive)"""
    # one bulk attempt over the first CHUNK; if it yields anything, take it
    head = accs[:CHUNK]
    parsed = get(head)
    if parsed:
        return parsed, len(head), len(head) >= len(accs)

    # bulk returned nothing -- probe individually, up to the cap
    probed, found = 0, {}
    for a in accs[:PROBE_LIMIT]:
        probed += 1
        one = get([a])
        if one:
            found.update(one)
    return found, probed, probed >= len(accs)


def run(limit=None):
    s = pd.read_csv(SAMPLES, sep="\t", dtype=str,
                    usecols=["study_accession", "sample_accession"])
    summ = pd.read_csv(SUMMARY, sep="\t")
    targets = summ[summ.status != "OK"].study_accession.tolist()
    if limit:
        targets = targets[:limit]
    print(f"Bounded re-fetch of {len(targets)} studies "
          f"(probe cap {PROBE_LIMIT}/study)\n")

    with open(XML_CHECKPOINT) as f:
        ck = json.load(f)
    acc_tags = ck["acc_tags"]

    rows, t0 = [], time.time()
    for i, study in enumerate(targets, 1):
        accs = s[s.study_accession == study].sample_accession.tolist()
        parsed, probed, exhaustive = probe_study(accs)

        store = acc_tags.setdefault(study, {})
        for _, pairs in parsed.items():
            for tag, val in pairs:
                store.setdefault(tag, {}).setdefault(val.lower(), val)

        status = ("RECOVERED" if parsed else
                  ("NOT_IN_ENA_XML" if exhaustive else "NOT_IN_ENA_XML_PROBED"))
        rows.append({"study_accession": study, "n_samples": len(accs),
                     "n_probed": probed, "n_recovered": len(parsed),
                     "n_tags": len(store), "probe_exhaustive": exhaustive,
                     "status": status})
        print(f"  [{i:2d}/{len(targets)}] {study:14s} samples={len(accs):4d} "
              f"probed={probed:3d} recovered={len(parsed):3d} tags={len(store):3d} "
              f"{status}")

    with open(XML_CHECKPOINT, "w") as f:
        json.dump({**ck, "acc_tags": acc_tags}, f)

    res = pd.DataFrame(rows)
    res.to_csv(OUT, sep="\t", index=False)
    print(f"\n{'='*64}\nBOUNDED RE-FETCH COMPLETE  ({time.time()-t0:.0f}s)\n{'='*64}")
    print(res.status.value_counts().to_string())
    print(f"\nstudies with recovered samples: {int((res.n_recovered > 0).sum())}")
    print(f"samples recovered              : {int(res.n_recovered.sum())}")
    print(f"studies where the probe cap stopped us (not exhaustive): "
          f"{int((~res.probe_exhaustive).sum())}")
    print(f"\nmerged into {XML_CHECKPOINT}")
    print(f"Saved -> {OUT}")
    return res


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--limit", type=int, default=None)
    a = p.parse_args()
    run(limit=a.limit)
