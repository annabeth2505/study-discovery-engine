"""
build_ena_two_columns.py
The COMPLETE set of ENA fields each study used, split so the standard/custom distinction
is explicit:

    ena_standard_fields  -- every STANDARDIZED ENA field the study populated
                            (sample + run + experiment + study level, from the Portal API)
    ena_custom_fields    -- every SUBMITTER-SPECIFIC tag with NO standardized equivalent
                            (from the sample XML SAMPLE_ATTRIBUTES, verbatim)

Together the two columns are everything the study used, with nothing dropped and nothing
counted twice.

WHY TWO SOURCES
---------------
The Portal API normalizes tag names and drops genuinely custom ones; the sample XML keeps
the submitter's verbatim vocabulary but carries only sample-level attributes. Neither is
complete alone. Verified on SAMD00891581 (minke whale):
    env_broad_scale     XML tag, standardized twin = broad_scale_environmental_context
                        -> already in column 1, so NOT repeated in column 2
    ref_biomaterial     no standardized twin -> genuinely custom -> column 2
    source_material_id  no standardized twin -> genuinely custom -> column 2

THE PARTITION RULE
------------------
An XML tag goes in column 2 only if it has no standardized twin. "Twin" is decided by:
  1. normalizing the tag (lowercase; every run of non-alphanumerics -> "_") and testing it
     against ENA's returnFields vocabulary, then
  2. testing an explicit ALIASES table for checklist tags whose name differs from the
     portal field entirely (geo_loc_name -> country, lat_lon -> location, ...).
When neither matches, the tag is treated as CUSTOM. That is deliberate: erring toward
custom keeps a submitter tag visible, whereas erring toward standard would silently drop
it from both columns.

GRAMMAR (identical in both columns)
-----------------------------------
    field = v1; v2 || field2 = [N distinct] || field3 = value…[truncated]

  " || "        separates FIELDS. Chosen because ENA uses a bare "|" INSIDE values as its
                own multi-value separator, so a single "|" is NOT safe as a field separator.
  " = "         separates field from value(s)  (split on the FIRST occurrence)
  "; "          separates DISTINCT VALUES within one field
  [N distinct]  a field with more than HIGH_CARD (8) distinct values
  …[truncated]  a single value longer than MAX_VALUE_LEN (200) characters
  sentinels     "missing", "not applicable" etc. are shown LITERALLY, never blanked
  absent        a field no sample uses is omitted entirely

  ENA's internal "|" is SPLIT into separate distinct values (it means exactly what our
  "; " means), so no "|" survives inside a value. Residual ";" and "=" inside a value are
  replaced with "," and "-".

  CAVEAT: truncation is display-only. Distinctness is computed on FULL values, so two
  values differing only past the cutoff render identically but still count as two.

EXCLUSIONS (documented by name, never a silent filter)
------------------------------------------------------
Column 1, file plumbing:
  fastq_ftp fastq_bytes fastq_md5 fastq_aspera fastq_galaxy
  bam_ftp bam_bytes bam_md5 bam_aspera bam_galaxy
  sra_ftp sra_bytes sra_md5 sra_aspera sra_galaxy
  submitted_ftp submitted_bytes submitted_md5 submitted_aspera submitted_galaxy
  fastq_file_role bam_file_role submitted_file_role submitted_format submitted_read_type
Column 2, ENA/INSDC bookkeeping every sample carries:
  ENA-FIRST-PUBLIC  ENA-LAST-UPDATE  INSDC secondary accession
  BioSampleModel  NCBI submission package  ENA-CHECKLIST

Output: ../results/ena_field_summary.tsv  (study_accession + both columns + counts).
For a later deliberate merge into yes_catalog.tsv, not written there directly.
"""

import io
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
XML_CHECKPOINT = R + "ena_field_summary_xml_checkpoint.json"
OUT = R + "ena_field_summary.tsv"

RETURN_FIELDS = "https://www.ebi.ac.uk/ena/portal/api/returnFields"
ENA_SEARCH = "https://www.ebi.ac.uk/ena/portal/api/search"

BATCH = 25
HIGH_CARD = 8
MAX_VALUE_LEN = 200
FIELD_SEP = " || "
VALUE_SEP = "; "

CHECKSUM_PLUMBING = re.compile(r"(ftp|md5|bytes|_aspera|_galaxy)$")
FILE_ROLE_PLUMBING = {"fastq_file_role", "bam_file_role", "submitted_file_role",
                      "submitted_format", "submitted_read_type"}
ALWAYS_COUNT = {"sample_title", "sample_accession", "run_accession",
                "experiment_accession", "run_alias", "sample_alias",
                "experiment_alias", "secondary_sample_accession", "sample_name"}

BOOKKEEPING = {"ena-first-public", "ena-last-update", "insdc secondary accession",
               "biosamplemodel", "ncbi submission package", "ena-checklist"}

# Checklist tags whose standardized twin has a DIFFERENT name. Normalization alone cannot
# connect these, so they are stated explicitly. Anything not listed and not matching a
# returnField is treated as custom.
ALIASES = {
    "env_broad_scale": "broad_scale_environmental_context",
    "env_local_scale": "local_environmental_context",
    "env_medium": "environmental_medium",
    # MIxS short forms whose standardized twin is spelled out in full
    "env_biome": "environment_biome",
    "env_feature": "environment_feature",
    "env_material": "environment_material",
    "latitude": "lat",
    "longitude": "lon",
    "tissue": "tissue_type",
    "environment_biome": "environment_biome",
    "environment_feature": "environment_feature",
    "environment_material": "environment_material",
    "geo_loc_name": "country",
    "geographic_location_country_and_or_sea": "country",
    "geographic_location_latitude": "lat",
    "geographic_location_longitude": "lon",
    "geographic_location_region_and_locality": "region",
    "lat_lon": "location",
    "organism": "scientific_name",
    "common_name": "scientific_name",
    "host_scientific_name": "host_scientific_name",
    "host_common_name": "host",
    "host_subject_id": "host_subject_id",
    "sequencing_method": "sequencing_method",
    "investigation_type": "investigation_type",
    "project_name": "project_name",
    "collection_date": "collection_date",
    "isolation_source": "isolation_source",
    "host": "host",
    "sample_name": "sample_alias",
    "sample_description": "sample_description",
    "sample_title": "sample_title",
    "tissue_type": "tissue_lib",
    "host_body_site": "host_body_site",
    "host_sex": "host_sex",
    "host_age": "host_age",
    "checklist": "checklist",
}

SAFE = {";": ",", "=": "-"}


def session():
    s = requests.Session()
    s.mount("https://", HTTPAdapter(max_retries=Retry(
        total=5, backoff_factor=0.5, status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=["GET", "POST"])))
    return s


SESSION = session()


def norm(name):
    """lowercase, every run of non-alphanumerics -> single underscore."""
    return re.sub(r"_+", "_", re.sub(r"[^a-z0-9]+", "_", str(name).lower())).strip("_")


def sanitize(v):
    s = re.sub(r"\s+", " ", str(v).strip())
    for bad, good in SAFE.items():
        s = s.replace(bad, good)
    if len(s) > MAX_VALUE_LEN:
        s = s[:MAX_VALUE_LEN] + "…[truncated]"
    return s


def get_portal_fields():
    r = SESSION.get(RETURN_FIELDS, params={"result": "read_run", "format": "json"},
                    timeout=60)
    r.raise_for_status()
    names = [x["columnId"] for x in r.json()]
    checksum = [n for n in names if CHECKSUM_PLUMBING.search(n)]
    roles = [n for n in names if n in FILE_ROLE_PLUMBING]
    keep = [n for n in names
            if n not in set(checksum) | set(roles) and n != "study_accession"]
    print(f"  Portal fields: {len(names)} available, {len(keep)} requested")
    print(f"    excluded checksum/URL/byte ({len(checksum)}): {', '.join(sorted(checksum))}")
    print(f"    excluded file-role/format ({len(roles)}): {', '.join(sorted(roles))}")
    return keep, set(names)


def add_values(store, field, raw):
    """Split ENA's internal '|' into distinct values; dedupe case-insensitively."""
    for piece in str(raw).split("|"):
        piece = piece.strip()
        if piece and piece.lower() != "nan":
            store.setdefault(field, {}).setdefault(piece.lower(), piece)


def format_fields(store):
    parts = []
    for f in sorted(store, key=str.lower):
        seen = store[f]
        if not seen:
            continue
        if f in ALWAYS_COUNT or len(seen) > HIGH_CARD:
            parts.append(f"{f} = [{len(seen)} distinct]")
        else:
            parts.append(f"{f} = " + VALUE_SEP.join(sanitize(v) for v in seen.values()))
    return FIELD_SEP.join(parts)


def fetch_portal(accs, fields):
    """-> {study: {field: {lc: value}}}"""
    out = {}
    for i in range(0, len(accs), BATCH):
        chunk = accs[i:i + BATCH]
        q = " OR ".join(f'study_accession="{a}"' for a in chunk)
        r = SESSION.post(ENA_SEARCH, data={
            "result": "read_run", "query": q,
            "fields": "study_accession," + ",".join(fields),
            "format": "tsv", "limit": 0}, timeout=600)
        if r.status_code != 200 or not r.text.strip():
            continue
        df = pd.read_csv(io.StringIO(r.text), sep="\t", dtype=str, keep_default_na=False)
        for study, sub in df.groupby("study_accession"):
            store = out.setdefault(study, {})
            for f in fields:
                if f in sub.columns:
                    for v in sub[f]:
                        if v:
                            add_values(store, f, v)
        if (i // BATCH) % 10 == 0:
            print(f"    portal {min(i+BATCH, len(accs))}/{len(accs)} studies")
    return out


def run(studies=None, out=OUT):
    t0 = time.time()
    print("Loading ENA vocabulary...")
    fields, all_portal_names = get_portal_fields()
    standard_norm = {norm(n) for n in all_portal_names}
    # second-chance match ignoring underscores entirely: the tag "host taxid" normalizes
    # to host_taxid, which is not host_tax_id, yet they are plainly the same field.
    standard_tight = {norm(n).replace("_", "") for n in all_portal_names}

    cat = pd.read_csv(CATALOG, sep="\t", low_memory=False)
    accs = cat.study_accession.dropna().astype(str).unique().tolist()
    if studies:
        accs = [a for a in accs if a in set(studies)]

    print(f"\nFetching Portal data for {len(accs)} studies...")
    portal = fetch_portal(accs, fields)

    print("\nLoading sample XML attributes...")
    with open(XML_CHECKPOINT) as f:
        xml_tags = json.load(f)["acc_tags"]
    print(f"  {len(xml_tags)} studies in the XML checkpoint")

    rows, custom_seen, twin_seen = [], {}, {}
    for acc in accs:
        std = portal.get(acc, {})
        custom = {}
        for tag, vals in (xml_tags.get(acc) or {}).items():
            if tag.lower() in BOOKKEEPING:
                continue
            n = norm(tag)
            twin = (n in standard_norm
                    or n.replace("_", "") in standard_tight
                    or ALIASES.get(n) in all_portal_names)
            if twin:
                twin_seen[tag] = twin_seen.get(tag, 0) + 1
                continue                      # already represented in column 1
            custom[tag] = vals
            custom_seen[tag] = custom_seen.get(tag, 0) + 1

        rows.append({
            "study_accession": acc,
            "ena_standard_fields": format_fields(std),
            "ena_custom_fields": format_fields(custom),
            "n_standard_fields": len(std),
            "n_custom_fields": len(custom),
        })

    res = pd.DataFrame(rows)
    res.to_csv(out, sep="\t", index=False)

    print(f"\n{'='*66}\nTWO-COLUMN ENA FIELD SUMMARY\n{'='*66}")
    print(f"Studies: {len(res)}   elapsed {time.time()-t0:.0f}s")
    for col, n in (("ena_standard_fields", "n_standard_fields"),
                   ("ena_custom_fields", "n_custom_fields")):
        L = res[col].str.len()
        print(f"\n{col}:")
        print(f"  fields per study: median {int(res[n].median())}  max {int(res[n].max())}"
              f"  studies with none {int((res[n]==0).sum())}")
        print(f"  string length   : median {int(L.median())}  max {int(L.max())}")

    print(f"\nXML tags routed to column 1 (standardized twin exists): {len(twin_seen)} distinct")
    print("  " + ", ".join(sorted(twin_seen, key=lambda k: -twin_seen[k])[:18]))
    print(f"\nDISTINCT CUSTOM TAGS: {len(custom_seen)}")
    for tag, n in sorted(custom_seen.items(), key=lambda kv: -kv[1])[:45]:
        print(f"  {n:5d} studies  {tag}")
    print(f"\nSaved -> {out}")
    return res


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--studies", type=str, default=None)
    p.add_argument("--out", type=str, default=OUT)
    a = p.parse_args()
    run(studies=[x.strip() for x in a.studies.split(",")] if a.studies else None,
        out=a.out)
