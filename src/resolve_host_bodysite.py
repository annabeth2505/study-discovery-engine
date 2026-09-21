"""
resolve_host_bodysite.py
Re-resolve each study's HOST and BODY SITE from the evidence that already answers them.

Driven by the 75-study manual validation (results/error_rate_gmtol.xlsx). Nearly every
host error and body-site "uncertain" traced to the SAME gap -- and it turned out to be an
input gap, not a reasoning one: 11 of the 13 failing studies are among the 181 that
enrich_catalog_from_ena.py never ran on, so the audit judged them with biome_labels,
isolation_sources and sample_titles_cache all BLANK. samples.tsv holds that evidence for
every study; this script uses it directly, plus the ENA study title/description, which
no step had ever been given.

HOST -- first source that yields a real animal wins:
  1. ena_host_field  ENA sample host / host_scientific_name (structured), tax_id-collapsed
  2. biome_label     animal named in scientific_name: "rat gut metagenome" -> rat
  3. ena_study_text  animal named in the ENA study title/description (Claude reads it)
  4. paper           host from the paper title/abstract, reliable paper links only
  5. none            no animal named anywhere -> "unknown"
Multi-host is kept as a "; "-joined list, never collapsed to one animal.

Biome words map to a binomial ONLY where the common name is one species (human, mouse,
pig, chicken, horse, goat, sheep, rabbit). General words stay general -- "rat", "fish",
"bovine", "insect" -- per the rule that a general term is never expanded into a species.
Environmental prefixes (soil, wastewater, synthetic, food...) are never a host.

BODY SITE -- scanned across scientific_name, isolation_source, sample_title and
host_body_site for every sample:
  GUT        gut terms, no non-gut terms
  GUT+OTHER  gut AND non-gut terms (e.g. gut + blood) -- a mixed study stays visible
  NON_GUT    only non-gut terms
  UNCERTAIN  nothing either way in the sample fields OR the ENA study text
The ENA study text is only a FALLBACK for body site: titles often mention sites they did
not sample ("unlike oral studies..."), so they must not override sample evidence.
Control samples (negative/positive control, blank, mock) are ignored for body site.

New columns: ena_study_title, ena_study_description, host_final, host_final_from,
             body_site_final, body_site_final_from, body_site_final_detail
Existing host_species / audit_host / body_site / audit_body_site / gut are NOT modified,
so before and after stay comparable.

Usage:
  python resolve_host_bodysite.py --studies A,B,...   # dry run, writes nothing to catalog
  python resolve_host_bodysite.py --merge             # all studies, merged into catalog
"""

import io
import os
import re
import sys
import json
import time
import shutil
import argparse
import importlib.util
import requests
import pandas as pd
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

R = "../results/"
CATALOG = R + "yes_catalog.tsv"
SAMPLES = R + "samples.tsv"
STUDY_TEXT = R + "ena_study_text.tsv"
LLM_CHECKPOINT = R + "host_studytext_verdicts.json"
OUT = R + "host_bodysite_resolved.tsv"
BACKUP = R + "yes_catalog.tsv.bak8"
VALIDATION = R + "error_rate_gmtol.xlsx"

ENA_SEARCH = "https://www.ebi.ac.uk/ena/portal/api/search"
MODEL = "claude-haiku-4-5-20251001"

NEW_COLS = ["ena_study_title", "ena_study_description", "host_final", "host_final_from",
            "body_site_final", "body_site_final_from", "body_site_final_detail"]

# ------------------------------------------------------------------ body site

GUT_RX = re.compile(
    r"\b(gut|guts|fa?ec(?:es|al)|intestin\w*|ca?ec(?:a|al|um)|colon|colonic|rumen|"
    r"ruminal|hindgut|foregut|midgut|digestive (?:tract|system)|stool|"
    r"gastro-?intestinal|digesta|ile(?:um|al)|jejun(?:um|al)|duoden(?:um|al)|"
    r"cloaca\w*|rect(?:um|al))\b", re.I)

NONGUT_RX = re.compile(
    r"\b(oral|saliva\w*|mouth|tongue|dental|plaque|skin|vagina\w*|blood|serum|plasma|"
    r"urine|urinary|respiratory|nasal|nasopharyn\w*|lung|sputum|bronch\w*|milk|breast|"
    r"semen|gills?|liver|spleen|kidney|brain|muscle|soil|sediment|seawater|freshwater|"
    r"wastewater|water|sludge|air|dust|indoor|surface)\b", re.I)

CONTROL_RX = re.compile(r"\b(control|blank|mock|negative|positive|extraction)\b", re.I)
SAMPLE_SITE_FIELDS = ["scientific_name", "isolation_source", "sample_title",
                      "host_body_site"]

# ------------------------------------------------------------------ biome -> host

BIOME_HOST = {   # binomial only where the common name is one species
    "human": "Homo sapiens", "mouse": "Mus musculus", "pig": "Sus scrofa",
    "chicken": "Gallus gallus", "horse": "Equus caballus", "goat": "Capra hircus",
    "sheep": "Ovis aries", "rabbit": "Oryctolagus cuniculus",
    # general words stay general
    "rat": "rat", "fish": "fish", "insect": "insect", "bovine": "bovine", "bee": "bee",
    "snake": "snake", "termite": "termite", "bird": "bird", "invertebrate": "invertebrate",
    "primate": "primate", "bat": "bat", "canine": "canine", "cetacean": "cetacean",
    "feline": "feline", "shrimp": "shrimp", "spider": "spider",
}
BIOME_SITE_WORDS = re.compile(
    r"\b(gut|feces|faeces|fecal|faecal|oral|skin|vaginal|intestinal|stool|rumen|blood|"
    r"milk|urine|nasal|lung|respiratory|upper respiratory tract)\b", re.I)

SENTINELS = {"missing", "na", "n/a", "not applicable", "none", "null", "unknown",
             "not collected", "not provided", "nan", ""}


def _clean(v):
    if v is None or (isinstance(v, float) and pd.isna(v)):
        return None
    s = str(v).strip()
    return None if s.lower() in SENTINELS else s


def biome_animals(labels):
    """'rat gut metagenome' -> {'rat'}; 'feces metagenome' / 'soil metagenome' -> set()."""
    found = []
    for lab in labels:
        for part in str(lab).split(";"):
            p = part.strip().lower()
            if not re.search(r"metagenomes?$", p):
                continue
            w = re.sub(r"metagenomes?$", "", p)
            w = re.sub(r"\s+", " ", BIOME_SITE_WORDS.sub("", w)).strip()
            if w in BIOME_HOST and BIOME_HOST[w] not in found:
                found.append(BIOME_HOST[w])
    return found


def body_site(sub, study_text):
    gut, other = set(), set()
    for f in SAMPLE_SITE_FIELDS:
        if f not in sub.columns:
            continue
        for v in sub[f].dropna().astype(str).unique():
            if CONTROL_RX.search(v):
                continue
            gut |= {m.lower() for m in GUT_RX.findall(v)}
            other |= {m.lower() for m in NONGUT_RX.findall(v)}
    src = "sample_fields"
    if not gut and not other and study_text:          # fallback only
        gut = {m.lower() for m in GUT_RX.findall(study_text)}
        src = "ena_study_text" if gut else "none"
    if not gut and not other:
        return "UNCERTAIN", "none", ""
    if gut and other:
        label = "GUT+OTHER"
    elif gut:
        label = "GUT"
    else:
        label = "NON_GUT"
    detail = "; ".join(filter(None, [
        ("gut: " + ", ".join(sorted(gut))) if gut else "",
        ("other: " + ", ".join(sorted(other))) if other else ""]))
    return label, src, detail


# ------------------------------------------------------------------ structured host

def load_pattern_fn():
    spec = importlib.util.spec_from_file_location(
        "psp", os.path.join(os.path.dirname(__file__), "parse_sample_patterns.py"))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m.study_pattern


def structured_hosts(sub, study_pattern):
    """Only the ENA host fields -- never the scientific_name fallback, which can hold
    bacteria or biome labels. tax_id-collapsed so RAT / Mus musculus variants merge."""
    s = sub[sub.host_resolved_from.isin(["host_scientific_name", "host"])].copy()
    if s.empty:
        return []
    # submitters sometimes put per-animal IDs in the host field ("Gorilla1".."Gorilla30",
    # "calves 1".."calves 44"). Only a 1-2 word name directly followed by a number is
    # stripped, so strain names like "Mus musculus C57BL/6" are left alone.
    s["host_resolved"] = s.host_resolved.astype(str).str.replace(
        r"^([A-Za-z]+(?: [A-Za-z]+)?)[\s_\-]*\d+$", r"\1", regex=True)
    lst = study_pattern(s)["distinct_host_list"]
    return [h for h in lst.split("; ") if h] if lst else []


# ------------------------------------------------------------------ ENA study text

def session():
    s = requests.Session()
    s.mount("https://", HTTPAdapter(max_retries=Retry(
        total=5, backoff_factor=0.5, status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=["GET", "POST"])))
    return s


def fetch_study_text(accs):
    """Bulk fetch study_title / study_description; cached so reruns never refetch."""
    have = pd.read_csv(STUDY_TEXT, sep="\t", dtype=str) if os.path.exists(STUDY_TEXT) \
        else pd.DataFrame(columns=["study_accession", "study_title", "study_description"])
    todo = sorted(set(accs) - set(have.study_accession))
    if todo:
        sess, rows = session(), []
        for i in range(0, len(todo), 100):
            batch = todo[i:i + 100]
            q = " OR ".join(f'study_accession="{a}"' for a in batch)
            r = sess.post(ENA_SEARCH, data={
                "result": "study", "query": q,
                "fields": "study_accession,study_title,study_description",
                "format": "tsv", "limit": 0}, timeout=180)
            if r.status_code == 200 and r.text.strip():
                df = pd.read_csv(io.StringIO(r.text), sep="\t", dtype=str)
                # an umbrella project also returns its children -- keep what was asked
                rows.append(df[df.study_accession.isin(batch)])
        if rows:
            have = pd.concat([have] + rows, ignore_index=True).drop_duplicates(
                "study_accession")
            have.to_csv(STUDY_TEXT, sep="\t", index=False)
        print(f"  ENA study text: fetched {len(todo)}, cached {len(have)}")
    return have.set_index("study_accession")


# ------------------------------------------------------------------ Claude, step 3

SYSTEM = (
    "You read an ENA sequencing project's TITLE and DESCRIPTION and name the host "
    "ANIMAL the sequenced samples came from. Rules: "
    "(1) The host is the animal the SAMPLES came from, not every organism mentioned. "
    "Microbes, bacteria, viruses and phages are NEVER the host. "
    "(2) MODEL ORGANISMS: if the samples came from an animal modelling a human disease "
    "(mice, rats, germ-free, gnotobiotic), the host is that animal, not human. "
    "(3) Prefer the Latin binomial. If the inputs give one for the host (including the "
    "SAMPLE ORGANISM field), return it. If a common name denotes exactly ONE species, "
    "return that species ('Atlantic cod' -> Gadus morhua, 'honey bee' -> Apis mellifera). "
    "If the name is general ('pigs', 'fish', 'ruminants', 'slug'), keep it general -- "
    "never invent a species from a general term. "
    "(4) Several animals sampled -> list them all. "
    "(5) If no host animal is named or clearly implied, return an empty list. "
    "Return ONLY JSON.")

PROMPT = """Project {acc}

TITLE: {title}

DESCRIPTION: {desc}

SAMPLE ORGANISM field (ENA scientific_name, non-metagenome values; may be the host
animal, or a microbe -- judge which): {organisms}

Return ONLY:
{{"hosts": ["..."], "is_model_organism": false, "note": "one short sentence on where the samples came from"}}"""


def call_llm(client, prompt, retries=6):
    import anthropic
    for a in range(retries):
        try:
            m = client.messages.create(model=MODEL, max_tokens=400, system=SYSTEM,
                                       messages=[{"role": "user", "content": prompt}])
            t = m.content[0].text.strip().replace("```json", "").replace("```", "")
            s, e = t.find("{"), t.rfind("}")
            return json.loads(t[s:e + 1])
        except anthropic.RateLimitError:
            time.sleep(13 * (a + 1))
        except (json.JSONDecodeError, ValueError):
            return {"hosts": None, "note": "PARSE_ERROR"}
        except Exception as ex:
            return {"hosts": None, "note": f"ERROR: {str(ex)[:100]}"}
    return {"hosts": None, "note": "ERROR: max retries"}


def llm_hosts(need, text, organisms):
    """need: list of accessions still unresolved after steps 1-2."""
    from pathlib import Path
    from dotenv import load_dotenv
    import anthropic
    load_dotenv(Path(__file__).resolve().parent.parent / ".env")
    client = anthropic.Anthropic(api_key=os.getenv("ANTHROPIC_API_KEY"))

    done = json.load(open(LLM_CHECKPOINT)) if os.path.exists(LLM_CHECKPOINT) else {}
    todo = [a for a in need if a not in done and a in text.index
            and (_clean(text.loc[a, "study_title"]) or _clean(text.loc[a, "study_description"]))]
    print(f"  Claude on ENA study text: {len(todo)} to call ({len(done)} cached)")
    for i, a in enumerate(todo):
        done[a] = call_llm(client, PROMPT.format(
            acc=a, title=_clean(text.loc[a, "study_title"]) or "(none)",
            desc=(_clean(text.loc[a, "study_description"]) or "(none)")[:3000],
            organisms=organisms.get(a) or "(none)"))
        if i % 10 == 0:
            json.dump(done, open(LLM_CHECKPOINT, "w"))
        time.sleep(0.3)
    json.dump(done, open(LLM_CHECKPOINT, "w"))
    return done


# ------------------------------------------------------------------ main

def resolve(accs):
    cat = pd.read_csv(CATALOG, sep="\t", low_memory=False).set_index("study_accession")
    s = pd.read_csv(SAMPLES, sep="\t", dtype=str)
    s = s[s.study_accession.isin(accs)]
    groups = dict(tuple(s.groupby("study_accession")))
    study_pattern = load_pattern_fn()
    text = fetch_study_text(accs)

    rows, need = {}, []
    for a in accs:
        sub = groups.get(a, pd.DataFrame(columns=s.columns))
        title = _clean(text.loc[a, "study_title"]) if a in text.index else None
        desc = _clean(text.loc[a, "study_description"]) if a in text.index else None
        site, site_from, detail = body_site(sub, " ".join(filter(None, [title, desc])))

        hosts, frm = structured_hosts(sub, study_pattern), "ena_host_field"
        if not hosts:
            hosts, frm = biome_animals(sub.scientific_name.dropna().unique()), "biome_label"
        if not hosts:
            need.append(a)
            frm = None
        rows[a] = {"study_accession": a, "ena_study_title": title or "",
                   "ena_study_description": desc or "",
                   "host_final": "; ".join(hosts), "host_final_from": frm,
                   "body_site_final": site, "body_site_final_from": site_from,
                   "body_site_final_detail": detail}

    organisms = {a: "; ".join(sorted({v for v in groups[a].scientific_name.dropna()
                                      if "metagenome" not in v.lower()})[:10])
                 for a in need if a in groups}
    verdicts = llm_hosts(need, text, organisms) if need else {}
    for a in need:
        v = verdicts.get(a) or {}
        # same vocabulary as the biome step: 'pig'/'pigs' -> Sus scrofa, 'rat' stays rat
        hosts = []
        for h in (v.get("hosts") or []):
            h = _clean(h)
            if h:
                key = h.lower().rstrip("s") if h.lower().rstrip("s") in BIOME_HOST else h.lower()
                hosts.append(BIOME_HOST.get(key, h))
        if hosts:
            rows[a].update(host_final="; ".join(hosts), host_final_from="ena_study_text")
            continue
        r = cat.loc[a] if a in cat.index else None
        paper = None
        if r is not None and bool(r.get("paper_host_reliable", True)):
            paper = _clean(r.get("host_abstract_value")) or _clean(r.get("host_title_value"))
        if paper:
            rows[a].update(host_final=paper, host_final_from="paper")
        else:
            note = "NOT_CHECKED" if str(v.get("note", "")).startswith(("PARSE", "ERROR")) else ""
            rows[a].update(host_final="unknown", host_final_from=note or "none")
    return pd.DataFrame(rows.values())[["study_accession"] + NEW_COLS]


def merge(res):
    cat = pd.read_csv(CATALOG, sep="\t", low_memory=False)
    n, cols = len(cat), len(cat.columns)
    shutil.copy(CATALOG, BACKUP)
    cat = cat.drop(columns=[c for c in NEW_COLS if c in cat.columns])
    clash = sorted(set(res.columns) & set(cat.columns) - {"study_accession"})
    assert not clash, f"collision: {clash}"
    m = cat.merge(res, on="study_accession", how="left")
    assert len(m) == n and not [c for c in m.columns if c.endswith(("_x", "_y"))]
    m.to_csv(CATALOG, sep="\t", index=False)
    print(f"\nmerged: backup {BACKUP}; columns {cols} -> {len(m.columns)}; rows {n}")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--studies", type=str, default=None)
    p.add_argument("--merge", action="store_true")
    a = p.parse_args()
    all_accs = pd.read_csv(CATALOG, sep="\t", usecols=["study_accession"]).study_accession.tolist()
    accs = [x.strip() for x in a.studies.split(",")] if a.studies else all_accs
    res = resolve(accs)
    res.to_csv(OUT if not a.studies else R + "host_bodysite_resolved_test.tsv",
               sep="\t", index=False)
    print("\nhost_final_from:", res.host_final_from.value_counts().to_dict())
    print("body_site_final:", res.body_site_final.value_counts().to_dict())
    if a.merge and not a.studies:
        merge(res)
