"""
classify_llm.py
Confidence-tiered semantic classifier for animal gut WGS metagenomics studies.

Replaces keyword matching with an LLM judgment that:
  - reads ALL available evidence (paper title, ENA study title, biome label,
    host, body site, sample title, isolation source, library fields, abstract)
  - returns a verdict WITH a confidence level AND the specific evidence that
    drove it, so every call is auditable
  - is forced to distinguish EXPLICIT evidence ("sample_title = cecal content")
    from INFERENCE ("host is a penguin, so probably gut") and mark the latter
    low-confidence instead of presenting a guess as a fact

Inputs (all on disk, so this is resumable):
  - ../results/pubmed_new_classified.json  (the fetch results: which accessions
    returned runs + their library fields, from the prior step)
  - ../results/batch*.csv                  (Alan's tool output: paper titles/pmids)

Output:
  - ../results/llm_verdicts.json           (checkpoint, resumable)
  - ../results/llm_classified.tsv          (final table with evidence + tier)

Rate-limit safe: exponential backoff on 429. Set ADD_CREDITS_FOR_SPEED note below.
"""

import os
import io
import json
import time
import glob
import re
import requests
import pandas as pd
import anthropic
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

MODEL = "claude-haiku-4-5-20251001"   # cheap; ~$0.35 for ~275 studies
ENA_URL = "https://www.ebi.ac.uk/ena/portal/api/search"
FETCH_JSON = "../results/pubmed_new_classified.json"
CHECKPOINT = "../results/llm_verdicts.json"
OUT_PATH = "../results/llm_classified.tsv"

client = anthropic.Anthropic()


def make_session():
    s = requests.Session()
    retry = Retry(total=5, backoff_factor=0.5,
                  status_forcelist=[429, 500, 502, 503, 504], allowed_methods=["GET"])
    s.mount("https://", HTTPAdapter(max_retries=retry))
    return s


SESSION = make_session()

SYSTEM = (
    "You are a metagenomics data curator. You decide whether a sequencing study is a "
    "SHOTGUN/WGS METAGENOMIC study of an ANIMAL's GUT microbiome (gut, intestinal, "
    "fecal/faecal, cecal, rumen, hindgut, digesta; cloacal/rectal count as gut-proxy). "
    "Base your decision ONLY on the evidence given. If body site or method is NOT "
    "explicitly stated and you infer it from the host organism or typical practice, you "
    "MUST set confidence='low' and body_site_signal='inferred' and say so. Never present "
    "an inference as a stated fact. If evidence is insufficient, return UNCERTAIN rather "
    "than guessing. Return ONLY valid JSON, no other text."
)

PROMPT = """Evidence for study {acc}:
- Paper title: {paper_title}
- ENA study title: {ena_title}
- Biome label(s) (sample scientific_name): {biomes}
- Host organism(s): {hosts}
- Host body site: {body_sites}
- Sample title(s): {sample_titles}
- Isolation source(s): {isolation}
- library_source: {lib_source}  (METAGENOMIC=community DNA; GENOMIC=single-organism genome; TRANSCRIPTOMIC=RNA)
- library_strategy: {lib_strategy}  (WGS or METAGENOMIC = shotgun; AMPLICON = 16S, NOT shotgun)
- Abstract: {abstract}

Classify:
- verdict: one of
    GUT_WGS         = shotgun metagenomics of an animal gut/intestinal/fecal microbiome
    NOT_GUT         = metagenomics but of a non-gut site (skin/oral/soil/water/blood/etc.)
    NOT_METAGENOMIC = host genome, transcriptome, or 16S/amplicon (not shotgun community)
    UNCERTAIN       = evidence insufficient to decide
- confidence: high | medium | low
- body_site_signal: explicit | inferred | none
- evidence: quote the SPECIFIC field value(s) that drove the decision
- reasoning: one sentence

Return ONLY:
{{"verdict":"...","confidence":"...","body_site_signal":"...","evidence":"...","reasoning":"..."}}"""


def gather_ena_evidence(acc):
    """Pull the descriptive sample fields that indicate body site / biome."""
    params = {
        "result": "read_run",
        "query": f'study_accession="{acc}"',
        "fields": ("scientific_name,host_scientific_name,host_body_site,"
                   "sample_title,isolation_source,library_source,library_strategy"),
        "format": "tsv", "limit": 0,
    }
    try:
        r = SESSION.get(ENA_URL, params=params, timeout=60)
        if r.status_code == 200 and r.text.strip():
            return pd.read_csv(io.StringIO(r.text), sep="\t")
    except Exception:
        pass
    return pd.DataFrame()


def _distinct(df, col):
    if col not in df.columns:
        return ""
    vals = sorted(set(str(x) for x in df[col].dropna() if str(x).strip()))
    return "; ".join(vals[:8]) if vals else "(none)"


def call_llm(prompt, max_retries=6):
    for attempt in range(max_retries):
        try:
            msg = client.messages.create(
                model=MODEL, max_tokens=300, system=SYSTEM,
                messages=[{"role": "user", "content": prompt}],
            )
            text = msg.content[0].text.strip().replace("```json", "").replace("```", "").strip()
            return json.loads(text)
        except anthropic.RateLimitError:
            wait = 13 * (attempt + 1)
            print(f"    rate limited, waiting {wait}s...")
            time.sleep(wait)
        except json.JSONDecodeError:
            return {"verdict": "PARSE_ERROR", "confidence": "", "body_site_signal": "",
                    "evidence": "", "reasoning": "LLM returned non-JSON"}
        except Exception as e:
            return {"verdict": "ERROR", "confidence": "", "body_site_signal": "",
                    "evidence": "", "reasoning": str(e)[:120]}
    return {"verdict": "ERROR", "confidence": "", "body_site_signal": "",
            "evidence": "", "reasoning": "max retries exceeded"}


def run():
    # --- load fetch results: only classify accessions that returned runs ---
    with open(FETCH_JSON) as f:
        fetched = json.load(f)
    ok_accs = [a for a, v in fetched.items() if v.get("status") == "ok"]
    print(f"Studies to classify (had runs): {len(ok_accs)}")

    # --- load paper titles from Alan's tool output ---
    mmc = pd.concat([pd.read_csv(f) for f in sorted(glob.glob("../results/batch*.csv"))],
                    ignore_index=True)
    title_map, abstract_map = {}, {}
    for _, r in mmc.iterrows():
        for acc in re.findall(r'PRJ[ENDB][A-Z]\d+', str(r.get("accessions", ""))):
            title_map[acc] = str(r.get("title", "") or "")

    # --- resume ---
    verdicts = {}
    if os.path.exists(CHECKPOINT):
        with open(CHECKPOINT) as f:
            verdicts = json.load(f)
        print(f"Resuming — {len(verdicts)} already classified")

    todo = [a for a in ok_accs if a not in verdicts]
    print(f"Classifying {len(todo)} studies with {MODEL}...\n")

    for i, acc in enumerate(todo):
        if i % 10 == 0:
            print(f"  {i}/{len(todo)}")
            with open(CHECKPOINT, "w") as f:
                json.dump(verdicts, f)

        ena = gather_ena_evidence(acc)
        fields = fetched[acc]
        prompt = PROMPT.format(
            acc=acc,
            paper_title=title_map.get(acc, "(none)"),
            ena_title="(fetched below)" if ena.empty else _distinct(ena, "sample_title"),
            biomes=_distinct(ena, "scientific_name"),
            hosts=_distinct(ena, "host_scientific_name") or fields.get("host_species", ""),
            body_sites=_distinct(ena, "host_body_site") or fields.get("body_site", "(none)"),
            sample_titles=_distinct(ena, "sample_title"),
            isolation=_distinct(ena, "isolation_source"),
            lib_source=fields.get("library_source", "(none)"),
            lib_strategy=fields.get("library_strategy", "(none)"),
            abstract=abstract_map.get(acc, "(none available)"),
        )
        v = call_llm(prompt)
        v["host_species"] = fields.get("host_species", "")
        v["n_samples"] = fields.get("n_samples", 0)
        v["library_source"] = fields.get("library_source", "")
        v["library_strategy"] = fields.get("library_strategy", "")
        v["paper_title"] = title_map.get(acc, "")
        verdicts[acc] = v
        time.sleep(0.2)

    with open(CHECKPOINT, "w") as f:
        json.dump(verdicts, f)

    # --- assemble output + confidence tiers ---
    rows = [{"study_accession": a, **v} for a, v in verdicts.items()]
    out = pd.DataFrame(rows)

    def tier(r):
        if r["verdict"] == "GUT_WGS" and r["confidence"] in ("high", "medium"):
            return "CONFIRMED_GUT"
        if r["verdict"] == "GUT_WGS":                       # low-confidence gut
            return "LIKELY_GUT_REVIEW"
        if r["verdict"] == "UNCERTAIN":
            return "REVIEW"
        if r["verdict"] in ("NOT_GUT", "NOT_METAGENOMIC"):
            return "EXCLUDED"
        return "ERROR"

    out["tier"] = out.apply(tier, axis=1)
    out.to_csv(OUT_PATH, sep="\t", index=False)

    print(f"\n{'='*52}")
    print("LLM CLASSIFICATION — CONFIDENCE TIERS")
    print(f"{'='*52}")
    print(out["tier"].value_counts().to_string())
    print(f"\nVerdict x confidence:")
    print(pd.crosstab(out["verdict"], out["confidence"]).to_string())
    print(f"\nbody_site_signal (how many were INFERRED vs stated):")
    print(out["body_site_signal"].value_counts().to_string())

    conf = out[out["tier"] == "CONFIRMED_GUT"]
    print(f"\nCONFIRMED GUT (high/medium confidence): {len(conf)}")
    print(f"LIKELY GUT — needs review (low conf/inferred): {(out['tier']=='LIKELY_GUT_REVIEW').sum()}")
    print(f"UNCERTAIN — needs review: {(out['tier']=='REVIEW').sum()}")
    print(f"EXCLUDED: {(out['tier']=='EXCLUDED').sum()}")

    print(f"\nSaved -> {OUT_PATH}")
    return out


if __name__ == "__main__":
    run()