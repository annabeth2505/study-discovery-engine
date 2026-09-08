"""
extract_host_sources.py
For each study, record WHERE its host can be found -- the ENA deposit, the paper title,
the paper abstract -- and WHICH host each source names.

The question this answers: for NO_HOST_DATA studies (no sample-level ENA host), is the
host RECOVERABLE from the paper, is it a human-gut study whose host is implicit, or is
it genuinely unknown?

GRAIN: ENA host is per-sample and already tracked in samples.tsv via host_resolved_from.
Title and abstract are per-STUDY facts (one paper per study), so those columns live here,
at study level. Only ~561 studies have a linked paper; the rest are deposit-only and get
title/abstract = N/A rather than a fabricated blank.

Input : ../results/yes_catalog.tsv            (paper_title, abstract, pmid, host_species,
                                               biome_labels)
        ../results/study_sample_patterns.tsv  (sample-level host completeness)
Output: ../results/host_sources.tsv           (one row per study)

Claude extracts the host NAMED in the title and in the abstract, reusing the
false-positive guards from detect_discrepancies.py (model organisms, general-vs-specific).
The three-source synthesis is done in CODE, not by the model -- the model only reads text.
Checkpointed + resumable.

Usage:
  python extract_host_sources.py --limit 20 --out ../results/host_sources_test.tsv
  python extract_host_sources.py
"""

import os
import json
import time
import argparse
import pandas as pd

CATALOG = "../results/yes_catalog.tsv"
PATTERNS = "../results/study_sample_patterns.tsv"
OUT = "../results/host_sources.tsv"
CHECKPOINT = "../results/host_source_verdicts.json"

MODEL = "claude-haiku-4-5-20251001"

SENTINELS = {"missing", "na", "n/a", "not applicable", "none", "null", "unknown",
             "not collected", "not provided", "nan", ""}

HUMAN_TERMS = ("human", "homo sapiens", "patient", "infant", "adult", "child",
               "men", "women", "volunteer", "participant")

COLS = ["study_accession", "has_paper",
        "host_ena_referenced", "host_ena_value",
        "host_title_referenced", "host_title_value",
        "host_abstract_referenced", "host_abstract_value",
        "is_model_organism", "model_organism_disease", "human_by_biome_label",
        "host_source_summary", "host_source_note",
        "sample_host_pattern", "n_host_populated", "biome_labels"]


def _clean(v):
    if v is None or pd.isna(v):
        return None
    s = str(v).strip()
    return None if not s or s.lower() in SENTINELS else s


def _is_biome_label(v):
    """'gut metagenome' is a biome, not a host organism."""
    s = (v or "").lower()
    return "metagenome" in s or "microbiome" in s or s in ("metagenomes", "feces")


def _looks_human(v):
    s = (v or "").lower()
    return any(t in s for t in HUMAN_TERMS)


# ----------------------------- Claude extraction -------------------------------

SYSTEM = (
    "You extract the HOST ORGANISM that a metagenomics paper studied, from its title and "
    "its abstract, separately. Report the host NAME, not a yes/no. "
    "CRITICAL RULES: "
    "(1) THE HOST IS WHERE THE SEQUENCED SAMPLES CAME FROM -- not every animal the paper "
    "mentions. This is the single most important rule. Many papers sequence a HUMAN cohort "
    "and then run a SEPARATE downstream experiment in mice to test a mechanism. Those mice "
    "are NOT the host: the metagenomes came from the humans. Signals that an animal is only "
    "a downstream validation and NOT the host: 'finally, we show ... in mouse models', "
    "'further validated by FMT in mice', 'we confirmed the effect by administering X to "
    "mice', mice appearing only in the last sentence after a human cohort was described. "
    "In that case report the HUMAN cohort as the host and set is_model_organism FALSE. "
    "(2) MODEL ORGANISMS: when the SEQUENCED SAMPLES THEMSELVES come from the model animal "
    "(e.g. 'we constructed three mouse models and performed metagenomic sequencing', "
    "'we employed germ-free and SPF murine models ... to investigate the microbiota'), the "
    "HOST IS THE MODEL ANIMAL, not the human whose disease is modelled. Set "
    "is_model_organism true, put the DISEASE BEING MODELLED in model_organism_disease "
    "(e.g. 'colitis', 'acute pancreatitis', 'IBD', 'hepatocellular carcinoma'), and write "
    "the host as 'mouse (model of <disease>)' -- always name the specific disease if the "
    "title or abstract gives it, never the generic 'model of human disease'. "
    "(3) If the paper sequences BOTH a human cohort AND animals, report both, human first, "
    "e.g. 'human; mouse (model of type 2 diabetes)'. "
    "(4) GENERAL vs SPECIFIC: record the host at EXACTLY the specificity the paper uses. "
    "If the paper says 'ruminant', record 'ruminant' -- do NOT invent a species, and do "
    "NOT treat a general term as multiple species. If it says 'Bubalus bubalis', record "
    "that. Never expand a general term into a list. "
    "(5) NO GUESSING: if that field names no host organism, return an empty string for it. "
    "An empty string is a correct answer; a guess is not. "
    "(6) The host is the ANIMAL the samples came FROM, never the microbes in the sample. "
    "Bacterial or viral names are never the host. "
    "(7) In-vitro / bioreactor / environmental studies with no animal have no host: empty. "
    "Return ONLY valid JSON."
)

PROMPT = """Study {acc}.

=== TITLE ===
{title}

=== ABSTRACT ===
{abstract}

Extract the host organism named in the TITLE, and the host organism named in the
ABSTRACT, independently. If a field names no host, return "" for it -- do not guess,
and do not copy an answer across from the other field.

Apply the model-organism and general-vs-specific rules from the system prompt.

Before answering, ask yourself: did the SEQUENCED SAMPLES come from this animal, or is
the animal only a downstream experiment used to test a mechanism? Only the organism the
samples came from is the host.

If it IS a model-organism study, name the DISEASE being modelled -- "mouse (model of
colitis)", not "mouse (model of human disease)".

Return ONLY:
{{"title_host":"","abstract_host":"","is_model_organism":false,
"model_organism_disease":"","note":"one short sentence saying WHERE the sequenced samples came from"}}"""


def parse_json(text):
    """Pull the JSON object out of a reply that may carry reasoning around it.

    Asking the model to weigh 'sequenced samples vs downstream experiment' makes it
    reason out loud before answering, so the response is not always bare JSON.
    """
    t = text.strip().replace("```json", "").replace("```", "").strip()
    try:
        return json.loads(t)
    except json.JSONDecodeError:
        pass
    start, end = t.find("{"), t.rfind("}")
    if start != -1 and end > start:
        return json.loads(t[start:end + 1])   # raises if this is not JSON either
    raise json.JSONDecodeError("no JSON object in reply", t, 0)


def call_llm(client, prompt, retries=6):
    import anthropic
    for a in range(retries):
        try:
            m = client.messages.create(model=MODEL, max_tokens=900, system=SYSTEM,
                                       messages=[{"role": "user", "content": prompt}])
            return parse_json(m.content[0].text)
        except anthropic.RateLimitError:
            time.sleep(13 * (a + 1))
        except json.JSONDecodeError:
            return {"title_host": None, "abstract_host": None, "is_model_organism": None,
                    "model_organism_disease": None, "note": "PARSE_ERROR"}
        except Exception as e:
            return {"title_host": None, "abstract_host": None, "is_model_organism": None,
                    "model_organism_disease": None, "note": f"ERROR: {str(e)[:110]}"}
    return {"title_host": None, "abstract_host": None, "is_model_organism": None,
            "model_organism_disease": None, "note": "ERROR: max retries"}


def trunc(x, n):
    s = _clean(x) or "(none)"
    return s[:n] + "..." if len(s) > n else s


# ------------------------------ synthesis (code) -------------------------------

def synthesize(row, verdict):
    """Combine the three sources into one picture. Deliberately NOT done by the LLM."""
    # ---- ENA ----
    hs = _clean(row.get("host_species"))
    ena_val = hs if (hs and not _is_biome_label(hs)) else None
    n_pop = row.get("n_host_populated")
    n_pop = 0 if pd.isna(n_pop) else int(n_pop)
    # sample-level hosts are ENA-derived too, and are the better value when present
    if ena_val is None and n_pop > 0:
        dl = _clean(row.get("distinct_host_list"))
        if dl:
            ena_val = dl
    ena_ref = bool(ena_val)

    # ---- paper ----
    has_paper = bool(_clean(row.get("paper_title")) or _clean(row.get("abstract")))
    if verdict is None:
        t_val = a_val = None
        model_org = False
        disease = ""
        note = "" if has_paper else "no linked paper -- deposit-only study"
    else:
        t_val = _clean(verdict.get("title_host"))
        a_val = _clean(verdict.get("abstract_host"))
        model_org = bool(verdict.get("is_model_organism"))
        disease = _clean(verdict.get("model_organism_disease")) or ""
        note = _clean(verdict.get("note")) or ""
        if verdict.get("note") in ("PARSE_ERROR",) or str(verdict.get("note", "")).startswith("ERROR"):
            note = f"NOT_CHECKED: {verdict.get('note')}"

    paper_named = bool(t_val or a_val)

    # ---- human-by-biome-label ----
    # A model-organism study is NOT human-by-biome: its host is the model animal.
    biome = _clean(row.get("biome_labels")) or ""
    human_signal = (not model_org) and (
        "human" in biome.lower() or _looks_human(a_val) or _looks_human(t_val))
    human_by_biome = bool((not ena_ref) and human_signal)

    # ---- summary ----
    if human_by_biome:
        summary = "human-by-biome-label"
    elif ena_ref and paper_named:
        summary = "ENA + paper"
    elif paper_named:
        summary = "paper only"
    elif ena_ref:
        summary = "ENA only"
    else:
        summary = "unknown"

    if model_org and note:
        tag = f"model of {disease}" if disease else "model organism, disease unnamed"
        note = f"[{tag}] {note}"

    return {
        "study_accession": row["study_accession"],
        "has_paper": has_paper,
        "host_ena_referenced": ena_ref,
        "host_ena_value": ena_val or "",
        "host_title_referenced": bool(t_val),
        "host_title_value": t_val or ("N/A" if not has_paper else ""),
        "host_abstract_referenced": bool(a_val),
        "host_abstract_value": a_val or ("N/A" if not has_paper else ""),
        "is_model_organism": model_org,
        "model_organism_disease": disease,
        "human_by_biome_label": human_by_biome,
        "host_source_summary": summary,
        "host_source_note": note,
        "sample_host_pattern": row.get("host_pattern_confirmed"),
        "n_host_populated": n_pop,
        "biome_labels": biome,
    }


# ----------------------------------- run ---------------------------------------

def run(limit=None, studies=None, out=OUT, checkpoint=CHECKPOINT):
    import anthropic
    from pathlib import Path
    from dotenv import load_dotenv
    load_dotenv(Path(__file__).resolve().parent.parent / ".env")
    client = anthropic.Anthropic(api_key=os.getenv("ANTHROPIC_API_KEY"))

    cat = pd.read_csv(CATALOG, sep="\t", low_memory=False)
    pat = pd.read_csv(PATTERNS, sep="\t")
    df = cat.merge(pat[["study_accession", "host_pattern_confirmed",
                        "n_host_populated", "distinct_host_list"]],
                   on="study_accession", how="left")

    if studies:
        df = df[df.study_accession.isin(studies)]
    if limit:
        df = df.head(limit)

    # only studies with a paper go to the LLM; deposit-only studies have nothing to read
    need = df[df.apply(lambda r: bool(_clean(r.get("paper_title")) or
                                      _clean(r.get("abstract"))), axis=1)]

    done = {}
    if os.path.exists(checkpoint):
        with open(checkpoint) as f:
            done = json.load(f)
        print(f"Resuming -- {len(done)} papers already extracted")

    todo = [r for _, r in need.iterrows() if r.study_accession not in done]
    print(f"Studies: {len(df)}   with a paper: {len(need)}   to extract: {len(todo)}")

    for i, r in enumerate(todo):
        if i % 10 == 0:
            print(f"  {i}/{len(todo)}")
            with open(checkpoint, "w") as f:
                json.dump(done, f)
        done[r.study_accession] = call_llm(client, PROMPT.format(
            acc=r.study_accession, title=trunc(r.get("paper_title"), 400),
            abstract=trunc(r.get("abstract"), 3500)))
        time.sleep(0.3)

    with open(checkpoint, "w") as f:
        json.dump(done, f)

    rows = [synthesize(r, done.get(r.study_accession)) for _, r in df.iterrows()]
    res = pd.DataFrame(rows)[COLS].sort_values("study_accession").reset_index(drop=True)
    res.to_csv(out, sep="\t", index=False)

    print(f"\n{'='*58}\nHOST SOURCE EXTRACTION COMPLETE\n{'='*58}")
    print(f"Studies: {len(res)}")
    print("\nhost_source_summary:")
    for k, v in res.host_source_summary.value_counts().items():
        print(f"  {k:22s} {v:5d}")
    print(f"\nmodel-organism studies      : {int(res.is_model_organism.sum())}")
    print(f"human-by-biome-label        : {int(res.human_by_biome_label.sum())}")
    nc = int(res.host_source_note.astype(str).str.startswith("NOT_CHECKED").sum())
    print(f"NOT_CHECKED (LLM failure)   : {nc}")

    rescue = res[(~res.host_ena_referenced) & (res.host_source_summary.isin(
        ["paper only", "human-by-biome-label"]))]
    print(f"\nRESCUED (no ENA host, host found elsewhere): {len(rescue)}")
    print(f"\nSaved -> {out}")
    return res


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--limit", type=int, default=None)
    p.add_argument("--studies", type=str, default=None)
    p.add_argument("--out", type=str, default=OUT)
    p.add_argument("--checkpoint", type=str, default=CHECKPOINT)
    a = p.parse_args()
    run(limit=a.limit,
        studies=[s.strip() for s in a.studies.split(",")] if a.studies else None,
        out=a.out, checkpoint=a.checkpoint)
