import os
import re
import time
from pathlib import Path

import pandas as pd
import anthropic
from dotenv import load_dotenv

from src.ena_fetcher import fetch_pubmed_abstract, fetch_study_origin
from src.fetcher import configure_entrez


# (case-insensitive regex, canonical scientific binomial)
# common-name patterns map back to a single canonical species so title and
# abstract mentions can be intersected even when phrased differently.
SPECIES_PATTERNS = [
    (r"\bHomo sapiens\b|\bhumans?\b|\bpatients?\b|\binfants?\b|\bchild(?:ren)?\b|\badults?\b|\bvolunteers?\b", "Homo sapiens"),
    (r"\bSus scrofa(?: domesticus)?\b|\bpigs?\b|\bswine\b|\bporcine\b|\bpiglets?\b|\bsows?\b|\bboars?\b|\bhogs?\b", "Sus scrofa"),
    (r"\bGallus gallus(?: domesticus)?\b|\bchickens?\b|\bbroilers?\b|\bhens?\b|\bpoultry\b|\blaying hens?\b", "Gallus gallus"),
    (r"\bMus musculus\b|\bmouse\b|\bmice\b|\bmurine\b", "Mus musculus"),
    (r"\bBos taurus\b|\bcattle\b|\bcows?\b|\bbovine\b|\bcalves\b|\bcalf\b|\bbulls?\b|\bheifers?\b|\bsteers?\b|\bdairy (?:cows?|cattle)\b|\bbeef cattle\b", "Bos taurus"),
    (r"\bRattus norvegicus\b|\brats?\b", "Rattus norvegicus"),
    (r"\bSalmo salar\b|\bAtlantic salmon\b|\bsalmon\b", "Salmo salar"),
    (r"\bDanio rerio\b|\bzebrafish\b", "Danio rerio"),
    (r"\bFelis catus\b|\bcats?\b|\bfeline\b|\bkittens?\b", "Felis catus"),
    (r"\bCanis lupus familiaris\b|\bdogs?\b|\bcanine\b|\bpuppies\b|\bpuppy\b", "Canis lupus familiaris"),
    (r"\bMeleagris gallopavo\b|\bturkeys?\b", "Meleagris gallopavo"),
    (r"\bEquus caballus\b|\bhorses?\b|\bequine\b|\bfoals?\b|\bmares?\b|\bstallions?\b", "Equus caballus"),
    (r"\bCapra hircus\b|\bgoats?\b|\bcaprine\b", "Capra hircus"),
    (r"\bOvis aries\b|\bsheep\b|\bovine\b|\blambs?\b|\bewes?\b|\brams?\b", "Ovis aries"),
    (r"\bOryctolagus cuniculus\b|\brabbits?\b", "Oryctolagus cuniculus"),
    (r"\bAiluropoda melanoleuca\b|\bgiant pandas?\b", "Ailuropoda melanoleuca"),
    (r"\bMicropterus salmoides\b|\blargemouth bass\b", "Micropterus salmoides"),
    (r"\bOchotona curzoniae\b|\bplateau pikas?\b|\bpikas?\b", "Ochotona curzoniae"),
    (r"\bAnas platyrhynchos\b|\bducks?\b|\bmallards?\b", "Anas platyrhynchos"),
    (r"\bApis mellifera\b|\bhoney ?bees?\b", "Apis mellifera"),
    (r"\bDrosophila melanogaster\b|\bfruit fl(?:y|ies)\b", "Drosophila melanogaster"),
    (r"\bOncorhynchus mykiss\b|\brainbow trout\b", "Oncorhynchus mykiss"),
    (r"\bCyprinus carpio\b|\bcommon carp\b", "Cyprinus carpio"),
    (r"\bOreochromis niloticus\b|\bNile tilapia\b|\btilapia\b", "Oreochromis niloticus"),
    (r"\bLitopenaeus vannamei\b|\bPenaeus vannamei\b|\bwhiteleg shrimp\b", "Litopenaeus vannamei"),
    (r"\bPan troglodytes\b|\bchimpanzees?\b", "Pan troglodytes"),
    (r"\bGorilla gorilla\b|\bgorillas?\b", "Gorilla gorilla"),
    (r"\bMacaca mulatta\b|\brhesus (?:monkeys?|macaques?)\b", "Macaca mulatta"),
]


def find_species_in_text(text):
    """Return set of canonical scientific names matched in text."""
    if not isinstance(text, str) or not text:
        return set()
    found = set()
    for pattern, canonical in SPECIES_PATTERNS:
        if re.search(pattern, text, flags=re.IGNORECASE):
            found.add(canonical)
    return found


def derive_via_regex(title, abstract):
    """
    Returns (species, source) or (None, None) if regex can't confirm.

    source values:
      - 'regex_confirmed'              both title and abstract contain the same single species
      - 'regex_title_only_no_abstract' abstract unavailable, title has exactly one species
    """
    title_set = find_species_in_text(title)
    abs_set = find_species_in_text(abstract)

    has_abstract = isinstance(abstract, str) and abstract.strip()

    if not has_abstract:
        if len(title_set) == 1:
            return (next(iter(title_set)), "regex_title_only_no_abstract")
        return (None, None)

    intersection = title_set & abs_set
    if len(intersection) == 1:
        return (next(iter(intersection)), "regex_confirmed")
    return (None, None)


_LLM_SYSTEM = """You are a careful biologist. Given a study title and abstract, identify the SINGLE primary host species whose microbiome or biology is being characterized.

Rules:
- Return ONLY the scientific binomial (e.g. "Sus scrofa", "Gallus gallus", "Homo sapiens").
- If multiple hosts are studied equally, pick the most prominent one.
- Ignore species mentioned only as comparison, model system, reference, or pathogen — focus on the host being sampled/sequenced.
- If no specific animal host can be confidently identified, return exactly: UNKNOWN
- Output one line, no explanation, no punctuation, no quotes."""


def derive_via_llm(title, abstract, client, model="claude-haiku-4-5"):
    title = title or ""
    abstract_text = abstract if isinstance(abstract, str) and abstract.strip() else "(no abstract available)"
    user_msg = f"Title: {title}\n\nAbstract: {abstract_text}"

    resp = client.messages.create(
        model=model,
        max_tokens=32,
        system=[{"type": "text", "text": _LLM_SYSTEM, "cache_control": {"type": "ephemeral"}}],
        messages=[{"role": "user", "content": user_msg}],
    )
    text = resp.content[0].text.strip().splitlines()[0].strip().strip('"').strip("'")
    if not text or text.upper() == "UNKNOWN":
        return (None, "llm_unknown")
    return (text, "llm")


OUTPUT_COLUMNS = [
    "study_accession", "host_species", "derived_host_species", "derivation_source",
    "host_tax_id", "body_site", "country", "n_samples", "library_strategy", "is_animal",
    "title", "abstract",
]


def _setup_clients():
    load_dotenv(Path(__file__).resolve().parent.parent / ".env")
    configure_entrez()
    return anthropic.Anthropic(api_key=os.getenv("ANTHROPIC_API_KEY"))


def enrich_one(accession, client, existing_title=None):
    """Enrich a single accession. Returns dict with title, abstract, derived_host_species, derivation_source."""
    title = existing_title
    if not isinstance(title, str) or not title.strip():
        origin = fetch_study_origin(accession)
        title = origin.get("title")

    abstract = fetch_pubmed_abstract(accession)

    species, source = derive_via_regex(title, abstract)
    if species is None:
        try:
            species, source = derive_via_llm(title, abstract, client)
        except Exception as e:
            print(f"  [{accession}] LLM error: {e}")
            species, source = (None, "llm_error")

    return {
        "title": title,
        "abstract": abstract,
        "derived_host_species": species,
        "derivation_source": source,
    }


def enrich_accessions(accessions, output_csv=None, resume=True, checkpoint_every=25):
    """Enrich a list of bare accession codes. Returns DataFrame; optionally writes CSV."""
    client = _setup_clients()

    output_path = Path(output_csv) if output_csv else None
    if output_path:
        output_path.parent.mkdir(parents=True, exist_ok=True)

    rows = []
    done = set()
    if resume and output_path and output_path.exists():
        existing = pd.read_csv(output_path)
        rows = existing.to_dict("records")
        done = set(existing["study_accession"])
        print(f"Resuming: {len(done)}/{len(accessions)} already processed")

    for i, acc in enumerate(accessions):
        if acc in done:
            continue
        result = enrich_one(acc, client)
        rows.append({"study_accession": acc, **result})
        print(f"[{i+1}/{len(accessions)}] {acc} -> {result['derived_host_species']} ({result['derivation_source']})")
        if output_path and (i + 1) % checkpoint_every == 0:
            _write(rows, output_path)
        time.sleep(0.1)

    out_df = pd.DataFrame(rows)
    if output_path:
        _write(rows, output_path)
        print(f"\nWrote {len(rows)} rows -> {output_path}")
    return out_df


def enrich_catalog(input_tsv, output_csv, resume=True, checkpoint_every=25):
    """Enrich a TSV catalog (preserves all input columns; reuses existing 'title' column when present)."""
    df = pd.read_csv(input_tsv, sep="\t")
    output_csv = Path(output_csv)
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    client = _setup_clients()

    if resume and output_csv.exists():
        existing = pd.read_csv(output_csv)
        done = set(existing["study_accession"])
        rows = existing.to_dict("records")
        print(f"Resuming: {len(done)}/{len(df)} already processed")
    else:
        done = set()
        rows = []

    for i, row in df.iterrows():
        acc = row["study_accession"]
        if acc in done:
            continue

        result = enrich_one(acc, client, existing_title=row.get("title"))
        out = row.to_dict()
        out.update(result)
        rows.append(out)

        print(f"[{i+1}/{len(df)}] {acc} -> {result['derived_host_species']} ({result['derivation_source']})")
        if (i + 1) % checkpoint_every == 0:
            _write(rows, output_csv)
        time.sleep(0.1)

    _write(rows, output_csv)
    print(f"\nWrote {len(rows)} rows -> {output_csv}")
    return rows


def _write(rows, path):
    out_df = pd.DataFrame(rows)
    cols = [c for c in OUTPUT_COLUMNS if c in out_df.columns] + \
           [c for c in out_df.columns if c not in OUTPUT_COLUMNS]
    out_df[cols].to_csv(path, index=False)


def _parse_accessions_arg(raw_tokens):
    """Accepts list of strings; splits each on whitespace/commas/newlines."""
    out = []
    for tok in raw_tokens:
        for part in re.split(r"[\s,]+", tok.strip()):
            if part:
                out.append(part)
    return out


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Enrich study accessions with derived host species (regex + Claude Haiku).",
    )
    parser.add_argument("accessions", nargs="*",
                        help="Study accession codes (PRJNA*, PRJEB*, etc.). Comma- or space-separated.")
    parser.add_argument("--input-tsv", help="Path to a TSV catalog (alternative to listing accessions).")
    parser.add_argument("--output", default="results/host_enhanced.csv", help="Output CSV path.")
    parser.add_argument("--no-resume", action="store_true", help="Ignore any existing output CSV.")
    args = parser.parse_args()

    resume = not args.no_resume

    if args.input_tsv:
        enrich_catalog(args.input_tsv, args.output, resume=resume)
    elif args.accessions:
        accs = _parse_accessions_arg(args.accessions)
        if not accs:
            parser.error("No accession codes parsed from arguments.")
        enrich_accessions(accs, args.output, resume=resume)
    else:
        # default: full animal WGS catalog
        enrich_catalog("results/animal_wgs_catalog.tsv", "results/host_enhanced_animals.csv", resume=resume)
