"""
classify_sweep_gut.py
Classify the 64 sweep studies as gut / not-gut / needs-abstract.

Signal per sample (run): the biome label (scientific_name), the host_body_site,
and the sample_title. A study's verdict aggregates its samples:
  - KEEP           : ANY sample looks gut/feces
  - DROP           : EVERY sample is clearly non-gut (skin/oral/soil/viral/...)
  - NEEDS_ABSTRACT : no clear gut sample, but at least one ambiguous
                     ("bird metagenome", bare species name) -> decide from abstract

Rationale for the asymmetric rule: a study that sampled both penguin gut and
penguin skin is still a legitimate gut study, so any gut sample => keep.

Outputs: results/sweep_gut_classified.tsv (one row per study, with evidence)
"""

import io
import re
import time
import requests
import pandas as pd
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

ENA_URL = "https://www.ebi.ac.uk/ena/portal/api/search"
SWEEP_PATH = "../results/sweep_new_studies_v2.tsv"
OUT_PATH = "../results/sweep_gut_classified.tsv"


# ---- resilient session (fixes the connection-reset drops) ----
def make_session():
    s = requests.Session()
    retry = Retry(total=5, backoff_factor=0.5,
                  status_forcelist=[429, 500, 502, 503, 504],
                  allowed_methods=["GET"])
    s.mount("https://", HTTPAdapter(max_retries=retry))
    return s


SESSION = make_session()


# ---- gut vs not-gut term patterns (word-boundary aware) ----
GUT_PATTERNS = [
    r"\bgut\b", r"\bfaec", r"\bfec", r"\bfeces\b", r"\bintestin", r"\bcec",
    r"\bcaec", r"\bcolon\b", r"\brumen\b", r"\brumin", r"\bdigesta\b",
    r"\bgastrointestinal\b", r"\bgastro-intestinal\b", r"\bstool\b",
    r"\bhindgut\b", r"\bmidgut\b", r"\bforegut\b", r"\bdung\b", r"\bmanure\b",
    r"\brectal\b", r"\brectum\b", r"\bcloaca", r"\banus\b", r"\bdigestive\b",
    r"\bcoprolite\b", r"\bcopro",
]
NOT_GUT_PATTERNS = [
    r"\bskin\b", r"\bderm", r"\boral\b", r"\bmouth\b", r"\btongue\b",
    r"\bbuccal\b", r"\bsaliva", r"\bnasal\b", r"\bnostril\b", r"\bnares\b",
    r"\bnasopharyng", r"\bfeather", r"\bplumage\b", r"\bwing\b", r"\bfur\b",
    r"\bpelage\b", r"\bsoil\b", r"\bsediment\b", r"\bseawater\b",
    r"\bfreshwater\b", r"\bmarine\b", r"\bblood\b", r"\bplasma\b", r"\bserum\b",
    r"\bvirome\b", r"\bviral\b", r"\bvirus\b", r"\bphage\b", r"\brespiratory\b",
    r"\blung\b", r"\bocular\b", r"\beye\b", r"\bgenital\b", r"\bvaginal\b",
]
GUT_RE = re.compile("|".join(GUT_PATTERNS), re.I)
NOTGUT_RE = re.compile("|".join(NOT_GUT_PATTERNS), re.I)


def classify_sample(*texts):
    """Classify one sample from its biome label / body site / title."""
    blob = " ".join(str(t) for t in texts if pd.notna(t)).lower()
    gut = bool(GUT_RE.search(blob))
    notgut = bool(NOTGUT_RE.search(blob))
    if gut:
        return "gut"          # gut wins even if other terms also present
    if notgut:
        return "not_gut"
    return "ambiguous"


def fetch_runs(study_acc):
    params = {
        "result": "read_run",
        "query": f'study_accession="{study_acc}"',
        "fields": "study_accession,host_scientific_name,host_body_site,scientific_name,sample_title",
        "format": "tsv", "limit": 0,
    }
    try:
        r = SESSION.get(ENA_URL, params=params, timeout=60)
        if r.status_code == 200 and r.text.strip():
            return pd.read_csv(io.StringIO(r.text), sep="\t")
    except Exception as e:
        print(f"  fetch failed {study_acc}: {e}")
    return pd.DataFrame()


def run():
    sweep = pd.read_csv(SWEEP_PATH, sep="\t")
    clade_map = dict(zip(sweep["study_accession"], sweep.get("_clade", "")))

    results = []
    accs = sweep["study_accession"].unique()
    print(f"Classifying {len(accs)} studies...")

    for i, acc in enumerate(accs):
        if i % 10 == 0:
            print(f"  {i}/{len(accs)}")
        runs = fetch_runs(acc)
        time.sleep(0.25)
        if len(runs) == 0:
            results.append({"study_accession": acc, "gut_verdict": "NO_DATA",
                            "clade": clade_map.get(acc, ""), "n_runs": 0,
                            "biome_labels": "", "body_sites": "", "host": ""})
            continue

        sample_verdicts = [
            classify_sample(r.get("scientific_name"), r.get("host_body_site"),
                            r.get("sample_title"))
            for _, r in runs.iterrows()
        ]

        if "gut" in sample_verdicts:
            verdict = "KEEP"
        elif "ambiguous" in sample_verdicts:
            verdict = "NEEDS_ABSTRACT"
        else:
            verdict = "DROP"

        biomes = sorted(set(str(x) for x in runs["scientific_name"].dropna()))
        sites = sorted(set(str(x) for x in runs["host_body_site"].dropna()))
        hosts = sorted(set(str(x) for x in runs["host_scientific_name"].dropna()))

        results.append({
            "study_accession": acc,
            "gut_verdict": verdict,
            "clade": clade_map.get(acc, ""),
            "host": "; ".join(hosts),
            "n_runs": len(runs),
            "n_gut_samples": sample_verdicts.count("gut"),
            "n_notgut_samples": sample_verdicts.count("not_gut"),
            "n_ambiguous_samples": sample_verdicts.count("ambiguous"),
            "biome_labels": "; ".join(biomes),
            "body_sites": "; ".join(sites),
        })

    out = pd.DataFrame(results)
    out.to_csv(OUT_PATH, sep="\t", index=False)

    print(f"\n{'='*50}")
    print("GUT CLASSIFICATION OF 64 SWEEP STUDIES")
    print(f"{'='*50}")
    print(out["gut_verdict"].value_counts().to_string())

    keep = out[out["gut_verdict"] == "KEEP"]
    print(f"\nCONFIRMED GUT: {len(keep)} studies")
    print(f"\nConfirmed-gut by clade:")
    print(keep["clade"].value_counts().to_string())

    print(f"\nNEEDS_ABSTRACT (decide from literature):")
    na = out[out["gut_verdict"] == "NEEDS_ABSTRACT"]
    print(na[["study_accession", "clade", "host", "biome_labels"]].to_string(index=False))

    print(f"\nDROPPED (clearly non-gut):")
    drop = out[out["gut_verdict"] == "DROP"]
    print(drop[["study_accession", "clade", "host", "biome_labels"]].to_string(index=False))

    print(f"\nSaved -> {OUT_PATH}")
    return out


if __name__ == "__main__":
    run()