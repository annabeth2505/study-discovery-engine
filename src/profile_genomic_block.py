"""
profile_genomic_block.py
Characterize the library_source=GENOMIC samples in the sample-level catalog.

Question for Sam: GENOMIC samples sit inside metagenomics deposits. Are they cultured
ISOLATE genomes (a different data type, arguably out of scope), or something else?
This is a PROFILING step only -- nothing is dropped, reclassified, or written back
into samples.tsv / yes_catalog.tsv.

Each GENOMIC sample is bucketed by what its scientific_name actually names:
  HOST_GENOME       -- names the study's own host animal (Salmo salar in a salmon
                       study) => host DNA sequencing, not gut metagenome, not isolate
  METAGENOME_LABEL  -- "gut metagenome", "mouse gut metagenome" => a metagenome sample
                       whose library_source is simply labelled GENOMIC
  MICROBIAL_ISOLATE -- a specific named microbe (Escherichia coli, Bacteroides spp.)
                       => the cultured-isolate hypothesis, where it holds
  UNIDENTIFIED      -- "unidentified" / blank

Output: ../results/genomic_block_profile.tsv (one row per study containing GENOMIC samples)
"""

import pandas as pd

SAMPLES = "../results/samples.tsv"
OUT = "../results/genomic_block_profile.tsv"

SENTINELS = {"missing", "na", "n/a", "none", "null", "unknown", "not collected",
             "not provided", "unidentified", "uncultured organism", ""}


def _clean(v):
    if v is None or pd.isna(v):
        return None
    s = str(v).strip()
    return None if not s or s.lower() in SENTINELS else s


def _is_metagenome_label(v):
    s = (v or "").lower()
    return "metagenome" in s or "microbiome" in s


def bucket_row(sci, host_sci, study_hosts):
    """What does this GENOMIC sample's scientific_name actually name?"""
    sci = _clean(sci)
    if sci is None:
        return "UNIDENTIFIED"
    if _is_metagenome_label(sci):
        return "METAGENOME_LABEL"
    # the animal the study is about -> host DNA, not a microbial isolate
    if sci == _clean(host_sci) or sci in study_hosts:
        return "HOST_GENOME"
    return "MICROBIAL_ISOLATE"


def joined(series, n=6):
    vals = [v for v in series.dropna().astype(str) if v.strip()]
    uniq = sorted(set(vals))
    out = "; ".join(uniq[:n])
    return out + (f" (+{len(uniq) - n} more)" if len(uniq) > n else "")


def mix(series):
    vc = series.dropna().value_counts()
    return "; ".join(f"{k}:{v}" for k, v in vc.items())


def run():
    d = pd.read_csv(SAMPLES, sep="\t", dtype=str)
    g = d[d["library_source"] == "GENOMIC"].copy()

    # a study's own hosts, used to tell host DNA apart from a cultured microbe
    hosts_by_study = (d[d["host_resolved"].notna()]
                      .groupby("study_accession")["host_resolved"]
                      .apply(lambda s: set(x for x in s if _clean(x))).to_dict())

    g["genomic_kind"] = [
        bucket_row(sci, hs, hosts_by_study.get(acc, set()))
        for sci, hs, acc in zip(g["scientific_name"], g["host_scientific_name"],
                                g["study_accession"])
    ]

    rows = []
    for acc, sub in g.groupby("study_accession"):
        study = d[d["study_accession"] == acc]
        src = study["library_source"].value_counts()
        n_meta = int(src.get("METAGENOMIC", 0))
        kinds = sub["genomic_kind"].value_counts()
        dominant = kinds.index[0]

        rows.append({
            "study_accession": acc,
            "n_samples_total": len(study),
            "n_genomic": len(sub),
            "pct_genomic": round(100 * len(sub) / len(study), 1),
            "n_metagenomic": n_meta,
            "deposit_type": "GENOMIC_ONLY" if n_meta == 0 else "MIXED",
            "genomic_dominant_kind": dominant,
            "n_host_genome": int(kinds.get("HOST_GENOME", 0)),
            "n_metagenome_label": int(kinds.get("METAGENOME_LABEL", 0)),
            "n_microbial_isolate": int(kinds.get("MICROBIAL_ISOLATE", 0)),
            "n_unidentified": int(kinds.get("UNIDENTIFIED", 0)),
            "genomic_strategy_mix": mix(sub["library_strategy"]),
            "genomic_scientific_names": joined(sub["scientific_name"]),
            "study_seq_type_mix": mix(study["library_source"]),
            "example_sample_titles": joined(sub["sample_title"], n=3),
        })

    prof = (pd.DataFrame(rows)
            .sort_values("n_genomic", ascending=False)
            .reset_index(drop=True))
    prof.to_csv(OUT, sep="\t", index=False)

    # ---------------- summary ----------------
    tot = len(g)
    W = 60
    print("=" * W)
    print("GENOMIC BLOCK PROFILE")
    print("=" * W)
    print(f"GENOMIC samples: {tot} across {g.study_accession.nunique()} studies "
          f"({100*tot/len(d):.1f}% of the {len(d)} sample catalog)")
    print()

    print("-- concentration --")
    top5 = prof.head(5)
    print(f"top 5 studies hold {top5.n_genomic.sum()} samples "
          f"({100*top5.n_genomic.sum()/tot:.0f}% of the block)")
    print(prof.head(10)[["study_accession", "n_genomic", "n_samples_total",
                         "pct_genomic", "deposit_type",
                         "genomic_dominant_kind"]].to_string(index=False))
    print()

    print("-- what the GENOMIC samples actually are --")
    kc = g["genomic_kind"].value_counts()
    for k, v in kc.items():
        print(f"  {k:18s} {v:6d}  ({100*v/tot:4.1f}%)")
    print()

    print("-- deposit type --")
    dt = prof["deposit_type"].value_counts()
    for k, v in dt.items():
        n = prof.loc[prof.deposit_type == k, "n_genomic"].sum()
        print(f"  {k:14s} {v:3d} studies, {n:5d} GENOMIC samples")
    print()

    print("-- library_strategy within the GENOMIC block --")
    print("  " + mix(g["library_strategy"]))
    print()

    iso = prof[prof.genomic_dominant_kind == "MICROBIAL_ISOLATE"]
    print("-- studies that look like genuine ISOLATE-genome deposits --")
    if iso.empty:
        print("  none dominant; isolates appear only as a minority tail")
    else:
        print(iso[["study_accession", "n_genomic", "n_samples_total", "deposit_type",
                   "genomic_scientific_names"]].to_string(index=False))
    print(f"\nTotal MICROBIAL_ISOLATE samples anywhere: "
          f"{int(prof.n_microbial_isolate.sum())}")
    print(f"\nSaved -> {OUT}")
    return prof


if __name__ == "__main__":
    run()
