"""
pubmed_search.py
Run the animal-gut-metagenomics PubMed query and return NEW papers only.

The query is the postdoc's, with ONE expansion: `metagenom*` added to the
method block (block 1). That rescues papers describing their method as
"metagenomic sequencing" / "shotgun metagenomics" / "metagenome" without ever
saying "shotgun" or "WGS" -- the change the postdoc predicted would 2-3x results.

Flow this file covers: esearch -> PMIDs -> dedupe against existing.
Downstream (separate steps): new PMIDs -> Alan's MMC tool -> accessions ->
fetch runs -> library_source/library_strategy WGS classification -> catalog diff.

Note: adding metagenom* also pulls in more 16S/amplicon papers (people say
"metagenomic" loosely). That's fine under keep-more-filter-later; the WGS vs 16S
call happens downstream at classification, not here.
"""

import time
import json
import glob
import pandas as pd
from Bio import Entrez

# ---------------------------------------------------------------------------
# The query. Block 1 (method) now includes metagenom*.  Everything else is the
# postdoc's query verbatim.
# ---------------------------------------------------------------------------
QUERY = """
(
  shotgun[Title/Abstract]
  OR "whole-genome sequencing"[Title/Abstract]
  OR "whole genome"[Title/Abstract]
  OR WGS[Title/Abstract]
  OR metagenom*[Title/Abstract]
)
AND
(
  "gut microbiom*"[Title/Abstract] OR "gut microbiot*"[Title/Abstract]
  OR "intestinal microbiom*"[Title/Abstract] OR "intestinal microbiot*"[Title/Abstract]
  OR "fecal microbiom*"[Title/Abstract] OR "fecal microbiot*"[Title/Abstract]
  OR "faecal microbiom*"[Title/Abstract] OR "faecal microbiot*"[Title/Abstract]
)
AND
(
  animals[MeSH] OR animal[Title/Abstract] OR animals[Title/Abstract]
  OR vertebrate*[Title/Abstract] OR invertebrate*[Title/Abstract] OR wildlife[Title/Abstract]
  OR mammal*[Title/Abstract] OR Mammalia[Title/Abstract]
  OR bird*[Title/Abstract] OR avian[Title/Abstract] OR Aves[Title/Abstract]
  OR insect*[Title/Abstract] OR Insecta[Title/Abstract]
  OR bee[Title/Abstract] OR bees[Title/Abstract] OR bumblebee*[Title/Abstract]
  OR wasp*[Title/Abstract] OR hornet*[Title/Abstract] OR ant[Title/Abstract] OR ants[Title/Abstract]
  OR beetle*[Title/Abstract] OR Coleoptera[Title/Abstract]
  OR butterfl*[Title/Abstract] OR moth[Title/Abstract] OR moths[Title/Abstract] OR Lepidoptera[Title/Abstract] OR caterpillar*[Title/Abstract]
  OR fly[Title/Abstract] OR flies[Title/Abstract] OR Diptera[Title/Abstract] OR mosquito*[Title/Abstract] OR midge*[Title/Abstract]
  OR cockroach*[Title/Abstract] OR termite*[Title/Abstract] OR Blattodea[Title/Abstract]
  OR cricket*[Title/Abstract] OR grasshopper*[Title/Abstract] OR locust*[Title/Abstract] OR Orthoptera[Title/Abstract]
  OR aphid*[Title/Abstract] OR hemipter*[Title/Abstract] OR leafhopper*[Title/Abstract] OR planthopper*[Title/Abstract]
  OR drosophila[Title/Abstract]
  OR dragonfl*[Title/Abstract] OR damselfl*[Title/Abstract] OR cicada*[Title/Abstract] OR lacewing*[Title/Abstract] OR earwig*[Title/Abstract]
  OR fish[Title/Abstract] OR fishes[Title/Abstract] OR teleost*[Title/Abstract] OR Actinopterygii[Title/Abstract]
  OR zebrafish[Title/Abstract]
  OR tuna[Title/Abstract] OR mackerel[Title/Abstract]
  OR cod[Title/Abstract] OR haddock[Title/Abstract] OR perch[Title/Abstract] OR Perciform*[Title/Abstract]
  OR eel[Title/Abstract] OR eels[Title/Abstract] OR Anguilla[Title/Abstract]
  OR flounder[Title/Abstract] OR halibut[Title/Abstract]
  OR herring[Title/Abstract] OR sardine*[Title/Abstract] OR anchov*[Title/Abstract]
  OR cichlid*[Title/Abstract] OR seahorse*[Title/Abstract]
  OR reptile*[Title/Abstract] OR Reptilia[Title/Abstract]
  OR lizard*[Title/Abstract] OR gecko*[Title/Abstract] OR iguana*[Title/Abstract] OR skink*[Title/Abstract] OR chameleon*[Title/Abstract] OR anole*[Title/Abstract] OR Squamata[Title/Abstract]
  OR snake*[Title/Abstract] OR serpent*[Title/Abstract] OR python*[Title/Abstract] OR boa[Title/Abstract] OR viper*[Title/Abstract] OR cobra*[Title/Abstract] OR rattlesnake*[Title/Abstract] OR colubrid*[Title/Abstract]
  OR turtle*[Title/Abstract] OR tortoise*[Title/Abstract] OR terrapin*[Title/Abstract] OR Testudines[Title/Abstract] OR chelonian*[Title/Abstract]
  OR crocodil*[Title/Abstract] OR alligator*[Title/Abstract] OR caiman*[Title/Abstract] OR gharial*[Title/Abstract] OR Crocodylia[Title/Abstract]
  OR tuatara*[Title/Abstract] OR Sphenodon[Title/Abstract]
  OR amphibian*[Title/Abstract] OR Amphibia[Title/Abstract]
  OR frog*[Title/Abstract] OR toad*[Title/Abstract] OR Anura[Title/Abstract] OR bullfrog*[Title/Abstract] OR treefrog*[Title/Abstract] OR xenopus[Title/Abstract]
  OR salamander*[Title/Abstract] OR newt*[Title/Abstract] OR axolotl*[Title/Abstract] OR Caudata[Title/Abstract] OR Urodela[Title/Abstract]
  OR caecilian*[Title/Abstract] OR Gymnophiona[Title/Abstract]
  OR shark*[Title/Abstract] OR Chondrichthy*[Title/Abstract] OR elasmobranch*[Title/Abstract]
  OR stingray*[Title/Abstract] OR manta[Title/Abstract] OR batoid*[Title/Abstract] OR chimaera*[Title/Abstract] OR ratfish[Title/Abstract]
  OR lamprey*[Title/Abstract] OR Petromyzon*[Title/Abstract] OR Hyperoartia[Title/Abstract]
  OR hagfish*[Title/Abstract] OR Myxini[Title/Abstract]
  OR lancelet*[Title/Abstract] OR amphioxus[Title/Abstract] OR Branchiostoma[Title/Abstract] OR Leptocardii[Title/Abstract] OR Cephalochordata[Title/Abstract]
  OR gastropod*[Title/Abstract] OR Gastropoda[Title/Abstract]
  OR snail*[Title/Abstract] OR slug*[Title/Abstract] OR conch*[Title/Abstract] OR limpet*[Title/Abstract] OR periwinkle*[Title/Abstract] OR whelk*[Title/Abstract] OR nudibranch*[Title/Abstract]
  OR Clitellata[Title/Abstract] OR oligochaet*[Title/Abstract] OR earthworm*[Title/Abstract] OR lumbric*[Title/Abstract] OR leech*[Title/Abstract] OR Hirudinea[Title/Abstract]
  OR bivalve*[Title/Abstract] OR Bivalvia[Title/Abstract] OR cockle*[Title/Abstract]
  OR Malacostraca[Title/Abstract] OR crustacean*[Title/Abstract]
  OR crab[Title/Abstract] OR crabs[Title/Abstract] OR lobster*[Title/Abstract] OR crayfish[Title/Abstract] OR krill[Title/Abstract] OR isopod*[Title/Abstract] OR amphipod*[Title/Abstract]
  OR holothurian*[Title/Abstract] OR Holothuroidea[Title/Abstract] OR "sea cucumber"[Title/Abstract] OR "sea cucumbers"[Title/Abstract]
  OR starfish[Title/Abstract] OR "sea star"[Title/Abstract] OR "sea stars"[Title/Abstract] OR Asteroidea[Title/Abstract]
  OR cephalopod*[Title/Abstract] OR Cephalopoda[Title/Abstract]
  OR octopus[Title/Abstract] OR octopi[Title/Abstract] OR octopuses[Title/Abstract] OR squid*[Title/Abstract] OR cuttlefish[Title/Abstract] OR nautilus*[Title/Abstract]
  OR arachnid*[Title/Abstract] OR Arachnida[Title/Abstract]
  OR spider*[Title/Abstract] OR scorpion*[Title/Abstract] OR tick[Title/Abstract] OR ticks[Title/Abstract] OR mite[Title/Abstract] OR mites[Title/Abstract] OR harvestman[Title/Abstract] OR harvestmen[Title/Abstract]
  OR Echinoidea[Title/Abstract] OR echinoid*[Title/Abstract] OR "sea urchin"[Title/Abstract] OR "sea urchins"[Title/Abstract] OR "sand dollar*"[Title/Abstract]
  OR Anthozoa[Title/Abstract] OR coral[Title/Abstract] OR corals[Title/Abstract] OR "sea anemone"[Title/Abstract] OR "sea anemones"[Title/Abstract]
  OR Sagittoidea[Title/Abstract] OR chaetognath*[Title/Abstract] OR "arrow worm"[Title/Abstract] OR "arrow worms"[Title/Abstract]
  OR Demospongiae[Title/Abstract] OR sponge[Title/Abstract] OR sponges[Title/Abstract] OR Porifera[Title/Abstract]
  OR Sipuncula[Title/Abstract] OR sipunculid*[Title/Abstract] OR "peanut worm"[Title/Abstract] OR "peanut worms"[Title/Abstract]
)
NOT humans[MeSH:noexp]
"""


def run_esearch(query=QUERY, retmax=100000):
    """Return the full list of PMIDs matching the query."""
    handle = Entrez.esearch(db="pubmed", term=query, retmax=retmax, retmode="xml")
    rec = Entrez.read(handle)
    handle.close()
    pmids = rec.get("IdList", [])
    total = int(rec.get("Count", len(pmids)))
    print(f"PubMed reports {total} total matches; retrieved {len(pmids)} PMIDs")
    if total > len(pmids):
        print(f"  (raise retmax above {retmax} to get all)")
    return pmids


def load_existing_pmids():
    """Collect PMIDs already processed, from the MMC result CSVs + postdoc CSV + acc_to_pmid."""
    existing = set()

    # MMC / pmid-ena result CSVs in results/
    for path in glob.glob("../results/pmid-ena-results*.csv"):
        try:
            df = pd.read_csv(path)
            for col in df.columns:
                if col.lower() in ("pmid", "pmids"):
                    existing |= {str(int(float(x))) for x in df[col].dropna()}
        except Exception as e:
            print(f"  skip {path}: {e}")

    # postdoc CSV
    for path in glob.glob("../results/*studies*.csv") + glob.glob("/mnt/user-data/uploads/*Studies*.csv"):
        try:
            df = pd.read_csv(path)
            if "PMID" in df.columns:
                existing |= {str(int(float(x))) for x in df["PMID"].dropna()}
        except Exception:
            pass

    # acc_to_pmid.json
    try:
        with open("../results/acc_to_pmid.json") as f:
            acc_to_pmid = json.load(f)
        existing |= {str(int(float(v))) for v in acc_to_pmid.values() if v}
    except Exception:
        pass

    print(f"Existing PMIDs already in pipeline: {len(existing)}")
    return existing


def run():
    pmids = run_esearch()
    pmid_set = {str(p) for p in pmids}

    existing = load_existing_pmids()
    new_pmids = sorted(pmid_set - existing)

    print(f"\n{'='*50}")
    print("PUBMED QUERY (metagenom* expanded)")
    print(f"{'='*50}")
    print(f"Query matches:        {len(pmid_set)}")
    print(f"Already processed:    {len(pmid_set & existing)}")
    print(f"NEW papers to fetch:  {len(new_pmids)}")

    with open("../results/pubmed_new_pmids.txt", "w") as f:
        f.write("\n".join(new_pmids))
    print(f"\nSaved {len(new_pmids)} new PMIDs -> pubmed_new_pmids.txt")
    return new_pmids


if __name__ == "__main__":
    run()