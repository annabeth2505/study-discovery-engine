import requests
import pandas as pd
import time 
from Bio import Entrez

ENA_PORTAL_URL = "https://www.ebi.ac.uk/ena/portal/api/search"

GENERIC_SCIENTIFIC_NAMES = [
    "metagenome", "gut metagenome", "microbial mat metagenome",
    "soil metagenome", "food metagenome", "environmental metagenome",
    "marine metagenome", "freshwater metagenome", "human gut metagenome",
    "mouse gut metagenome", "bovine gut metagenome"
]

def search_ena_studies(query, max_results=100):
    """
    Search ENA for studies matching a query.
    Returns list of unique study accessions.
    """
    params = {
        "result": "read_run",
        "query": query,
        "fields": "study_accession",
        "limit": max_results * 10,  # request more runs to get enough unique studies
        "format": "json"
    }

    try:
        response = requests.get(ENA_PORTAL_URL, params=params, timeout=30)
        if response.status_code != 200:
            print(f"ENA search failed: {response.status_code}")
            return []

        data = response.json()
        if not data:
            return []

        accessions = list({record["study_accession"] for record in data})
        print(f"Found {len(data)} runs, {len(accessions)} unique studies")
        time.sleep(0.3)
        return accessions[:max_results]

    except Exception as e:
        print(f"Error searching ENA: {e}")
        return []
    
def fetch_runs_for_study(study_accession):
    """  
    Fetch all runs for a study from ENA. 

    Args:
        study_accession: ENA/NCBI study accession 
    Returns:
        DataFrame of runs or none if fetch failed
    """
   
    params = {
        "result": "read_run",
        "query": f'study_accession="{study_accession}"',
        "fields": "run_accession,study_accession,sample_accession,scientific_name,host,library_strategy,library_source,fastq_ftp,\
        host_scientific_name,host_tax_id,host_body_site,disease,country,lat,lon,collection_date\
        ,library_strategy",
        "limit": 1000,
        "format": "json"
    }

    try:
        response = requests.get(ENA_PORTAL_URL, params = params, timeout = 30)
        if response.status_code != 200:
            print(f'ENA fetch failed for {study_accession}: {response.status_code}')
            return None 

        data = response.json()
        if not data:
            print(f'No runs found for {study_accession}')
            return None 

        df = pd.DataFrame(data)
        time.sleep(0.3)
        return df
    except Exception as e:
        print(f'Error fetching {study_accession}: {e}')
        return None




def resolve_host_species(runs_df):
    """
    Resolve the best available host species from a runs DataFrame.
    Checks host_scientific_name first, then host, then scientific_name.
    """
    generic_terms = [
        "metagenome", "gut metagenome", "microbial mat metagenome",
        "soil metagenome", "food metagenome", "environmental metagenome",
        "marine metagenome", "freshwater metagenome", "human gut metagenome",
        "mouse gut metagenome", "bovine gut metagenome"
    ]

    # try host_scientific_name first
    if "host_scientific_name" in runs_df.columns:
        hosts = runs_df["host_scientific_name"].dropna()
        hosts = hosts[hosts != ""]
        if len(hosts) > 0:
            return hosts.mode()[0]

    # try host field next if not found in host_scientific_name 
    if "host" in runs_df.columns:
        hosts = runs_df["host"].dropna()
        hosts = hosts[hosts != ""]
        if len(hosts) > 0:
            return hosts.mode()[0]

    # fall back to scientific_name if not generic
    sci_names = runs_df["scientific_name"].dropna()
    sci_name = sci_names.mode()[0] if len(sci_names) > 0 else None
    if sci_name and sci_name.lower() not in [g.lower() for g in generic_terms]:
        return sci_name

    return None

def fetch_pubmed_id(study_accession):
    """
    Search PubMed for a paper associated with a study accession.
    """
    try:
        handle = Entrez.esearch(
            db="pubmed",
            term=f"{study_accession}[All Fields]",
            retmax=1
        )
        record = Entrez.read(handle)
        handle.close()
        
        if record["IdList"]:
            return record["IdList"][0]
        return None
        
    except Exception as e:
        print(f"Error fetching PubMed ID for {study_accession}: {e}")
        return None
    
    
def fetch_study_origin(study_accession):
    """
    Fetches the origin of the study.
    """
    origin = {
        "accession": study_accession,
        "source": "ENA" if study_accession.startswith("PRJEB") else "NCBI",
        "title": None,
        "description": None,
    }

    params = {
        "result": "study",
        "query": f'study_accession="{study_accession}"',
        "fields": "study_accession,study_title,study_description",
        "limit": 1,
        "format": "json"
    }

    try:
        response = requests.get(ENA_PORTAL_URL, params=params, timeout=30)
        if response.status_code == 200 and response.json():
            data = response.json()[0]
            origin["title"] = data.get("study_title")
            origin["description"] = data.get("study_description")
        time.sleep(0.3)

    except Exception as e:
        print(f"Error fetching ENA metadata for {study_accession}: {e}")


    # fallback to NCBI if title still None
    if origin["title"] is None:
        try:
            handle = Entrez.esearch(
                db="bioproject",
                term=f"{study_accession}[Project Accession]",
                retmax=1
            )
            record = Entrez.read(handle)
            handle.close()

            if record["IdList"]:
                fetch_handle = Entrez.efetch(
                    db="bioproject",
                    id=record["IdList"][0],
                    rettype="xml",
                    retmode="xml"
                )
                import xml.etree.ElementTree as ET
                xml_data = fetch_handle.read()
                fetch_handle.close()
                root = ET.fromstring(xml_data)

                title_elem = root.find(".//ProjectDescr/Title") 

                # look specifically in ProjectDescr for the title
                title = root.find(".//ProjectDescr/Title")
                description = root.find(".//ProjectDescr/Description")

                if title is not None:
                    origin["title"] = title.text
                if description is not None:
                    origin["description"] = description.text
                
            time.sleep(0.3)

        except Exception as e:
            print(f"Error fetching NCBI title for {study_accession}: {e}")

    return origin

    
def fetch_pubmed_abstract(study_accession):
    """ 
    Fetches the abstract of a PubMed paper given its study accession. 
    """
    
    try:
        pubmed_id = fetch_pubmed_id(study_accession)
        
        if pubmed_id is None:
            return None
        
        handle = Entrez.efetch(
            db = 'pubmed',
            id = pubmed_id,
            rettype = 'abstract',
            retmode = 'text'
        )
        abstract = handle.read()
        handle.close()
        return abstract
    except Exception as e:
        print(f"No abstract found for PubMed ID {pubmed_id}: {e}")
    return None 


def fetch_pubmed_abstract_by_title(title):
    """
    Search PubMed by study title and return abstract.
    """

    if not title or not title.strip():
        return None
    
    try:
        handle = Entrez.esearch(
            db ="pubmed", 
            term = f'"{title}"[Title]', 
            retmax = 1
        )
        record = Entrez.read(handle)
        handle.close()
        time.sleep(0.3)

        if not record["IdList"]:
            return None
        
        pubmed_id = record["IdList"][0]

        # grab abstract using pubmed id 
        handle = Entrez.efetch(
            db = 'pubmed', 
            id = pubmed_id, 
            rettype = 'abstract', 
            retmode = 'text'
        )
        abstract = handle.read()
        handle.close()
        time.sleep(0.3)

        return abstract 
    except Exception as e:
        print(f"Error: {e}")
        return None 
    

def fetch_abstract_from_pmid(pubmed_id):
    """ 
    Fetches the abstract of a PubMed paper given its PMID. 
    """
    
    try:
        if pubmed_id is None:
            return None
        
        handle = Entrez.efetch(
            db = 'pubmed',
            id = pubmed_id,
            rettype = 'abstract',
            retmode = 'text'
        )
        abstract = handle.read()
        handle.close()
        return abstract
    except Exception as e:
        print(f"No abstract found for PubMed ID {pubmed_id}: {e}")
    return None 

def get_taxonomy(tax_id):
    """
    Get full taxonomy for a host_tax_id.
    Returns dict with kingdom, phylum, class, order, family, genus, species.
    """
    from Bio import Entrez as _Entrez
    _Entrez.email = Entrez.email
   
    if tax_id is None or str(tax_id) == 'nan':           
        return {"kingdom": None, "phylum": None, "class": None, 
            "order": None, "family": None, "genus": None, "species": None}
        
    try:
        converted_id = str(int(float(tax_id)))

        handle = Entrez.efetch(
            db = 'taxonomy', 
            id = converted_id,
            retmode = 'xml'
        )
        records = Entrez.read(handle)
        handle.close()
        time.sleep(0.3)

        taxonomy = {"kingdom": None, "phylum": None, "class": None,
                   "order": None, "family": None, "genus": None, "species": None}

        if records:
            taxonomy["species"] = records[0].get("ScientificName")
            lineage = records[0].get('LineageEx', [])
            for node in lineage:
                rank = node.get("Rank")
                name = node.get("ScientificName")
                if rank in taxonomy:
                    taxonomy[rank] = name
                
        return taxonomy
    
    except Exception as e:
        print(f"Error fetching taxonomy for tax_id {tax_id}: {e} ")
        return {"kingdom": None, "phylum": None, "class": None,
                "order": None, "family": None, "genus": None, "species": None}