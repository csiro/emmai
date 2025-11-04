#!/usr/bin/env python
import os
import re
import cobra
import logging
import requests
import pandas as pd
import pubchempy as pcp
from chemspipy import ChemSpider
from Bio import SeqIO
from Bio.SeqUtils import molecular_weight
from concurrent.futures import ThreadPoolExecutor, as_completed
import threading
import yaml
from typing import List

DEBUG = False

# pubchem_columns = ["cid", "cmpdname", "cmpdsynonym", "smiles"]
# pubchemdb = pd.read_csv("PubChem_compound_all_pathways.csv", usecols=pubchem_columns)
smiles_ref_db_columns: list[str] = [
    "kegg_metabolite_ID",
    "Metabolite_aliases",
    "BiGG_metabolite_name",
    "SMILES",
]
smiles_db = pd.read_csv("SMILES_reference_DB.csv", usecols=smiles_ref_db_columns)

inputs_path = ""

if not inputs_path:
    inputs_path = os.getenv("INPUTS")  # From your env.sh file
    if inputs_path is None:
        raise ValueError("The INPUTS environment variable is not set.")

inputs_file = os.path.join(inputs_path, "inputs.yml")
if not os.path.isfile(inputs_file):
    raise FileNotFoundError(
        f"The 'inputs.yml' file could not be found at {inputs_file}."
    )

with open(inputs_file, "r") as file:
    data = yaml.safe_load(file)

strain = data["strain"]
species = data["species"]
chem_spider_key = data["chem_spider_key"]
cs = ChemSpider(chem_spider_key)

# Input and output file paths
sbml_model = os.path.join(inputs_path, data["sbml_model"])
output_file_path = os.path.join(inputs_path, data["output_file_path"])
os.makedirs(output_file_path, exist_ok=True)
if data["protein_file_path"]:
    protein_file_path = os.path.join(inputs_path, data["protein_file_path"])
else:
    protein_file_path = None

cofactors = data["cofactors"]

logging.getLogger("cobra").setLevel(logging.ERROR)
model = cobra.io.read_sbml_model(sbml_model)

# Define a locks for accessing data structures
batch_updates_lock = threading.Lock()
processed_values_lock = threading.Lock()


def batch_loop_setup(file_to_update, df_columns, column_key):
    # Check if checkpointing activated
    if os.path.exists(file_to_update):
        # Read in data, and update list of already processed
        values_df = pd.read_csv(file_to_update)
        processed_values = set(values_df[column_key].unique())
    else:
        # Set up empty data frame and empty set of already processed
        values_df = pd.DataFrame(columns=df_columns)
        processed_values = set()

    # Process in batches
    batch_updates = []
    batch_size = 200
    return values_df, processed_values, batch_updates, batch_size


"""Reusable function to run in parallel for various IO tasks"""


def process_futures(
    futures,
    values_df,
    processed_values,
    column_key,
    batch_updates,
    batch_size,
    file_to_update,
):
    try:
        for future in as_completed(futures):
            result = future.result()
            if result:
                with batch_updates_lock:
                    batch_updates.append(result)
                with processed_values_lock:
                    processed_values.add(result[column_key])

                # Check if batch size is reached
                if len(batch_updates) >= batch_size:
                    with batch_updates_lock:
                        if len(batch_updates) >= batch_size:
                            # Convert batch updates to DataFrame
                            batch_df = pd.DataFrame(batch_updates)
                            # Append batch updates to the main DataFrame and save
                            values_df = pd.concat(
                                [values_df.astype(batch_df.dtypes), batch_df],
                                ignore_index=True,
                            )
                            values_df.to_csv(file_to_update, index=False)
                            # Clear batch updates after saving
                            batch_updates = []
    except Exception as e:
        print(f"Error processing {result[column_key]}: {e}")

    # After the loop, check if there are any remaining updates not yet saved
    if batch_updates:
        with batch_updates_lock:
            if batch_updates:
                batch_df = pd.DataFrame(batch_updates)
                values_df = pd.concat(
                    [values_df.astype(batch_df.dtypes), batch_df], ignore_index=True
                )
                values_df.to_csv(file_to_update, index=False)


def get_smiles_from_csv_apis(name):
    try:
        # Get direct match if possible
        mask = smiles_db["BiGG_metabolite_name"] == name
        if mask.any():
            result = smiles_db[mask]
        else:
            # Search in 'cmpdsynonym' column using vectorized operations for speed
            mask = smiles_db["Metabolite_aliases"].str.contains(
                rf"(?:^|\|){re.escape(name)}(?:\||$)", na=False, regex=True
            )
            result = smiles_db[mask]
        # If a match is found, return the 'smiles' for the match
        if not result.empty:
            try:
                keggid = result.iloc[0]["kegg_metabolite_ID"]
                smiles = result.iloc[0]["SMILES"]
                if DEBUG:
                    print(
                        f"DEBUG: SYNONYM keggid: {keggid} name: {name} smile: {smiles}"
                    )
                return smiles
            except Exception as e:
                print(f"Error while processing 'cmpdsynonym' match: {e}")
    except Exception as e:
        print(f"Error while searching 'cmpdsynonym': {e}")

    try:
        # Looking for corresponding metabolite smiles using the PubChem API
        compounds = pcp.get_compounds(name, "name")
        if len(compounds) > 0:
            try:
                smiles = compounds[0].isomeric_smiles
                if DEBUG:
                    print(f"DEBUG: API name: {name} smile: {smiles}")
                return smiles
            except Exception as e:
                print(f"Error while processing PubChem API response: {e}")
    except Exception as e:
        print(f"Error while querying PubChem API: {e}")

    # try:
    #     # Looking for corresponding metabolite smiles using the ChemSpider API
    #     simple_name = remove_characters_within_brackets(name)
    #     results = cs.search(simple_name)
    #     if results:
    #         try:
    #             smiles = results[0].smiles
    #             if DEBUG:
    #                 print(f"SPIDER name: {name} smile: {smiles}")
    #             return smiles
    #         except Exception as e:
    #             print(f"Error while processing ChemSpider API response: {e}")
    #     print(f"SMILE NOT FOUND FOR name: {name} simple name: {simple_name}")
    # except Exception as e:
    #     print(f"Error while querying ChemSpider API: {e}")

    # If no match is found
    return "Compound not found"


def remove_characters_within_brackets(text):
    # Pattern to match content within brackets (including the brackets themselves)
    pattern = r"\s*\([^)]*\)"
    # Replace matched content with an empty string
    cleaned_text = re.sub(pattern, "", text)
    return cleaned_text


def get_accession(query, target_id):
    url = "https://rest.uniprot.org/uniprotkb/stream"
    params = {
        "query": query,
        "format": "json",
        "fields": "accession,ec,mass,gene_names,lineage,organism_name,sequence",
    }

    try:
        uniprot_response = requests.get(url, params=params).json()
        if not uniprot_response.get("results"):
            return None, None, None, None

        for result in uniprot_response["results"]:
            if "genes" not in result:
                continue

            # Check if this result contains our target gene
            found = False
            for gene in result["genes"]:
                # Check main gene name
                if "geneName" in gene and gene["geneName"]["value"] == target_id:
                    found = True
                    break

                # Check synonyms
                if "synonyms" in gene:
                    for synonym in gene["synonyms"]:
                        if synonym["value"] == target_id:
                            found = True
                            break
                if found:
                    break

            if not found:
                continue

            # If we found the target gene, extract all needed data from this result
            accession = result["primaryAccession"]

            # Get molecular weight and sequence
            mass = result["sequence"]["molWeight"] if "sequence" in result else None
            seq = result["sequence"]["value"] if "sequence" in result else None

            # Extract EC number
            ec = None
            if "proteinDescription" in result:
                protein_desc = result["proteinDescription"]
                if (
                    "recommendedName" in protein_desc
                    and "ecNumbers" in protein_desc["recommendedName"]
                ):
                    ec_numbers = protein_desc["recommendedName"]["ecNumbers"]
                    if ec_numbers:
                        ec = ec_numbers[0]["value"]
                elif (
                    "includes" in protein_desc
                    and "recommendedName" in protein_desc["includes"]
                ):
                    includes_rec = protein_desc["includes"]["recommendedName"]
                    if "ecNumbers" in includes_rec and includes_rec["ecNumbers"]:
                        ec = includes_rec["ecNumbers"][0]["value"]

            return accession, mass, ec, seq

    except Exception as e:
        print(f"Error processing response: {e}")

    return None, None, None, None


# extra = ['umpH', 'ldtB', 'ldtD', 'ldtC', 'pfo', 'glsB',
#          'ldtE', 'ldtA', 'wbbH', 'wzxB', 'umpG', 'fau',
#          'gpmM'
#         ]


def process_uniprot_gene(gene):
    gene_name = gene.name

    pattern = r"^G_.*_\d+$"

    if re.match(pattern, gene_name):
        # Remove the 'G_' prefix and replace '_<integer>' with '.<integer>'
        gene_name = re.sub(r"^G_(.*)_(\d+)$", r"\1.\2", gene_name)

    if gene_name not in processed_genes:
        try:
            # Try with strain first
            query = f'(organism_name:"{species}" AND strain:"{strain}") AND {gene_name}'
            organism = strain
            accession, mass, ec, seq = get_accession(query, gene_name)
            if seq is None:
                print(f"Sequence for strain {organism} and {gene_name} is none")
                # Then try with species
                query = f'(organism_name:"{species}") AND {gene_name}'
                organism = species
                accession, mass, ec, seq = get_accession(query, gene_name)
                if seq is None:
                    print(f"Sequence for strain {organism} and {gene_name} is none")
                    return None
            return {
                "Gene ID": gene.id,
                "Gene name": gene.name,
                "Accession": accession,
                "Sequence": seq,
                "Mass": mass,
                "EC number": ec,
                "Organism": organism,
                "Gene reactions": [r.id for r in gene.reactions],
            }
        except Exception as e:
            print(f"Error processing {gene_name}: {e}")


def process_metabolite_model(m):
    if m.id not in processed_metabolites:
        name = m.name
        if m.formula is not None and m.formula in name:
            name = name.replace(m.formula, "")
        try:
            # Looking for corresponding metabolite smiles in PubChem CSV / API and ChemSpider
            smiles = get_smiles_from_csv_apis(name)
            if DEBUG:
                print(f"DEBUG: type(smiles)={type(smiles)}, smiles={smiles}")
            # Need to cater for different returns
            smiles = "Compound not found" if pd.isna(smiles) else str(smiles).strip()
            if DEBUG:
                print(f"DEBUG: name {name} smiles {smiles}")

            # Collect data for batch update
            return {"metabolite_id": m.id, "name": m.name, "smiles": smiles}

        except Exception as e:
            print(f"Error processing {m.id}: {e}")
            # Handle error (e.g., log it, attempt to recover, skip this metabolite, etc.)


if protein_file_path:
    """Gene-sequence retrieval from fasta file"""
    file_to_update = os.path.join(output_file_path, "gene_sequence_data.csv")
    columns: List[str] = [
        "Gene ID",
        "Gene name",
        "Accession",
        "Sequence",
        "Mass",
        "EC number",
        "Organism",
        "Gene reactions",
    ]
    genes_df = pd.DataFrame(columns)
    records = []
    for record in SeqIO.parse(protein_file_path, "fasta"):
        records.append(record)
    genes = []
    for gene in model.genes:
        gene.id = gene.id.replace("_", ".")
        genes.append(gene)
    intersections = [
        (record, gene)
        for record in records
        if any(record.id == gene.id for gene in genes)
    ]

    data = []
    for r_g_pair in intersections:
        record = r_g_pair[0]
        gene = r_g_pair[1]
        sequence = str(record.seq)
        mass = molecular_weight(sequence, seq_type="protein")
        data.append(
            {
                "Gene ID": record.id,
                "Gene name": record.id,
                "Accession": None,
                "Sequence": sequence,
                "Mass": mass,
                "EC number": None,
                "Organism": None,
                "Gene reactions": [r.id for r in gene.reactions],
            }
        )
    genes_df = pd.DataFrame(data)
    genes_df.to_csv(file_to_update, index=False)

else:
    """Gene-sequence retrieval from UniProt"""
    # Initialise checkpointing
    file_to_update = os.path.join(output_file_path, "gene_sequence_data.csv")
    df_columns = [
        "Gene ID",
        "Gene name",
        "Accession",
        "Sequence",
        "Mass",
        "EC number",
        "Organism",
        "Gene reactions",
    ]
    column_key = "Gene name"
    genes_df, processed_genes, batch_updates, batch_size = batch_loop_setup(
        file_to_update, df_columns, column_key
    )
    batch_size = 100
    num_cpus = min(os.cpu_count(), 16)
    with ThreadPoolExecutor(max_workers=num_cpus) as executor:
        futures = {
            executor.submit(process_uniprot_gene, gene): gene for gene in model.genes
        }
        process_futures(
            futures,
            genes_df,
            processed_genes,
            column_key,
            batch_updates,
            batch_size,
            file_to_update,
        )


"""Metabolite-SMILES retrieval"""
for _ in range(1):  # Loop runs twice (0 and 1)
    # Paths
    file_to_update = os.path.join(output_file_path, "metabolite_smiles_data.csv")
    df_columns = ["metabolite_id", "name", "smiles"]
    column_key = "metabolite_id"
    metabolites_df, processed_metabolites, batch_updates, batch_size = batch_loop_setup(
        file_to_update, df_columns, column_key
    )
    API_counter = [0]
    API_counter[0] = 0
    API_test = "10-Formyltetrahydrofolate"

    batch_size = 100
    with ThreadPoolExecutor(max_workers=1) as executor:
        futures = {
            executor.submit(process_metabolite_model, m): m for m in model.metabolites
        }
        process_futures(
            futures,
            metabolites_df,
            processed_metabolites,
            column_key,
            batch_updates,
            batch_size,
            file_to_update,
        )


"""
Pair sequences with adequate substrate SMILES
"""

# metabolites that should not be considered main substrates of reactions
metabolite_smiles_path = os.path.join(output_file_path, "metabolite_smiles_data.csv")
gene_sequence_path = os.path.join(output_file_path, "gene_sequence_data.csv")

# Load data and ensure files are readable
try:
    metabolites_df = pd.read_csv(metabolite_smiles_path, index_col="metabolite_id")
    genes_df = pd.read_csv(gene_sequence_path, index_col="Gene ID")
except Exception as e:
    print(f"Error loading files: {e}")
    raise

seqs_smiles = []
missing_gene_ids = []
missing_metabolite_ids = []


for gene in model.genes:
    if gene.id == "spontaneous":
        continue
    gene_id = gene.id

    # Regular expression pattern
    pattern = r"^G_.*_\d+$"

    # Check if the variable matches the pattern
    if re.match(pattern, gene_id):
        # Remove the 'G_' prefix and replace '_<integer>' with '.<integer>'
        gene_id = re.sub(r"^G_(.*)_(\d+)$", r"\1.\2", gene_id)

    try:
        sequence = genes_df.loc[gene_id, "Sequence"]
        mass = genes_df.loc[gene_id, "Mass"]
        ec = genes_df.loc[gene_id, "EC number"]
    except KeyError:
        print(f"Gene ID {gene_id} is missing in gene_sequence_data.csv")
        missing_gene_ids.append(gene_id)
        continue

    for r in gene.reactions:
        reactants = [(m.name, m.id) for m in r.reactants if m.name not in cofactors]
        for i in reactants:
            try:
                smiles = metabolites_df.loc[i[1], "smiles"]
            except KeyError:
                print(f"Metabolite ID {i[1]} is missing in metabolite_smiles_data.csv")
                missing_metabolite_ids.append(i[1])
                continue

            seqs_smiles.append(
                {
                    "Gene ID": gene_id,
                    "Gene name": gene.name,
                    "Sequence": sequence,
                    "Reaction ID": r.id,
                    "Reaction name": r.name,
                    "Reaction": r.reaction,
                    "Direction": "Forward",
                    "Substrate Name": i[0],
                    "Substrate ID": i[1],
                    "Substrate Smiles": smiles,
                    "Kcat": "",
                }
            )

        if r.reversibility:
            products = [(m.name, m.id) for m in r.products if m.name not in cofactors]
            for i in products:
                try:
                    smiles = metabolites_df.loc[i[1], "smiles"]
                except KeyError:
                    print(
                        f"Metabolite ID {i[1]} is missing in metabolite_smiles_data.csv"
                    )
                    missing_metabolite_ids.append(i[1])
                    continue

                seqs_smiles.append(
                    {
                        "Gene ID": gene_id,
                        "Gene name": gene.name,
                        "Sequence": sequence,
                        "Reaction ID": r.id,
                        "Reaction name": r.name,
                        "Reaction": r.reaction,
                        "Direction": "Reverse",
                        "Substrate Name": i[0],
                        "Substrate ID": i[1],
                        "Substrate Smiles": smiles,
                        "Kcat": "",
                    }
                )

# Convert to DataFrame and save to CSV
seqs_smiles_df = pd.DataFrame(seqs_smiles)
seqs_smiles_df.to_csv(
    os.path.join(output_file_path, "sequences_smiles.csv"), index=False
)

# Log missing gene and metabolite IDs
if missing_gene_ids:
    print(
        f"Warning: The following gene IDs were not found in gene_sequence_data.csv: {set(missing_gene_ids)}"
    )

if missing_metabolite_ids:
    print(
        f"Warning: The following metabolite IDs were not found in metabolite_smiles_data.csv: {set(missing_metabolite_ids)}"
    )
