import json
import os
import re
import subprocess
from typing import Dict, Any, List
import pandas as pd
import requests
from Bio import SeqIO
from Bio.Seq import Seq
from tqdm import tqdm

MRNAS_FOLDER = ""
CHROMOSOMES_FOLDER = ""

#RNAUP_BINARY_NAME = "/Users/netanelerlich/ViennaRNA-2.6.4/src/bin/RNAup"
#WINDOW_SIZE = 23

RNAUP_BINARY_NAME = "/workspace/miniconda3/envs/myenv/bin/RNAup"
CONTEXT_WINDOW_SIZE = 50
WINDOWS_COUNT = 1000

def get_gene_name(gene_id: str) -> str:
    response = requests.post("https://biit.cs.ut.ee/gprofiler/api/convert/convert", json={
        "organism": "hsapiens",
        "query": gene_id,
        "aresolve": {},
        "target": "ENSG",
        "numeric_namespace": "GENECARDS_ACC",
        "output": "app_json"
    })
    result = response.json()["results"][0]
    gene_name = result["name"]
    return gene_name


def get_gene_transcript_ids(gene_id: str) -> List[str]:
    response = requests.get(f"https://mygene.info/v3/gene/{gene_id}")
    ensemble_data = response.json()["ensembl"]

    if isinstance(ensemble_data, dict):
        transcript_ids = ensemble_data["transcript"]
    else:
        # type is a list of dicts
        transcript_ids = [] # default in case no match is found

        for item in ensemble_data:
            if item["gene"] == gene_id:
                transcript_ids = item["transcript"]
                break

    if isinstance(transcript_ids, list):
        return transcript_ids

    return [transcript_ids]


def get_transcript_sequence(transcript_id: str) -> str:
    response = requests.get(f"https://rest.ensembl.org/sequence/id/{transcript_id}?content-type=application/json")
    return response.json()["seq"].replace("T", "U")


def get_gene_mRNA_info(gene_name) -> Dict[str, Any]:
    mrna_info_file = open(os.path.join(MRNAS_FOLDER, f"{gene_name}_annotation.json"))
    return json.load(mrna_info_file)


def get_chromosome(gene_chromosome_name) -> str:
    chromosome_file_name = os.path.join(CHROMOSOMES_FOLDER,
                                        f"Homo_sapiens.GRCh38.dna_sm.chromosome.{gene_chromosome_name}.fa")
    return list(SeqIO.parse(open(chromosome_file_name), 'fasta'))[0].seq


def run_rna_up(sequence1, sequence2):
    try:
        return subprocess.check_output(
            [RNAUP_BINARY_NAME, "-d", "2", "--noLP", "--include_both", "-c" "'S'"],
            universal_newlines=True,
            input=f"{sequence1}&{sequence2}",
            text=True
        )
    except Exception as e:
        print(f"Error while running RNAup on {sequence1} and {sequence2}. Error: {e}")
        raise


def get_mRNA_opening_mfe(trigger, mRNA) -> float:
    rna_up_output = run_rna_up(trigger, mRNA)

    try:
        regex_result = re.search("\([0-9-.]+ = [0-9-.]+ \+ [0-9-+.]+ \+ ([0-9-+.]+\)|inf)",
                                 rna_up_output).group()
    except Exception as e:
        print(trigger)
        raise e

    string_values = regex_result.replace("(", "").replace(")", "").replace("=", "").replace("+", "").split()
    values = [float(value) for value in string_values]
    mRNA_opening_energy = values[2]
    return mRNA_opening_energy


def get_potential_windows_scores(mRNA: str, window_size: int, context_window_size: int) -> pd.DataFrame:
    windows_with_scores = []

    for i in tqdm(range(len(mRNA))):
        if i + window_size > len(mRNA):
            break

        window_sequence = mRNA[i:i + window_size]
        reverse_complement = Seq(window_sequence).reverse_complement_rna()
        mfe = get_mRNA_opening_mfe(reverse_complement,
                                   mRNA[max(i - context_window_size, 0):i + window_size + context_window_size])
        windows_with_scores.append((window_sequence, i, mfe))

    return pd.DataFrame(windows_with_scores, columns=["window_sequence", "window_index", 'mfe_score'])


def get_gene_top_ranked_windows(mRNA: str, window_size: int = 23) -> pd.DataFrame:

    windows_df = pd.DataFrame(columns=["window_sequence", "window_index", "transcript_id", "mfe_score"])
    potential_windows_scores_df = get_potential_windows_scores(mRNA, window_size, CONTEXT_WINDOW_SIZE)
    windows_df = pd.concat([
        windows_df, potential_windows_scores_df.sort_values(by="mfe_score", ascending=True)
        ])

    return windows_df


if __name__ == '__main__':
    # gene_id = sys.argv[1]
    # sorted_windows = get_gene_top_ranked_windows(gene_id)
    # sorted_windows.head(WINDOWS_COUNT).to_excel(
    #     f"/tamir2/peba/CancerToeholdPipeline2.0/igem_highly_expressed_genes_windows/{gene_id}_potential_windows.xlsx", index=False)
    """
    # GFP for eukaryotes
    mRNA = "ATGAGAAAAGGCGAAGAATTATTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGTCACAAATTTTCTGTCTCCGGTGAAGGTGAAGGTGATGCTACTTACGGTAAATTGACCTTAAAATTTATTTGTACTACTGGTAAATTGCCAGTTCCATGGCCAACCTTAGTCACTACTTTCGGTTATGGTGTTCAATGTTTTGCTAGATACCCAGATCATATGAAACAACATGACTTTTTCAAGTCTGCCATGCCAGAAGGTTATGTTCAAGAAAGAACTATTTTTTTCAAAGATGACGGTAACTACAAGACCAGAGCTGAAGTCAAGTTTGAAGGTGATACCTTAGTTAATAGAATCGAATTAAAAGGTATTGATTTTAAAGAAGATGGTAACATTTTAGGTCACAAATTGGAATACAACTATAACTCTCACAATGTTTACATCATGGCTGACAAACAAAAGAATGGTATCAAAGTTAACTTCAAAATTAGACACAACATTGAAGATGGTTCTGTTCAATTAGCTGACCATTATCAACAAAATACTCCAATTGGTGATGGTCCAGTCTTGTTACCAGACAACCATTACTTATCCACTCAATCTGCCTTATCCAAAGATCCAAACGAAAAGAGAGACCACATGGTCTTGTTAGAATTTGTTACTGCTGCTGGTATTACCCATGGTATGGATGAATTGTACAAA".replace("T", "U")
    potential_windows_scores_df = get_potential_windows_scores(mRNA, WINDOW_SIZE, CONTEXT_WINDOW_SIZE)
    potential_windows_scores_df.to_excel(
        f"/tamir2/peba/CancerToeholdPipeline2.0/igem_gfp_trigger_library/potential_windows_fixed_mfe.xlsx")
    """


