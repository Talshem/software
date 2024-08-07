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
from viennaRNA import rnaUp

WINDOW_SIZE = 23
CONTEXT_WINDOW_SIZE = 50
WINDOWS_COUNT = 1000

def run_rna_up(sequence1, sequence2):
    try:
        return subprocess.check_output(
            ["RNAup", "-d", "2", "--noLP", "--include_both", "-c", "'S'"],
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


def get_potential_windows_scores(mRNA: str, window_size: int = WINDOW_SIZE, context_window_size: int = CONTEXT_WINDOW_SIZE) -> pd.DataFrame:
    windows_with_scores = {}

    for i in tqdm(range(len(mRNA))):
        if i + window_size > len(mRNA):
            break

        window_sequence = mRNA[i:i + window_size]
        reverse_complement = Seq(window_sequence).reverse_complement_rna()
        mfe = get_mRNA_opening_mfe(reverse_complement,
                                   mRNA[max(i - context_window_size, 0):i + window_size + context_window_size])
        windows_with_scores[window_sequence] = mfe
    return windows_with_scores

if __name__ == '__main__':
    # gene_id = sys.argv[1]
    # sorted_windows = get_gene_top_ranked_windows(gene_id)
    # sorted_windows.head(WINDOWS_COUNT).to_excel(
    #     f"/tamir2/peba/CancerToeholdPipeline2.0/igem_highly_expressed_genes_windows/{gene_id}_potential_windows.xlsx", index=False)

    # GFP for eukaryotes
    mRNA = "ATGAGAAAAGGCGAAGAATTATTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGTCACAAATTTTCTGTCTCCGGTGAAGGTGAAGGTGATGCTACTTACGGTAAATTGACCTTAAAATTTATTTGTACTACTGGTAAATTGCCAGTTCCATGGCCAACCTTAGTCACTACTTTCGGTTATGGTGTTCAATGTTTTGCTAGATACCCAGATCATATGAAACAACATGACTTTTTCAAGTCTGCCATGCCAGAAGGTTATGTTCAAGAAAGAACTATTTTTTTCAAAGATGACGGTAACTACAAGACCAGAGCTGAAGTCAAGTTTGAAGGTGATACCTTAGTTAATAGAATCGAATTAAAAGGTATTGATTTTAAAGAAGATGGTAACATTTTAGGTCACAAATTGGAATACAACTATAACTCTCACAATGTTTACATCATGGCTGACAAACAAAAGAATGGTATCAAAGTTAACTTCAAAATTAGACACAACATTGAAGATGGTTCTGTTCAATTAGCTGACCATTATCAACAAAATACTCCAATTGGTGATGGTCCAGTCTTGTTACCAGACAACCATTACTTATCCACTCAATCTGCCTTATCCAAAGATCCAAACGAAAAGAGAGACCACATGGTCTTGTTAGAATTTGTTACTGCTGCTGGTATTACCCATGGTATGGATGAATTGTACAAA".replace("T", "U")
    potential_windows_scores_df = get_potential_windows_scores(mRNA, WINDOW_SIZE, CONTEXT_WINDOW_SIZE)
    potential_windows_scores_df.to_excel(
        f"/tamir2/peba/CancerToeholdPipeline2.0/igem_gfp_trigger_library/potential_windows_fixed_mfe.xlsx")
