
from switch_generator import SwitchGenerator
from prokaryotic_switch_generator import ProkaryoticSwitchGenerator
from window_folding_based_selection import get_gene_top_ranked_windows
from utils.send_mail import send_email_with_attachment as send
#from server import bucket

import json
import pickle
import logging
import pandas as pd
import io
import numpy as np
import pandas as pd
from fuzzysearch import find_near_matches
from joblib import Parallel, delayed
import time
import sys
from RNA import RNA
from nupack import *
config.threads = 8

EDIT_DIST = 7
TRIGGERS_BATCH = 15
SWITCH_BATCH = 2

E_COLI_DATA_PATH = "Escherichia_coli_ASM886v2.pkl"
HUMAN_DATA_PATH = "Homo_sapiens_GRCh38.pkl"
YEAST_DATA_PATH = "Saccharomyces_cerevisiae_S288C.pkl"
"""
DATA_PATHS = {
    "E.coli": E_COLI_DATA_PATH,
    "Homo sapiens": HUMAN_DATA_PATH,
    "Saccharomyces cerevisiae": YEAST_DATA_PATH
}
"""
"""
# for development
E_COLI_DATA_PATH = "/Users/netanelerlich/PycharmProjects/software/data/Escherichia_coli_ASM886v2.pkl"
HUMAN_DATA_PATH = "/Users/netanelerlich/PycharmProjects/software/data/Homo_sapiens_GRCh38.pkl"
YEAST_DATA_PATH = "/Users/netanelerlich/PycharmProjects/software/data/Saccharomyces_cerevisiae_S288C.pkl"

DATA_PATHS = {
    "E.coli": E_COLI_DATA_PATH,
    "Homo sapiens": HUMAN_DATA_PATH,
    "Saccharomyces cerevisiae": YEAST_DATA_PATH
    }
"""


logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

def extract_top_homology_sequences(triggers_homology_mapping):
    print(triggers_homology_mapping)
    if not isinstance(triggers_homology_mapping, list):
        raise ValueError("Expected a list for 'triggers_homology_mapping'")

    homo_dfs = []
    for trigger_homology in triggers_homology_mapping:
        if trigger_homology:
            trigger_all_matches = []
            for gene_homology in trigger_homology:
                # Check if gene_homology is a valid dictionary
                if isinstance(gene_homology, dict):
                    trigger_all_matches.extend(gene_homology.values())
                else:
                    continue

            if not trigger_all_matches or not trigger_all_matches[0]:
                # If no matches are found, create an empty DataFrame
                trig_res_df = pd.DataFrame(
                    columns=['distance', 'idx', 'sequence', 'gene', 'protein', 'homologous_trigger_mfe'])
            else:
                matches_df = pd.DataFrame(trigger_all_matches[0])

                # Ensure mfe calculation is safe
                mfe_dict = {'homologous_trigger_mfe': []}
                for homos in trigger_all_matches[0]:
                    if 'sequence' in homos:
                        homo_trigger = homos['sequence']
                        structure, mfe = RNA.fold(homo_trigger)
                        mfe_dict['homologous_trigger_mfe'].append(mfe)
                    else:
                        mfe_dict['homologous_trigger_mfe'].append(None)

                mfe_df = pd.DataFrame(mfe_dict)
                trig_res_df = pd.concat([matches_df, mfe_df], axis=1)
        else:
            trig_res_df = pd.DataFrame(
                columns=['distance', 'idx', 'sequence', 'gene', 'protein', 'homologous_trigger_mfe']
            )

        homo_dfs.append(trig_res_df)

    # RRF calculation with valid checks
    higher = ['homologous_trigger_mfe']
    lower = ['distance']
    print(homo_dfs)
    if not homo_dfs:
        return []

    rrf_rank = [RRF(trigger_homology_df, higher, lower, index='sequence') for trigger_homology_df in homo_dfs]
    return rrf_rank


def route_input(email, target_seq, trigger, reporter_gene, cell_type, user_trigger_boo, transcripts_list):
    # Basic input validation
    if not isinstance(email, str) or "@" not in email:
        raise ValueError("Invalid email provided")
    if not isinstance(target_seq, str) or not target_seq:
        raise ValueError("Invalid target sequence provided")

    results = pd.DataFrame()
    #blob_name = DATA_PATHS.get(cell_type)
    #blob = bucket.blob(blob_name)
    #bytes_data = blob.download_as_bytes()
    #file_like_obj = io.BytesIO(bytes_data)
    #cell_type_transcripts = pickle.load(file_like_obj)

    #Test
    cell_type_transcripts = [{'gene':'test', 'protein':'prot',
                              'sequence':"AGATAGATAGATAGATAGTAAGAATGAATGAATAGAATAGAATGATAAGATAGATAGATAGATAGTAAGAATGAATGAATAGAATAGAATGATAAGATAGATAGATAGATAGTAAGAATGAATGAATAGAATAGAATGATA"}]
    # Update transcripts_list
    if transcripts_list != 'EMPTY':
        transcripts_list = json.loads(transcripts_list)
        transcripts_list.extend(cell_type_transcripts)
    else:
        transcripts_list = cell_type_transcripts

    # Get optional triggers if user hasn't provided one
    triggers_with_mfe = np.array([[trigger, 1]])
    if user_trigger_boo == 'EMPTY':
        cell_window = 30 if cell_type == 'E.coli' else 23
        optional_triggers_df = get_gene_top_ranked_windows(target_seq, window_size=cell_window)
        n_triggers = min(TRIGGERS_BATCH, len(optional_triggers_df.index))
        triggers_with_mfe = np.array(optional_triggers_df.iloc[:n_triggers, :][["window_sequence", "mfe_score"]])

    results[['trigger_window', 'mfe_score']] = triggers_with_mfe
    triggers_seqs = triggers_with_mfe[:, 0]

    # Homology search with error handling
    try:
        s_h = time.time()
        homo_res = Parallel(n_jobs=-1)(delayed(find_homology)(trigger, transcripts_list) for trigger in triggers_seqs)
        e_h = time.time()
        print(f'Homology search time= {e_h - s_h} seconds, n_triggers = {len(triggers_seqs)}, cell= {cell_type}')
    except Exception as e:
        raise RuntimeError(f"Error during homology search: {e}")

    # Extract top homology sequences
    try:
        rrf_ranks = extract_top_homology_sequences(homo_res)
    except Exception as e:
        raise RuntimeError(f"Error extracting top homology sequences: {e}")

    # Filter triggers
    n_switch = min(SWITCH_BATCH, 1)

    homology_sequences = [ranked_df['sequence'].get(0) for ranked_df in rrf_ranks][:n_switch]
    triggers_seqs = triggers_seqs[:n_switch]

    # Generate switches safely
    try:
        s_switch = time.time()
        switch_res = Parallel(n_jobs=4)(delayed(generate_switch)(f_trigger, f_top_homology_sequence, reporter_gene, cell_type)
                                        for f_trigger, f_top_homology_sequence in zip(triggers_seqs, homology_sequences))
        e_switch = time.time()
        print(f'Switch generation time= {e_switch - s_switch} seconds, for {SWITCH_BATCH} switches')
    except Exception as e:
        raise RuntimeError(f"Error generating switches: {e}")
    print(switch_res)
    results[['switch',  "complex_concentration"]] = pd.DataFrame(switch_res, columns=['switch', "complex_concentration"])

    # Prepare and send the report
    prepare_and_send_report(results, rrf_ranks, email)
    return


def RRF(ranking_df, higher_is_better_cols, lower_is_better_cols, index, k=60):
    if not isinstance(higher_is_better_cols, list):
        raise ValueError("'higher_is_better_cols' should be a list")
    if not isinstance(lower_is_better_cols, list):
        raise ValueError("'lower_is_better_cols' should be a list")
    if not isinstance(ranking_df, pd.DataFrame):
        raise ValueError("'ranking_df' should be a pandas DataFrame")
    if not isinstance(index, str):
        raise ValueError("'index' should be a string")

    ranking_df = ranking_df.copy().reset_index(drop=True)

    for col in higher_is_better_cols:
        if col not in ranking_df.columns:
            raise KeyError(f"Column '{col}' not found in DataFrame")
        ranking_df[col + '_rank'] = ranking_df[col].rank(ascending=False)

    for col in lower_is_better_cols:
        if col not in ranking_df.columns:
            raise KeyError(f"Column '{col}' not found in DataFrame")
        ranking_df[col + '_rank'] = ranking_df[col].rank(ascending=True)

    ranked_columns = higher_is_better_cols + lower_is_better_cols
    ranked_columns = [f'{col}_rank' for col in ranked_columns]

    ranking_df[f'{index}_RRF'] = ranking_df[ranked_columns].apply(
        lambda row: sum(1 / (k + rank) for rank in row if pd.notna(rank)), axis=1)

    return ranking_df.sort_values(by=f'{index}_RRF', ascending=True)


def build_homology_map(trigger, seq, gene_name, protein_name):
    # Validate inputs
    if not isinstance(trigger, str) or not trigger:
        logging.error("Invalid trigger sequence provided.")
        return []
    if not isinstance(seq, str) or not seq:
        logging.error("Invalid sequence provided.")
        return []
    if not isinstance(gene_name, str):
        logging.error("Gene name must be a string.")
        gene_name = "Unknown"
    if not isinstance(protein_name, str):
        logging.error("Protein name must be a string.")
        protein_name = "Unknown"

    sequence_match_mapping = []

    try:
        # Find matches with error handling
        matches = find_near_matches(trigger, seq, max_insertions=0, max_deletions=0, max_l_dist=EDIT_DIST)
    except Exception as e:
        logging.error(f"Error finding matches: {e}")
        return sequence_match_mapping

    if not matches:
        logging.info("No matches found for the trigger in the sequence.")
        return sequence_match_mapping

    for match_obj in matches:
        locus_start = match_obj.start
        locus_end = match_obj.end
        sub_seq = seq[locus_start:locus_end]
        seq_location = (locus_start, locus_end)
        distance = match_obj.dist
        sequence_match_mapping.append({
            'distance': distance,
            'idx': seq_location,
            'sequence': sub_seq,
            'gene': gene_name,
            'protein': protein_name
        })

    return sequence_match_mapping


def find_homology(sequence, genome_data):
    # Validate inputs
    if not isinstance(sequence, str) or not sequence:
        logging.error("Invalid sequence provided for homology search.")
        return []
    if not isinstance(genome_data, list):
        logging.error("Invalid genome data provided.")
        return []

    genes_sub_sequences = []
    for gene_data_dict in genome_data:
        # Validate the structure of genome data
        if not all(key in gene_data_dict for key in ['sequence', 'gene', 'protein']):
            logging.error("Missing keys in gene data: 'sequence', 'gene', 'protein'")
            continue

        gene_sub_seqs = {}
        mRNA_seq = gene_data_dict.get('sequence', "")
        gene_name = gene_data_dict.get('gene', "")
        protein = gene_data_dict.get('protein', "")

        if not mRNA_seq:
            logging.warning(f"Skipping gene {gene_name} due to missing mRNA sequence.")
            continue

        seq_match_mapping = build_homology_map(sequence, mRNA_seq, gene_name, protein)
        if seq_match_mapping:
            gene_sub_seqs[gene_name] = seq_match_mapping
            genes_sub_sequences.append(gene_sub_seqs)

    return genes_sub_sequences


def generate_switch(trigger, homologous_sequence, reporter_gene, cell_type):
    # Validate inputs
    if not isinstance(trigger, str) or not trigger:
        logging.error("Invalid trigger sequence provided.")
        return None, None
    print(homologous_sequence)

    if not isinstance(reporter_gene, str) or not reporter_gene:
        logging.error("Invalid reporter gene provided.")
        return None, None

    try:
        if cell_type == "E.coli":
            switch_generator = ProkaryoticSwitchGenerator(translated_gene_sequence=reporter_gene)
        else:
            switch_generator = SwitchGenerator(translated_gene_sequence=reporter_gene)

        switch_designed_strand, ensemble_defect, complex_concentration, switch_mfe_structure = switch_generator.get_switch(
            trigger_sequence=trigger, healthy_sequence=homologous_sequence
        )

    except Exception as e:
        logging.error(f"Error generating switch: {e}")
        return None, None

    return switch_designed_strand, complex_concentration


def prepare_and_send_report(df_results, rrf_ranks, email):
    # Validate inputs
    if not isinstance(df_results, pd.DataFrame):
        logging.error("df_results must be a pandas DataFrame.")
        return
    if not isinstance(rrf_ranks, list) or not all(isinstance(item, pd.DataFrame) for item in rrf_ranks):
        logging.error("rrf_ranks must be a list of pandas DataFrames.")
        return
    if not isinstance(email, str) or "@" not in email:
        logging.error("Invalid email address provided.")
        return

    rename_dict = {
        "trigger_window": "Trigger",
        "switch": "Switch",
        "sequence": "Competitor Sequence",
        "mfe_score": "Trigger mfe Score",
        # Add regressor results here
        "sequence_RRF": "Competition Score",
        "distance": "Substitution Distance",
        "idx": "Competitor Location",
        "gene": "Competitor Gene",
        "protein": "Competitor Protein"
    }

    try:
        cols_rearrange = list(rename_dict.values())
        combined_rows = []
        rrf_ranks = rrf_ranks.copy()

        for index, row in df_results.iterrows():
            combined_df = pd.concat([pd.DataFrame([row]), rrf_ranks[index]], ignore_index=True)
            combined_rows.append(combined_df)

        final_df = pd.concat(combined_rows, ignore_index=True).rename(columns=rename_dict)
        final_df = final_df[cols_rearrange]

        # Send the report via the send() function
        send(final_df, email)
    except Exception as e:
        logging.error(f"Error preparing or sending report: {e}")


if __name__ == '__main__':
    # Get the arguments from the user form.
    """
    s_mail = sys.argv[1]
    s_target_seq = sys.argv[2]
    s_trigger = sys.argv[3]
    s_reporter_gene = sys.argv[4]
    s_cell_type = sys.argv[5]
    s_user_trigger_boo = sys.argv[6]
    s_transcripts_list = sys.argv[7]
    route_input(s_mail, s_target_seq, s_trigger, s_reporter_gene, s_cell_type, s_user_trigger_boo, s_transcripts_list)
    """

    route_input('erlichnet57@gmail.com', 'ATGCGTACGATGCGTACGATGCGTACGATGCGTACGATGCGTACGATGCGTACGATGCGTACGATGCGTACG',
                'EMPTY', 'ATGCGTACGATGCGTACGATGCGTACGATGCGTACGAGATGAATGATA',
                'E.coli', 'True', 'EMPTY')



