
from switch_generator import SwitchGenerator
from prokaryotic_switch_generator import ProkaryoticSwitchGenerator
from window_folding_based_selection import get_gene_top_ranked_windows
from utils.send_mail import send_email_with_attachment as send
from eukaryotic_score_calculator import EukaryoticScoreCalculator
from server import bucket

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


model_path = "/workspace/tool/files/webtool_model.txt"
feature_path = "/workspace/tool/files/model_features.txt"

E_COLI_DATA_PATH = "Escherichia_coli_ASM886v2.pkl"
HUMAN_DATA_PATH = "Homo_sapiens_GRCh38.pkl"
YEAST_DATA_PATH = "Saccharomyces_cerevisiae_S288C.pkl"


# for development
# model_path = "/Users/netanelerlich/Desktop/IGEM/webtool_model.txt"
# feature_path = "/Users/netanelerlich/Desktop/IGEM/model_features.txt"
# E_COLI_DATA_PATH = "/Users/netanelerlich/PycharmProjects/software/data/Escherichia_coli_ASM886v2.pkl"
# HUMAN_DATA_PATH = "/Users/netanelerlich/PycharmProjects/software/data/Homo_sapiens_GRCh38.pkl"
# YEAST_DATA_PATH = "/Users/netanelerlich/PycharmProjects/software/data/Saccharomyces_cerevisiae_S288C.pkl"


DATA_PATHS = {
    "E.coli": E_COLI_DATA_PATH,
    "Homo sapiens": HUMAN_DATA_PATH,
    "Saccharomyces cerevisiae": YEAST_DATA_PATH
}


logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
def extract_top_homology_sequences(triggers_homology_mapping):
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

                mfe_dict = {'homologous_trigger_mfe': []}
                for homos in trigger_all_matches[0]:
                    if 'sequence' in homos.keys():
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

    # Competition RRF calculation with valid checks
    higher = ['homologous_trigger_mfe']
    lower = ['distance']
    if not homo_dfs:
        return []
    rrf_rank = [RRF(trigger_homology_df, higher, lower, index='sequence') for trigger_homology_df in homo_dfs]
    return rrf_rank


def route_input(email, target_seq, trigger, reporter_gene, cell_type, user_trigger_boo, transcripts_list):
    # Basic input validation
    if not isinstance(email, str) or "@" not in email:
        raise ValueError("Invalid email provided")


    blob_name = DATA_PATHS.get(cell_type)
    blob = bucket.blob(blob_name)
    bytes_data = blob.download_as_bytes()
    file_like_obj = io.BytesIO(bytes_data)
    cell_type_transcripts = pickle.load(file_like_obj)


    # with open(blob_name, 'rb') as f:
    #     cell_type_transcripts = pickle.load(f)


    # Update transcripts_list
    if transcripts_list != 'EMPTY':
        transcripts_list = json.loads(transcripts_list)
        transcripts_list.extend(cell_type_transcripts)
    else:
        transcripts_list = cell_type_transcripts

    # Get optional triggers
    triggers_with_mfe_df = [[trigger, None]]
    n_switch = 1
    if user_trigger_boo == 'EMPTY':
        cell_window = 30 if cell_type == 'E.coli' else 23
        n_switch = SWITCH_BATCH
        n_triggers = TRIGGERS_BATCH

        optional_triggers_df = get_gene_top_ranked_windows(target_seq, window_size=cell_window)
        n_triggers = min(n_triggers, len(optional_triggers_df.index))
        triggers_with_mfe_df = (optional_triggers_df.iloc[:n_triggers, :][["window_sequence", "mfe_score"]].
                                reset_index(drop=True)).rename(columns={"window_sequence": "Trigger",
                                                                        "mfe_score": "Trigger MFE Score"})

    # Extract triggers
    triggers_with_mfe = np.array(triggers_with_mfe_df)
    triggers_sequences = triggers_with_mfe[:, 0]

    # Homology search
    homo_res = Parallel(n_jobs=-1)(delayed(find_homology)(trigger, transcripts_list) for trigger in triggers_sequences)

    # Extract top homology sequences
    rrf_ranks = extract_top_homology_sequences(homo_res)

    # Top (triggers, homology) -> switch design
    homology_sequences_final = [ranked_df['sequence'].iloc[0] if len(ranked_df) > 0 else None for ranked_df in rrf_ranks][:n_switch]
    homology_sequences_final_score = [ranked_df['sequence_RRF'].iloc[0] if len(ranked_df) > 0 else None for ranked_df in rrf_ranks][:n_switch]
    triggers_sequences_final = triggers_sequences[:n_switch]

    # Generate switches safely
    switch_res = Parallel(n_jobs=4)(delayed(generate_switch)(f_trigger, f_top_homology_sequence, reporter_gene, cell_type)
                                    for f_trigger, f_top_homology_sequence in
                                    zip(triggers_sequences_final, homology_sequences_final))
    # Model Score
    regg_results = {'Switch': [], 'Fold Change Toehold Score': []}
    score_calculator = get_score_calc(cell_type)
    switches = [res_tup[0] for res_tup in switch_res]

    for switch_model, trigger_model in zip(switches, triggers_sequences_final):
        score = score_calculator.get_score(switch_model, trigger_model)
        regg_results['Fold Change Toehold Score'].append(score)
        regg_results['Switch'].append(switch_model)
    regg_results_df = pd.DataFrame(regg_results)

    # Organize
    comp_df = pd.DataFrame(homology_sequences_final_score, columns=['Competition Score'])
    results = pd.concat([triggers_with_mfe_df.head(n_switch), regg_results_df, comp_df], axis=1)

    # RRF calculation- fold change toehold score, competition score, trigger mfe score
    higher = ["Fold Change Toehold Score"]
    lower = ["Trigger MFE Score", "Competition Score"]
    scored_tot = RRF(results, higher, lower, index='Switch').rename(columns={"Switch_RRF": "Combined Score"})

    # Prepare and send the report
    final_df = concat_tables(scored_tot, rrf_ranks[:n_switch])
    print(final_df.to_string())
    send(final_df, email)
    return

def get_score_calc(cell_type):
    with open(feature_path, 'r') as file:
        feature_f = file.read()
    feature_list = eval(feature_f)

    if cell_type == 'E.coli':
        prokaryotes = {
            "ideal_stem_structure": '..............((((((((((...((((((...........))))))...))))))))))',
            "trigger_binding_site_start": 3,
            "trigger_binding_site_end": 32,
            "stem_start": 14,
            "stem_end": 62,
            "loop_start": 33,
            "loop_end": 43,
            "stem_top_start": 27,
            "stem_top_end": 49
        }
        args = tuple(prokaryotes.values())
        score_gen = EukaryoticScoreCalculator(model_path, feature_list, *args)
    else:
        eukaryotes = {
            "ideal_stem_structure": '........(((((((((((((((............)))))))))))...',
            "trigger_binding_site_start": 0,
            "trigger_binding_site_end": 22,
            "stem_start": 8,
            "stem_end": 49,
            "loop_start": 23,
            "loop_end": 34,
            "stem_top_start": 17,
            "stem_top_end": 40
        }
        args = tuple(eukaryotes.values())
        score_gen = EukaryoticScoreCalculator(model_path, feature_list, *args)
    return score_gen


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
    cols = list(ranking_df.columns)


    for col in higher_is_better_cols:
        ranking_df[col + '_rank'] = ranking_df[col].rank(ascending=False)

    for col in lower_is_better_cols:
        ranking_df[col + '_rank'] = ranking_df[col].rank(ascending=True)

    ranked_columns = higher_is_better_cols + lower_is_better_cols
    ranked_columns = [f'{col}_rank' for col in ranked_columns]

    ranking_df[f'{index}_RRF'] = ranking_df[ranked_columns].apply(
        lambda row: sum(1 / (k + rank) for rank in row if pd.notna(rank)), axis=1)

    mask = cols + [f'{index}_RRF']
    return ranking_df[mask].sort_values(by=f'{index}_RRF', ascending=False)


def build_homology_map(trigger, seq, gene_name, protein_name):
    # Validate inputs
    if not isinstance(trigger, str) or not trigger:
        logging.error("Invalid trigger sequence provided.")
        return []
    if not isinstance(seq, str) or not seq:
        logging.error("Invalid sequence provided.")
        return []

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


def concat_tables(df_results, rrf_ranks):
    # Validate inputs
    combined_rows = []
    rrf_ranks = rrf_ranks.copy()
    if all([len(rrf_df) == 0 for rrf_df in rrf_ranks]):
        df_results['Competition Score'] = 'Competitor Not Found'
        return df_results


    for index, row in df_results.iterrows():
        combined_df = pd.concat([pd.DataFrame([row]), rrf_ranks[index]], axis=1, ignore_index=True)
        combined_rows.append(combined_df)

    final_df = pd.concat(combined_rows, ignore_index=True).drop(columns=[12])
    rename_dict = {
        0: 'Trigger',
        1: 'Trigger MFE Score',
        2: 'Switch',
        3: 'Fold Change Toehold Score',
        4: 'Competition Score',
        5: 'Combined Score',
        6: 'Substitution Distance',
        7: 'Competitor Location',
        8: 'Competitor Sequence',
        9: 'Competitor Gene',
        10: 'Competitor Protein',
        11: 'Competitor Trigger MFE'
    }

    # Rearrange the columns
    final_df = final_df.rename(columns=rename_dict)
    return final_df



if __name__ == '__main__':
    # Get the arguments from the user form.
    s_mail = sys.argv[1]
    s_target_seq = sys.argv[2]
    s_trigger = sys.argv[3]
    s_reporter_gene = sys.argv[4]
    s_cell_type = sys.argv[5]
    s_user_trigger_boo = sys.argv[6]
    s_transcripts_list = sys.argv[7]
    route_input(s_mail, s_target_seq, s_trigger, s_reporter_gene, s_cell_type, s_user_trigger_boo, s_transcripts_list)
    # for development generate inputs

    """
    s_mail = 'erlichnet57@gmail.com'
    s_target_seq = "TAGAATAGTAAGATAGATAGTAGAATAGTAAGATAGATAGTAGAATAGTAAGATAGATAGTAGAATAGTAAGATAGATAGTAGAATAGTAAGATAGATTAGAATAGTAAGATAGATAGTAGAATAGTAAGATAGATAGTAGAATAGTAAGATAGATAGTAGAATAGTAAGATAGATAGTAGAATAGTAAGATAGAT"
    s_trigger = 'EMPTY'
    s_reporter_gene = "ATGCGTAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGGCACAAATTTTCTGTCAGTGGAGAGGGTGAAGGTGATGCAACATACGGAAAACTTACCCTTAAATTTATTTGCACTACTGGAAAACTACCTGTTCCGTGGCCAACACTTGTCACTACTTTCGGTTATGGTGTTCAATGCTTTGCGAGATACCCAGATCACATGAAACAGCATGACTTTTTCAAGAGTGCCATGCCCGAAGGTTACGTACAGGAAAGAACTATATTTTTCAAAGATGACGGGAACTACAAGACACGTGCTGAAGTCAAGTTTGAAGGTGATACCCTTGTTAATAGAATCGAGTTAAAAGGTATTGATTTTAAAGAAGATGGAAACATTCTTGGACACAAATTGGAATACAACTATAACTCACACAATGTATACATCATGGCAGACAAACAAAAGAATGGAATCAAAGTTAACTTCAAAATTAGACACAACATTGAAGATGGAAGCGTTCAACTAGCAGACCATTATCAACAAAATACTCCGATTGGCGATGGCCCTGTCCTTTTACCAGACAACCATTACCTGTCCACACAATCTGCCCTTTCGAAAGATCCCAACGAAAAGAGAGACCACATGGTCCTTCTTGAGTTTGTAACCGCTGCTGGGATTACACATGGCATGGATGAACTATACAAA".replace('T', 'U')
    s_cell_type = 'E.coli'
    s_user_trigger_boo = 'EMPTY'
    s_transcripts_list = 'EMPTY'
    route_input(s_mail, s_target_seq, s_trigger, s_reporter_gene, s_cell_type, s_user_trigger_boo, s_transcripts_list)
    """

    # with open(YEAST_DATA_PATH,'rb') as f:
    #     cell_type_transcripts = pickle.load(f)
    # print(cell_type_transcripts[0])