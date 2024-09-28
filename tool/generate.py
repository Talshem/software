import sys
import pickle
import numpy as np
import pandas as pd
from fuzzysearch import find_near_matches
from google.cloud import storage
from joblib import Parallel, delayed
import time
from tool.switch_generator import SwitchGenerator
from tool.prokaryotic_switch_generator import ProkaryoticSwitchGenerator
from tool.window_folding_based_selection import get_gene_top_ranked_windows, get_mRNA_opening_mfe
from utils.send_mail import send_email_with_attachment as send
from utils.seqeunce_consts import GFP_GENE
from tool.window_folding_based_selection import WINDOW_SIZE
from RNA import RNA
from nupack import *
config.threads = 8


EDIT_DIST = 4
TRIGGERS_BATCH = 2
SWITCH_BATCH = 1


E_COLI_DATA_PATH = "Escherichia_coli_ASM886v2.pkl"
HUMAN_DATA_PATH = "Homo_sapiens_GRCh38.pkl"
YEAST_DATA_PATH = "Saccharomyces_cerevisiae_S288C.pkl"

DATA_PATHS = {
    "Prokaryote": E_COLI_DATA_PATH,
    "Homo sapiens": HUMAN_DATA_PATH,
    "Eukaryote": YEAST_DATA_PATH
}

"""
# for development
E_COLI_DATA_PATH = "/Users/netanelerlich/PycharmProjects/software/data/Escherichia_coli_ASM886v2.pkl"
HUMAN_DATA_PATH = "/Users/netanelerlich/PycharmProjects/software/data/Homo_sapiens_GRCh38.pkl"
YEAST_DATA_PATH = "/Users/netanelerlich/PycharmProjects/software/data/Saccharomyces_cerevisiae_S288C.pkl"

DATA_PATHS = {
    "Prokaryote": E_COLI_DATA_PATH,
    "Homo sapiens": HUMAN_DATA_PATH,
    "Eukaryote": YEAST_DATA_PATH
    }

"""

client = storage.Client()
bucket_name = 'protech_bucket'
bucket = client.get_bucket(bucket_name)
blobs = bucket.list_blobs()

def extract_top_homology_sequences(triggers_homology_mapping):
    homo_dfs = []
    for trigger_homology in triggers_homology_mapping:
        if trigger_homology:
            trigger_all_matches = []
            for gene_homology in trigger_homology:
                trigger_all_matches.extend(gene_homology.values())
            matches_df = pd.DataFrame(trigger_all_matches[0])

            mfe_dict = {'homologous_trigger_mfe': []}
            for homos in trigger_all_matches[0]:
                homo_trigger = homos['sequence']
                structure, mfe = RNA.fold(homo_trigger)
                mfe_dict['homologous_trigger_mfe'].append(mfe)
            mfe_df = pd.DataFrame(mfe_dict)

            trig_res_df = pd.concat([matches_df, mfe_df], axis=1)

        else:
            trig_res_df = pd.DataFrame(
                columns=['distance', 'idx',
                         'sequence',
                         'gene', 'protein',
                         'homologous_trigger_mfe'])

        homo_dfs.append(trig_res_df)

    # TODO: add expression levels
    # TODO: add opening mfe of trigger homology
    # RRF calc the best competitor:
    # higher is better - homology trigger mfe = open competing window, expression=higher more competing
    higher = ['homologous_trigger_mfe']
    # lower is better - distance= lower = more similar
    lower = ['distance']
    rrf_rank = [RRF(trigger_homology_df, higher, lower, index='sequence') for trigger_homology_df in homo_dfs]
    return rrf_rank


def route_input(email, target_seq, trigger, reporter_gene, cell_type, user_trigger_boo, transcripts_dict):
    results = pd.DataFrame()
    blob_name = DATA_PATHS[cell_type]

    blob = bucket.blob(blob_name)
    cell_type_transcripts = pickle.load(blob.download_as_bytes())

    # for development
    """""
    # Load transcripts for the specific cell type
    with open(blob_name, 'rb') as f:
        cell_type_transcripts = pickle.load(f)
    """""
    # Update transcripts dict
    if transcripts_dict != 'EMPTY':
        transcripts_dict = cell_type_transcripts | transcripts_dict
    else:
        transcripts_dict = cell_type_transcripts

    # Get optional triggers if the user hasn't provided one
    triggers_with_mfe = np.array([[trigger, 1]])
    if user_trigger_boo == 'EMPTY':
        cell_window = 30 if cell_type == 'Prokaryote' else 23
        optional_triggers_df = get_gene_top_ranked_windows(target_seq, window_size=cell_window)
        n_triggers = min(TRIGGERS_BATCH, len(optional_triggers_df.index))
        triggers_with_mfe = np.array(optional_triggers_df.iloc[:n_triggers, :][["window_sequence", "mfe_score"]])

    results[['trigger_window', 'mfe_score']] = triggers_with_mfe
    triggers_seqs = triggers_with_mfe[:, 0]

    # Homology search with parallelization
    s_h = time.time()
    homo_res = Parallel(n_jobs=-1)(delayed(find_homology)(trigger, transcripts_dict) for trigger in triggers_seqs)
    e_h = time.time()
    print(f'Homology search time= {e_h - s_h} seconds, n_triggers = {len(triggers_seqs)}, cell= {cell_type}')

    # Extract top homology sequences with parallelization
    rrf_ranks = extract_top_homology_sequences(homo_res)

    homology_sequences = [ranked_df['sequence'].get(0) for ranked_df in rrf_ranks]
    # Generate switches with parallelization
    s_switch = time.time()
    switch_res = Parallel(n_jobs=4)(delayed(generate_switch)(trigger, top_homology_sequence, reporter_gene, cell_type)
        for trigger, top_homology_sequence in zip(triggers_seqs, homology_sequences)
    )
    e_switch = time.time()
    print(f'Switch generation time= {e_switch - s_switch} seconds, for {SWITCH_BATCH} switches')
    results[['switch',  "complex_concentration"]] = pd.DataFrame(
        switch_res, columns=['switch', "complex_concentration"])

    # Prepare and send the report
    prepare_and_send_report(results, rrf_ranks, email)
    return


def RRF(ranking_df, higher_is_better_cols, lower_is_better_cols, index, k=60):
    isinstance(higher_is_better_cols, list)
    isinstance(lower_is_better_cols, list)
    isinstance(ranking_df, pd.DataFrame)
    isinstance(index, str)
    ranking_df = ranking_df.copy().reset_index(drop=True)

    for col in higher_is_better_cols:
        ranking_df[col + '_rank'] = ranking_df[col].rank(ascending=False)
    for col in lower_is_better_cols:
        ranking_df[col + '_rank'] = ranking_df[col].rank(ascending=True)

    ranked_columns = higher_is_better_cols + lower_is_better_cols
    ranked_columns = [f'{col}_rank' for col in ranked_columns]
    ranking_df[f'{index}_RRF'] = ranking_df[ranked_columns].apply(lambda row: sum(1 / (k + rank) for rank in row), axis=1)

    return ranking_df.sort_values(by=f'{index}_RRF', ascending=True)

def build_homology_map(trigger, seq, gene_name, protein_name):
    sequence_match_mapping = []
    matches = find_near_matches(trigger, seq, max_insertions=0, max_deletions=0, max_l_dist=EDIT_DIST)
    for match_obj in matches:
        locus_start = match_obj.start
        locus_end = match_obj.end
        sub_seq = seq[locus_start:locus_end]
        seq_location = (locus_start, locus_end)
        distance = match_obj.dist
        sequence_match_mapping.append({'distance': distance, 'idx': seq_location, 'sequence': sub_seq, 'gene': gene_name, 'protein': protein_name})
    return sequence_match_mapping


def find_homology(sequence, genome_data):
    genes_sub_sequences = []
    for gene_data_dict in genome_data:
        gene_sub_seqs = {}
        mRNA_seq = gene_data_dict['sequence']
        gene_name = gene_data_dict['gene']
        protein = gene_data_dict['protein']
        seq_match_mapping = build_homology_map(sequence, mRNA_seq, gene_name, protein)
        if seq_match_mapping:
            gene_sub_seqs[gene_name] = seq_match_mapping
            genes_sub_sequences.append(gene_sub_seqs)
    return genes_sub_sequences

def generate_switch(trigger, homologous_sequence, reporter_gene, cell_type):
    # blob = bucket.blob(cell_type)
    # content = blob.download_as_text()
    if cell_type == "Prokaryote":
        switch_generator = ProkaryoticSwitchGenerator(translated_gene_sequence=reporter_gene)
    else:
        switch_generator = SwitchGenerator(translated_gene_sequence=reporter_gene)
    switch_designed_strand, ensemble_defect, complex_concentration, switch_mfe_structure = switch_generator.get_switch(
        trigger_sequence=trigger, healthy_sequence=homologous_sequence)

    return switch_designed_strand, complex_concentration

def prepare_and_send_report(df_results, rrf_ranks, email):
    rename_dict = {
        "trigger_window": "Trigger",
        "switch": "Switch",
        "sequence": "Homologous Sequence",
        "mfe_score": "Trigger mfe Score",
        # add here regressor results
        "sequence_RRF": "Homology Score",
        "distance": "Homology Substitution Distance",
        "idx": "Homology Location",
        "gene": "Homology Gene",
        "protein": "Homology Protein"
    }

    cols_rearrange = rename_dict.values()
    combined_rows = []
    rrf_ranks = rrf_ranks.copy()
    for index, row in df_results.iterrows():
        combined_df = pd.concat([pd.DataFrame([row]), rrf_ranks[index]], ignore_index=True)
        combined_rows.append(combined_df)
    final_df = pd.concat(combined_rows, ignore_index=True).rename(columns= rename_dict)
    final_df = final_df[cols_rearrange]
    send(final_df, email)
    return


if __name__ == '__main__':
    # Get the arguments from the user form.

    s_mail = sys.argv[1]
    s_target_seq = sys.argv[2]
    s_trigger = sys.argv[3]
    s_reporter_gene = sys.argv[4]
    s_cell_type = sys.argv[5]
    s_user_trigger_boo = sys.argv[6]
    s_transcripts_dict = sys.argv[7]
    route_input(s_mail, s_target_seq, s_trigger, s_reporter_gene, s_cell_type, s_user_trigger_boo, s_transcripts_dict)
    """
    testing = {'gene': 'OR4F5',
     'protein': 'olfactory receptor 4F5',
     'sequence': 'ATGAAGAAGGTAACTGCAGAGGCTATTTCCTGGAATGAATCAACGAGTGAAACGAATAACTCTATGGTGACTGAATTCATTTTTCTGGGTCTCTCTGATTCTCAGGAACTCCAGACCTTCCTATTTATGTTGTTTTTTGTATTCTATGGAGGAATCGTGTTTGGAAACCTTCTTATTGTCATAACAGTGGTATCTGACTCCCACCTTCACTCTCCCATGTACTTCCTGCTAGCCAACCTCTCACTCATTGATCTGTCTCTGTCTTCAGTCACAGCCCCCAAGATGATTACTGACTTTTTCAGCCAGCGCAAAGTCATCTCTTTCAAGGGCTGCCTTGTTCAGATATTTCTCCTTCACTTCTTTGGTGGGAGTGAGATGGTGATCCTCATAGCCATGGGCTTTGACAGATATATAGCAATATGCAAGCCCCTACACTACACTACAATTATGTGTGGCAACGCATGTGTCGGCATTATGGCTGTCACATGGGGAATTGGCTTTCTCCATTCGGTGAGCCAGTTGGCGTTTGCCGTGCACTTACTCTTCTGTGGTCCCAATGAGGTCGATAGTTTTTATTGTGACCTTCCTAGGGTAATCAAACTTGCCTGTACAGATACCTACAGGCTAGATATTATGGTCATTGCTAACAGTGGTGTGCTCACTGTGTGTTCTTTTGTTCTTCTAATCATCTCATACACTATCATCCTAATGACCATCCAGCATCGCCCTTTAGATAAGTCGTCCAAAGCTCTGTCCACTTTGACTGCTCACATTACAGTAGTTCTTTTGTTCTTTGGACCATGTGTCTTTATTTATGCCTGGCCATTCCCCATCAAGTCATTAGATAAATTCCTTGCTGTATTTTATTCTGTGATCACCCCTCTCTTGAACCCAATTATATACACACTGAGGAACAAAGACATGAAGACGGCAATAAGACAGCTGAGAAAATGGGATGCACATTCTAGTGTAAAGTTTTAG'}
    route_input('erlichnet57@gmail.com','ATGAAGAAGGTAACTGCAGAGGCTATTTCCTGGAATGAATCAACGAGTGAAACGAATAACTCTATGGTGACTGAATTCATTTTTCTGGGTCTCTCTGATTCTCAGGAACTC','EMPTY','TCAGGAACTCCAGACCTTCCTATTTATGTTGTTTTTTGTATTCTATGGAGGAATCGT','Homo sapiens', 'EMPTY', 'EMPTY')
    """




