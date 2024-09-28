
import pickle
import numpy as np
import pandas as pd
from fuzzysearch import find_near_matches
from joblib import Parallel, delayed
from server import bucket
from tool.switch_generator import SwitchGenerator
from tool.prokaryotic_switch_generator import ProkaryoticSwitchGenerator
from tool.window_folding_based_selection import get_gene_top_ranked_windows, get_mRNA_opening_mfe
from utils.send_mail import send_email_with_attachment as send
from utils.seqeunce_consts import GFP_GENE


EDIT_DIST = 5
TRIGGERS_BATCH = 2
SWITCH_BATCH = 3

E_COLI_DATA_PATH = "Escherichia_coli_ASM886v2.pkl"
HUMAN_DATA_PATH = "Homo_sapiens_GRCh38.pkl"
YEAST_DATA_PATH = "Saccharomyces_cerevisiae_S288C.pkl"

DATA_PATHS = {
    "Prokaryote": E_COLI_DATA_PATH,
    "Homo sapiens": HUMAN_DATA_PATH,
    "Eukaryote": YEAST_DATA_PATH
}
def extract_top_homology_sequences(triggers_homology_mapping):
    homo_dfs = []
    for trigger_homology in triggers_homology_mapping:
        trigger_all_matches = []
        for gene_homology in trigger_homology:
            trigger_all_matches.extend(gene_homology.values())
        # TODO: optional- take only top closest distance

        if trigger_all_matches:
            top_dist_df = pd.DataFrame(
                sorted(trigger_all_matches[0], key=lambda single_homo_map: single_homo_map['distance']))
        else:
            top_dist_df = pd.DataFrame(columns=['distance', 'idx', 'sequence', 'gene'])
        homo_dfs.append(top_dist_df)

    # TODO: add expression levels
    # TODO: add opening mfe of trigger homology
    # RRF calc the best competitor:
    # higher is better - homology trigger mfe = open competing window, expression=higher more competing
    higher = []
    # lower is better - distance= lower = more similar
    lower = ['distance']
    rrf_rank = [RRF(trigger_homology_df, higher, lower, index='sequence') for trigger_homology_df in homo_dfs]
    return rrf_rank

def route_input(email, target_seq, trigger, reporter_gene, cell_type, user_trigger_boo, transcripts_dict):
    results = pd.DataFrame()
    blob_name = DATA_PATHS[cell_type]

    blob = bucket.blob(blob_name)
    cell_type_transcripts = pickle.load(blob.download_as_bytes())

    if transcripts_dict != 'EMPTY':
        transcripts_dict = cell_type_transcripts | transcripts_dict
    else:
        transcripts_dict = cell_type_transcripts

    # TODO: calc time for BATCH
    triggers_with_mfe = np.array([[trigger, 1]])
    if user_trigger_boo == 'EMPTY':
        optional_triggers_df = get_gene_top_ranked_windows(target_seq)
        n_triggers = min(TRIGGERS_BATCH, len(optional_triggers_df.index))
        triggers_with_mfe = np.array(optional_triggers_df.loc[:n_triggers, ["window_sequence", "mfe_score"]])
    results[['trigger_window', 'mfe_score']] = triggers_with_mfe

    triggers_seqs = triggers_with_mfe[:, 0]

    homo_res = Parallel(n_jobs=-1)(delayed(find_homology)(trigger, transcripts_dict) for trigger in triggers_seqs)
    rrf_ranks = extract_top_homology_sequences(homo_res)
    homology_sequences = [ranked_df['sequence'][0] for ranked_df in rrf_ranks]

    # TODO: organize input parameters after Peleg's changes
    switch_res = Parallel(n_jobs=SWITCH_BATCH)(
        delayed(generate_switch)(trigger, top_homology_sequence, reporter_gene, cell_type) for trigger, top_homology_sequence in zip(triggers_seqs[:SWITCH_BATCH], homology_sequences[:SWITCH_BATCH]))
    results[['switch', "ensemble_defect", "complex_concentration"]] = (
        pd.DataFrame(switch_res, columns=['switch', "ensemble_defect", "complex_concentration"]))
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

def build_homology_map(trigger, seq, gene_name):
    sequence_match_mapping = []
    matches = find_near_matches(trigger, seq, max_insertions=0, max_deletions=0, max_l_dist=EDIT_DIST)
    for match_obj in matches:
        locus_start = match_obj.start
        locus_end = match_obj.end
        sub_seq = seq[locus_start:locus_end]
        seq_location = (locus_start, locus_end)
        distance = match_obj.dist
        sequence_match_mapping.append({'distance': distance, 'idx': seq_location, 'sequence': sub_seq, 'gene': gene_name})
    return sequence_match_mapping

def find_homology(sequence, genome_data):
    genes_sub_sequences = []
    for gene_data_dict in genome_data:
        gene_sub_seqs = {}

        gene_name = gene_data_dict['gene']
        mRNA_seq = gene_data_dict['sequence']
        seq_match_mapping = build_homology_map(sequence, mRNA_seq, gene_name)
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

    return switch_designed_strand, ensemble_defect,complex_concentration

def prepare_and_send_report(df_results, rrf_ranks, email):
    #print(df_results)
    #print(rrf_ranks)
    #print(email)
    combined_rows = []
    rrf_ranks = rrf_ranks.copy()
    for index, row in df_results.iterrows():
        combined_df = pd.concat([pd.DataFrame([row]), rrf_ranks[index]], ignore_index=True)
        combined_rows.append(combined_df)
    final_df = pd.concat(combined_rows, ignore_index=True)
    print(final_df)
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
    print('sent to route_input')
    route_input(s_mail, s_target_seq, s_trigger, s_reporter_gene, s_cell_type, s_user_trigger_boo, s_transcripts_dict)
    
