import numpy as np
import pandas as pd
from fuzzysearch import find_near_matches
import seaborn as sns
from matplotlib import pyplot as plt
import pickle
from joblib import Parallel, delayed

EDIT_DIST = 15
def create_data(path_green,path_trans):
    table = pd.read_csv(path_green)
    transcripts = pd.read_csv(path_trans, sep='\t')
    table.rename(columns={'Unnamed: 0': 'id'}, inplace=True)
    data = []
    for seq_name, transcript in zip(transcripts['Transcription-Units'], transcripts['Sequence - DNA sequence']):
        d = {'gene': seq_name, 'sequence': transcript.replace('T', 'U')}
        data.append(d)
    return table, data

def build_seq_homology_dict(trigger, seq):
    return find_near_matches(trigger, seq, max_insertions=0, max_deletions=0, max_l_dist=EDIT_DIST)

def find_homology(trig, genome_data) -> dict[str:str]:
    """"function to find matches using fuzzy search
    :param trig: trigger sequence
    :param genome_data: dictionary of sequences key=gene name, value= rna sequence
    :return: dictionary of sequences that match the trigger sequence
    """""
    sub_seqs_dict = {}
    results = Parallel(n_jobs=8)(delayed(build_seq_homology_dict)(trig, gene_data_dict['sequence']) for gene_data_dict in genome_data)
    for gene_name, seq_match_mapping in zip([gene_data_dict['gene'] for gene_data_dict in genome_data], results):
        if seq_match_mapping:
            sub_seqs_dict[gene_name] = seq_match_mapping
    return sub_seqs_dict

def create_hist_time(table, data):
    sub_seqs = []
    for idx, row in table.iterrows():
        trigger = row['trigger']
        id = row['id']
        trig_dict = {"id": id, 'homo_dict': find_homology(trigger, data)}
        sub_seqs.append(trig_dict)
    return sub_seqs

def random_choice_trans(green_data):
    """RANDOM CHOICE OF TRIGGERS n=9000"""
    # Method 1 empirical dist of off_value
    #vals = green_data.sort_values('off_value').reset_index().rename(columns={'off_value': 'prob'})['prob']
    #emp_prob = vals.value_counts(normalize=True).reset_index().rename(columns={'index': 'off_value'})
    #greenp = green_data.merge(emp_prob, on='off_value')
    #idxs = greenp.sample(n=9000, weights='prob').index
    #rand_data = green_data.loc[idxs]

    # Method 2 - built in python func
    rand_data = green_data.sample(n=90000, random_state=42)
    return rand_data







if __name__ == "__main__":
    path_green = '/Users/netanelerlich/PycharmProjects/webTool/data/triggers_with_toeholds.csv'
    path_trans = "/Users/netanelerlich/PycharmProjects/webTool/data/Copy-of-All-transcription-units-of-E.-coli-BL21(DE3) - names.tsv"
    green_data, transcripts = create_data(path_green, path_trans)
    random_triggers = random_choice_trans(green_data)

    # SEARCH IN DATA
    """
    def run_slice(i):
        slice = random_triggers[i*100:(i+1)*100]
        results = Parallel(n_jobs=8)(delayed(create_hist_time)(slice[j*20:(j+1)*20], transcripts) for j in range(5))
        with open(f'data slices sample/sub_seqs_slice_{i}_df_sample.pkl', 'wb') as f:
            pickle.dump(results, f)
        return
    Parallel(n_jobs=8, verbose=51)(delayed(run_slice)(i) for i in range(438, 900))
    """

    # TRANSFORM DATA
    """
    data = {'id': [], 'count': [], 'mean_dist': []}
    for i in range(90):
        with (open(f'data slices sample/sub_seqs_slice_{i}_df_sample.pkl', 'rb') as f):
            results = pickle.load(f)
            for process_results in results:
                for triggs_matches_dict in process_results:# process_results is a list of lists, because of the parallelization.
                    trigg_id = triggs_matches_dict['id']
                    homo_dict = triggs_matches_dict['homo_dict']

                    count = 0
                    tot_dist = 0
                    for seq, match_list in homo_dict.items():
                        for match in match_list:
                            if match.dist <= 8:
                                tot_dist += match.dist
                                count += 1
                                
                    data['id'].append(trigg_id)
                    data['count'].append(count)
                    if count >= 1:
                        data['mean_dist'].append(tot_dist/count)
                    else:
                        data['mean_dist'].append(None)

    count_frame = pd.DataFrame.from_dict(data)
    count_frame.to_csv('count_frame_sample_15_subs.csv')
    """


    # CREATE HIST
    """
    hist = [0 for _ in range(20)]
    for i in range(100):
        with open(f'/Users/netanelerlich/PycharmProjects/webTool/utils/data_slices_emp_off/sub_seqs_slice_{i}_size_10.pkl', 'rb') as f:
            results = pickle.load(f)
            for process_results in results:
                # each process result is list of triggers matches dict [{"id": trigger_id, "homo_seqs" : {seq_x: [match1,match2..], seq_y:[match1,match2..]}}]
                for triggers_matches_dict in process_results:
                    trigger_dists_hist = set()
                    trigger_id = triggers_matches_dict['id']
                    homo_dict = triggers_matches_dict['homo_dict']
                    for seq_with_matches, match_list in homo_dict.items():
                        match_dis = set([match.dist for match in match_list])
                    trigger_dists_hist = trigger_dists_hist.union(match_dis)

                for num in trigger_dists_hist:
                    hist[num] += 1

    with open(f'hist_15_subs_off_emp.pkl', 'wb') as f:
        pickle.dump(hist, f)

    hist = []
    with open('hist_15_subs_off_emp.pkl', 'rb') as d:
        hist = pickle.load(d)
 
    plt.figure(figsize=(8, 6))
    bins = list(range(len(hist)))
    plt.bar(bins, hist, edgecolor='black')
    plt.grid()
    plt.yticks(range(1, hist[15]+10, 20))
    plt.title('Homology Count Over Subs')
    plt.ylabel('Count')
    plt.xlabel('Substitutions')
    plt.show()
    """


    #COUNT CORR ANALYSIS
    """
    import seaborn as sns
    import matplotlib.pyplot as plt
    count_frame = pd.read_csv('/Users/netanelerlich/PycharmProjects/webTool/utils/data slices sample/count_frame_sample_8_subs.csv')
    green_counts = green_data.merge(count_frame, on='id')
    green_counts['count'] = green_counts['count'].apply(lambda x: 0 if x is pd.NA else x)
    green_counts['mean_dist'] = green_counts['mean_dist'].apply(lambda x: 0 if x is pd.NA else x)
    green_counts['count'] = green_counts['count'].astype('int')/sum(green_counts['count'])
    corr_data = green_counts[['off_value', 'on_value', 'count', 'mean_dist']].sort_values('count')
    corr_data = corr_data.reset_index().drop(columns=['index'])
    corr_mat = corr_data.corr()
    sns.heatmap(corr_data)
    plt.show()

    """


