import pickle
import random
import time
import numpy as np
import pickle

from tool.generate import find_homology
from itertools import combinations
import matplotlib.pyplot as plt
from tool.generate import find_near_matches, build_homology_map


def homo_sapiens_search_plot(edit_distance, trigger_len=23):
    path = 'data/Homo_sapiens_GRCh38.pkl'
    with open(path, 'rb') as f:
        data = pickle.load(f)
    start = time.time()

    for i in range(10):
        trigger = ''.join([random.choice(['A', 'C', 'G', 'T']) for _ in range(trigger_len)])
        matches_over_genes = {}
        for dict in data:
            gene = dict['gene']
            transcript = dict['sequence']
            res = build_homology_map(trigger, transcript, gene, edit_distance)
            if res:
                matches_over_genes[gene] = res
    end = time.time()
    tot = end - start
    print(f" T={tot}")

    """
    plt.title('Transcriptome search time for different trigger lengths using fuzzy search')
    plt.xlabel('Length of trigger')
    plt.ylabel('Time in seconds')
    plt.show()
    """

def gen_mutations(trigger: str, sub: int):
    trig_len = len(trigger)
    indexes_permutation = list(combinations(range(trig_len), sub))
    nucleotides = ['A', 'C', 'G', 'T']
    mutated_triggers = {}
    for mutated_index in indexes_permutation:
        trig_char_list = list(trigger)
        for nuc_idx in mutated_index:
            choices = nucleotides.copy()
            nuc = trig_char_list[nuc_idx]
            choices.remove(nuc)
            trig_char_list[nuc_idx] = random.choice(choices)
        mutated_trigger = ''.join(trig_char_list)
        mutated_triggers[mutated_trigger] = mutated_index
    return mutated_triggers

def construct_dummy_seq(trigger_muts):
    trig_n = len(trigger_muts)
    n_list = random.sample(range(1, trig_n * 10), trig_n)
    fillers = []

    for x in n_list:
        filler_seq = ''.join(random.choice(['A', 'C', 'G', 'T']) for _ in range(x))
        fillers.append(filler_seq)

    dummy_seq = ''
    mut_trigger_mapping = {}
    for trig, mut_index in trigger_muts.items():
        dummy_seq += fillers.pop(0)
        start_index = len(dummy_seq)
        dummy_seq += trig
        end_index = len(dummy_seq)
        mut_trigger_mapping[trig] = (start_index, end_index)

    return dummy_seq, mut_trigger_mapping

def test_fuzzy_search(trig_len: int, sub: int):
    trigger = ''.join([random.choice(['A', 'C', 'G', 'T']) for _ in range(trig_len)])
    mutated_triggers = gen_mutations(trigger, sub)
    test_seq, trig_muts = construct_dummy_seq(mutated_triggers)
    test = find_near_matches(trigger, test_seq, max_substitutions=sub, max_deletions=0, max_insertions=0, max_l_dist=sub)
    test_res = {}
    for match in test:
        start = match.start
        end = match.end
        sub_seq = test_seq[start:end]
        test_res[sub_seq] = (start, end)

    return test_res, trig_muts, test_seq


if __name__ == "__main__":
    # Search test validation

    """
    random.seed(42)
    substitutions = 5
    triggers_len = 25
    res_array = []
    for trig_len in range(15, triggers_len, 2):
        for substitutions in range(5):
            test_res, mutations_mapping, dummy_seq = test_fuzzy_search(trig_len=trig_len, sub=substitutions)
            for mutated_trigger, mut_loc in mutations_mapping.items():
                for search_res_trigger, res_loc in test_res.items():
                    if mut_loc[0] == res_loc[0] and mut_loc[1] == res_loc[1]:
                        mutations_mapping[mutated_trigger] = all(np.array(list(mutated_trigger)) == np.array(list(search_res_trigger)))
            print(f'for trigger len={trig_len} with {substitutions} subs, all possible mutated triggers found: {all(mutations_mapping.values())}')
            res_array.append(f'for trigger len={trig_len} with {substitutions} subs, all possible mutated triggers found: {all(mutations_mapping.values())}')
    """

