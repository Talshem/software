# import itertools
# import pickle
# from datetime import datetime
# from multiprocessing import Pool
# from random import random


# import matplotlib.pyplot as plt
# import pandas as pd
# from Bio.SeqUtils import gc_fraction as GC
# import sys
# from Bio import pairwise2, SeqIO
# from Bio.pairwise2 import format_alignment
# import numpy as np
# import json
# import Levenshtein as lev
# from fuzzysearch import find_near_matches
# import re
# from Bio import Align
# from Bio.pairwise2 import format_alignment
# import os
# from tool.suffix_tree import SuffixTree


# def find_similar_sequences(trigger: str, sequences_dict: dict[str:str]) -> Align.Alignment:
#     """"function to find similar sequences to a trigger sequence
#     :param trigger: trigger sequence
#     :param sequences_dict: dictionary of sequences
#     :return: alignment of the trigger sequence and the most similar sequence
#     """""

#     aligner = Align.PairwiseAligner()
#     aligner.mode = 'local'
#     aligner.match_score = 2
#     aligner.mismatch_score = -1

#     alignments = aligner.align(trigger, sequences_dict)
#     return alignments[0]

# def chunks(data_iter, size):
#     """"function to split data into chunks for processing
#     :param data_iter: data to be split
#     :param size: size of the chunk
#     :return: chunked data
#     """""
#     it = iter(data_iter)
#     for i in range(0, len(data_iter), size):
#         yield {k: data_iter[k] for k in itertools.islice(it, size)}


# def run_chunks(data_path: str, func, args):
#     """""function to run large data sets in chunks using multiprocessing
#     :param data_path: path to the data file
#     :param func: function to be applied to the data
#     :param args: arguments to the function
#     :return: None
#     """""
#     file_path = data_path
#     with open(file_path, 'rb') as file:
#         file_data = pickle.load(file)

#     # Create enumeration for the data for reproducibility and file data indicator
#     current_time = datetime.now().strftime("%H:%M:%S")
#     enumerated_data = {key: i for i, key in enumerate(file_data.keys())}
#     with open(f'enumerated_data_{current_time}.pkl', 'wb') as f:
#         pickle.dump(enumerated_data, f)

#     # Process the data in chunks
#     n_process = 5
#     for i, chunk in enumerate(chunks(file_data, n_process)):
#         with Pool(n_process) as pool:
#             args = []
#             res = pool.map(func, args)

#         # Combine the results and save of to intermediate file
#         current_time = datetime.now().strftime("%H:%M:%S")
#         if res:
#             with open(f'intermediate_results_[{i * n_process}, {(i + 1) * n_process})_{current_time}.pkl',
#                       'wb') as f:
#                 pickle.dump(res, f)


# def find_matches_fuzzy(trigg: str, data_dict, edit_distance) -> dict[str:str]:
#     """"function to find matches using fuzzy search
#     :param edit_distance: maximum edit distance allowed
#     :param trigg: trigger sequence
#     :param data_dict: dictionary of sequences key=gene name, value= rna sequence
#     :return: dictionary of sequences that match the trigger sequence
#     """""
#     sub_seqs_dict = {}
#     for seq_name, seq in data_dict.items():
#         try:
#             matches = find_near_matches(trigg, seq, max_l_dist=edit_distance)
#             for match_obj in matches:
#                 locus_start = match_obj.start
#                 end_locus = match_obj.end
#                 sub_seq = seq[locus_start:end_locus]
#                 seq_match_mapping[sub_seq] = (locus_start, end_locus)
#         except Exception as e:
#             raise f"Error in finding matches for {seq_name} with error {e}"
#             sub_seqs_dict[seq_name] = seq_match_mapping
#     return sub_seqs_dict



# def find_with_suffix_tree(trigg: str, data_dict: dict[str:str]) -> dict[str:str]:
#     """"function to find matches using suffix tree
#     :param trigg: trigger sequence
#     :param data_dict: dictionary of sequences key=gene name/rna name, value= rna/dna sequence
#     :return: dictionary of sequences that match the trigger sequence
#     """""
#     sub_seqs_dict = {}
#     for seq_name, seq in data_dict.items():
#         tree = SuffixTree(seq, seq_name)
#         try:
#             matches = tree.search(trigg)
#             for match in matches:
#                 sub_seq = (seq[match:match + len(trigg)], match)
#                 sub_seqs_dict[seq_name] = sub_seq
#         except Exception as e:
#             raise f"Error in finding matches for {seq_name} with error {e}"
#     return sub_seqs_dict




# def main(gene, trigger, rna, cell_type, file_dict):
#     """"main function to combine Peleg's code with user request
#     output to be sent by email 
#     """""
#     print(f"Gene: {gene}, Trigger: {trigger}, RNA: {rna}, Cell Type: {cell_type}")
#     print(f"File: {file_dict}")
#     return


# #peleg's code


# if __name__ == '__main__':
#     # Get the arguments from the user form and send to main function

#     import sys
#     gene = sys.argv[1]
#     trigger = sys.argv[2]
#     rna = sys.argv[3]
#     cell_type = sys.argv[4]
#     email = sys.argv[5]
#     file_dict = sys.argv[6]



#     main(gene, trigger, rna, cell_type, file_dict)



#     """
#     import time
#     # randomize sequences for testing
#     seqs_random = [''.join(np.random.choice(['A', 'C', 'G', 'T']) for _ in range(i)) for i in range(100, 100000, 1000)]

#     times_random_fuzzy_1 = []
#     times_random_suffix_1 = []
    
#     times_random_fuzzy_2 = []
#     times_random_suffix_2 = []
    
#     # for each edit distance, find the time complexity of the fuzzy search and suffix tree search
#     # randomize the trigger sequence
#     trigger = random.randrange(20, len(seq))
    
#     for i in range(1,2):
#         for seq in seqs_random:        
#             start1 = time.time()
#             find_near_matches(trigger, seq, max_l_dist=i)
#             end1 = time.time()
#             times_random_fuzzy.append(end1 - start1)
            
            
#             //change the trigger sequence to a suffix tree
#             //suffix_tree search time complexity
            
            
#         plt.plot(range(100, 100000, 1000), times_random_fuzzy)

#         plt.title(f'Fuzzy search time complexity, edit distance= {i}')
#         plt.ylabel('Time(s)')
#         plt.plot()
#         plt.show()
#         """
