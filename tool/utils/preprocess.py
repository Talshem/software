import itertools
import pickle
from datetime import datetime
from multiprocessing import Pool

from Bio import Align
from .utils.suffix_tree import Suffixtree


def find_similar_sequences_bio(trigger: str, sequences_dict: dict[str:str]) -> Align.Alignment:
    """"function to find similar sequences to a trigger sequence
    :param trigger: trigger sequence
    :param sequences_dict: dictionary of sequences
    :return: alignment of the trigger sequence and the most similar sequence
    """""

    aligner = Align.PairwiseAligner()
    aligner.mode = 'local'
    aligner.match_score = 2
    aligner.mismatch_score = -1

    alignments = aligner.align(trigger, sequences_dict)
    return alignments[0]

def chunks(data_iter, size):
    """"function to split data into chunks for processing
    :param data_iter: data to be split
    :param size: size of the chunk
    :return: chunked data
    """""
    it = iter(data_iter)
    for i in range(0, len(data_iter), size):
        yield {k: data_iter[k] for k in itertools.islice(it, size)}


def run_chunks(data_path: str, func, args):
    """""function to run large data sets in chunks using multiprocessing
    :param data_path: path to the data file
    :param func: function to be applied to the data
    :param args: arguments to the function
    :return: None
    """""
    file_path = data_path
    with open(file_path, 'rb') as file:
        file_data = pickle.load(file)

    # Create enumeration for the data for reproducibility and file data indicator
    current_time = datetime.now().strftime("%H:%M:%S")
    enumerated_data = {key: i for i, key in enumerate(file_data.keys())}
    with open(f'enumerated_data_{current_time}.pkl', 'wb') as f:
        pickle.dump(enumerated_data, f)

    # Process the data in chunks
    n_process = 5
    for i, chunk in enumerate(chunks(file_data, n_process)):
        with Pool(n_process) as pool:
            args = []
            res = pool.map(func, args)

        # Combine the results and save of to intermediate file
        current_time = datetime.now().strftime("%H:%M:%S")
        if res:
            with open(f'intermediate_results_[{i * n_process}, {(i + 1) * n_process})_{current_time}.pkl',
                      'wb') as f:
                pickle.dump(res, f)


def find_with_suffix_tree(trigg: str, data_dict: dict[str:str]) -> dict[str:str]:
    """"function to find matches using suffix tree
    :param trigg: trigger sequence
    :param data_dict: dictionary of sequences key=gene name/rna name, value= rna/dna sequence
    :return: dictionary of sequences that match the trigger sequence
    """""
    sub_seqs_dict = {}
    for seq_name, seq in data_dict.items():
        tree = Suffixtree(seq, seq_name)
        try:
            matches = tree.search(trigg)
            for match in matches:
                sub_seq = (seq[match:match + len(trigg)], match)
                sub_seqs_dict[seq_name] = sub_seq
        except Exception as e:
            raise f"Error in finding matches for {seq_name} with error {e}"
    return sub_seqs_dict


if __name__ == '__main__':
    pass

