import numpy as np
import pickle
import pandas as pd

def is_within_mismatches(suffix, target, max_mismatches):
    if len(suffix) != len(target):
        return False
    mismatches = 0
    for s_char, t_char in zip(suffix, target):
        if s_char != t_char:
            mismatches += 1
            if mismatches > max_mismatches:
                return False
    return True

def is_in_genome(target, target_len, suffix_array, sequence, max_mismatches):
    positions = []
    low = 0
    high = len(suffix_array) - 1
    while low < high:
        mid = (low + high) // 2
        position = int(suffix_array[mid])
        suffix = sequence[position:position+target_len]
        if is_within_mismatches(suffix, target, max_mismatches):
            positions.append((position, position+target_len))
        elif target < suffix:
            high = mid - 1
        else:
            low = mid + 1
    return positions if len(positions) > 0 else False

def get_similar_transcripts(target, suffix_dict, max_mismatches=2):
    candidates = set()
    target_len = len(target)
    for row in suffix_dict:
        gene, protein, array, sequence = row['gene'], row['protein'], row['array'], row['sequence']
        positions = is_in_genome(target, target_len, array, sequence, max_mismatches)
        if positions:
            candidates.add({'gene': gene, 'protein': protein, 'sequence': sequence, 'positions': positions})
            break
    return candidates


if __name__ == '__main__':
    pass


