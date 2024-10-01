import itertools
from random import choice


class Suffixtree:
    def __init__(self, text, gene_name):
        self.name = gene_name
#        self.tree = STree.STree(text)

    def search(self, pattern):
        result = self.tree.find_all(pattern)
        return result



def induce_mutations(trigger, max_edit_distance):
    nucleotides = ['A', 'C', 'G', 'T']
    mutations = []
    seq_length = len(trigger)
    seq = list(trigger)
    for i in range(1, max_edit_distance + 1):
        for pos in itertools.combinations(range(seq_length), i):
            mutated_seq = seq.copy()
            # replace the nucleotide at the position with different nucleotide
            for p in pos:
                mutated_seq[p] = choice([n for n in nucleotides if n != mutated_seq[p]])
                mutations.append(''.join(mutated_seq))  # Add mutated sequence to the set
    return mutations

if __name__ == "__main__":

    seq = ''.join(choice(['A', 'C', 'G', 'T']) for _ in range(10))
    tree = Suffixtree(seq, 'Test')


