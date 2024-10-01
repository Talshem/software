import numpy as np
import re
import pickle


def build_suffix_array(elements):
    protein_pattern = r'\[protein=(.*?)\]'
    gene_pattern = r'\[gene=(.*?)\]'

    suffixes_dict = []
    for item in elements:
        try:
            sequence = item.split('[gbkey=CDS]')[1].replace('\n', '')
            protein = re.findall(protein_pattern, item)[0]
            gene = re.findall(gene_pattern, item)[0]
            suffixes = [(sequence[i:], i) for i in range(len(sequence))]
            suffixes.sort(key=lambda x: x[0])
            suffixes_dict.append({'gene': gene, 'protein': protein, 'array': [suffix[1] for suffix in suffixes], 'sequence': sequence})
        except:
            continue


    return suffixes_dict



if __name__ == '__main__':


    target = ''.join(np.random.choice(['A', 'C', 'G', 'T'], 4))
    with open('/webTool/data/Saccharomyces_cerevisiae_S288C.pkl', 'rb') as file:
        data = pickle.load(file)
    """suffixes = build_suffix_array(transcriptome)
       with open('Escherichia_coli_ASM886v2.pkl', 'wb') as f:
       pickle.dump(suffixes, f)
       file_path = "./webTool/search/cds_from_genomic.fna"
       with open(file_path, "r") as file:
       annotations = file.read()
       transcriptome = annotations.split('>lcl')[1:]
       sequence_pattern = re.compile(r'([ACGT]+)')"""





