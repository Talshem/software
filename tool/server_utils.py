import os
import re
from io import StringIO

from Bio import SeqIO
from flask import abort, current_app
from werkzeug.utils import secure_filename
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO



def process_file_stream(file_stream):
    genes_list = []
    # Decode the file stream to text if needed (depending on its format)
    file_content = file_stream.read().decode("utf-8")
    handle = StringIO(file_content)
    for record in SeqIO.parse(handle, "fasta"):
        seq_record = {}
        seq_record['gene'] = record.id
        try:
            seq_record['sequence'] = validate_sequence(str(record.seq))
        except ValueError as e:
            raise ValueError(f"Error in processing file: {e}")

        genes_list.append(seq_record)
    return genes_list

def validate_sequence(sequence):
    if not re.search(r"^[ACGT]", sequence):
        raise ValueError("Invalid sequence: must contain only A, C, G, or T nucleotides.")

    return sequence






def dict_to_fasta(genes_list, output_file):
    # Create a list of SeqRecord objects
    records = []

    for d in genes_list:
        gene_name = d['gene']
        sequence = d['sequence']
        seq_obj = Seq(sequence)
        record = SeqRecord(seq_obj, id=gene_name, description="")
        records.append(record)

    with open(output_file, "w") as fasta_out:
        SeqIO.write(records, fasta_out, "fasta")





