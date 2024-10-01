
import re
from io import StringIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO


def process_file_stream(file_stream, file_format):
    parse_inputs_dict = {
        ".fasta": "fasta",  # FASTA format
        ".fa": "fasta",  # FASTA format (alternative extension)
        ".fastq": "fastq",  # FASTQ format
        ".gb": "genbank",  # GenBank format
        ".gbk": "genbank",  # GenBank format (alternative extension)
        ".embl": "embl",  # EMBL format
        ".phy": "phylip",  # PHYLIP format
        ".phylip": "phylip",  # PHYLIP format (alternative extension)
        ".aln": "clustal",  # ClustalW format
        ".nex": "nexus",  # NEXUS format
        ".stockholm": "stockholm",  # Stockholm format
        ".tab": "tab",  # Tab-delimited format
        ".qual": "qual",  # FASTA quality scores
        ".abi": "abi",  # ABI chromatogram files
        ".sff": "sff",  # Standard Flowgram Format
        ".xml": "uniprot-xml",  # Uniprot XML format
        ".ig": "ig"  # IntelliGenetics format
    }
    parse_format = parse_inputs_dict[file_format]
    genes_list = []

    # Decode the file stream to text if needed (depending on its format)
    file_content = file_stream.read().decode("utf-8")
    handle = StringIO(file_content)

    for record in SeqIO.parse(handle, parse_format):
        seq_record = {}
        seq_record['gene'] = record.id
        try:
            seq_record['sequence'] = validate_sequence(str(record.seq))
        except ValueError as e:
             # TODO: send mail to user noteing the error?
             raise ValueError(f"Error in processing file: {e}")
        seq_record['protein'] = ''
        genes_list.append(seq_record)
    return genes_list

def validate_sequence(sequence):
    if not re.search(r"^[ACGUTactgu]", sequence):
        raise ValueError("Invalid sequence: must contain only A, C, G, or T, U , a ,c,t,g,u nucleotides.")
    return sequence



def validate_sequence_bool(sequence):
    return bool(re.fullmatch(r"[ACGUTacgut]+", sequence))





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





