import pandas as pd

"""
@author Jack Ringer
Date: 3/15/2023
Description:
    Various utility functions
"""


def load_codon_whl(csv_path: str):
    df = pd.read_csv(csv_path, sep='\t')
    df = df.rename({'    AminoAcid': "AminoAcid"}, axis="columns")
    return df


def load_fasta_seq(fasta_path: str):
    with open(fasta_path) as f:
        lines = f.readlines()
        lines = lines[1:]  # exclude header
    lines = list(map(lambda x: x.strip(), lines))  # remove \n chars
    genome = ''.join(lines)
    return genome


def acid_seq_to_genome(seq: str):
    """
    Convert a sequence of amino acids to a sequence of nucleotides.
    :param seq: str, sequence of amino acids
    :return: str, converted genome
    """
    codon_df = load_codon_whl("../data/codons.txt")
    map_fn = lambda acid: codon_df[codon_df["Letter"] == acid]["Codon"].iloc[0]
    nt_genome = list(map(map_fn, list(seq)))
    nt_genome = ''.join(nt_genome)
    return nt_genome
