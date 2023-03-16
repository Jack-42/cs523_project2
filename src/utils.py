"""
@author Jack Ringer
Date: 3/15/2023
Description:
    Various utility functions
"""


def load_fasta_seq(fasta_path: str):
    with open(fasta_path) as f:
        lines = f.readlines()
        lines = lines[1:]  # exclude header
    lines = list(map(lambda x: x.strip(), lines))  # remove \n chars
    genome = ''.join(lines)
    return genome
