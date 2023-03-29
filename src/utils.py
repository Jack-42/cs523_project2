import pandas as pd
import numpy as np

from binding_calculator import BindingCalculator

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


def gen_fasta_file(fpath: str, seqs: list, labels: list):
    """
    Generate a FASTA file with the given params.
    :param str fpath: path to save fasta file to
    :param list seqs: sequences to write
    :param labels: corresponding labels for each sequence
    """
    with open(fpath, "w") as ofile:
        for seq, label in zip(seqs, labels):
            ofile.write(">" + label + "\n" + seq + "\n")
    print("Successfully saved given seqs/ids to:", fpath)


def acid_seq_to_genome(seq: str, save_path: str = None, label: str = ""):
    """
    Convert a sequence of amino acids to a sequence of nucleotides.
    :param seq: str, sequence of amino acids
    :param save_path: str (optional), path to save full genome to
    :param label: str (optional), label for fasta file
    :return: str, converted genome
    """
    codon_df = load_codon_whl("../data/codons.txt")
    map_fn = lambda acid: codon_df[codon_df["Letter"] == acid]["Codon"].iloc[0]
    nt_genome = list(map(map_fn, list(seq)))
    nt_genome = ''.join(nt_genome)
    if save_path is not None:
        gen_fasta_file(save_path, [nt_genome], [label])
    return nt_genome


def generate_bloom_sites(save_path: str):
    """
    Generate list of sites where data is available for the bloom calculator for
    any given virus.
    :param save_path: str, path to save list to
    :return: list, list of sites
    """
    escape_data_csv = "../data/escape_calculator_data.csv"
    # antibodies = ['BA.2.75', 'BQ.1.1', 'BA.2', 'XBB', 'BA.1', 'BA.5', 'D614G']
    viruses = ['pre-Omicron SARS-CoV-2', 'SARS-CoV-1 then SARS-CoV-2',
               'pre-Omicron SARS-CoV-2 then BA.2',
               'pre-Omicron SARS-CoV-2 then BA.5', 'SARS-CoV-2',
               'pre-Omicron SARS-CoV-2 then BA.1']
    # use only sites with data
    valid_sites = set(np.arange(0, 10000))
    for vir in viruses:
        cal = BindingCalculator(csv_or_url=escape_data_csv,
                                eliciting_virus=vir)
        valid_sites = valid_sites & set(cal.sites)
    valid_sites = list(valid_sites)
    np.save(save_path, valid_sites)
    return valid_sites


if __name__ == "__main__":
    generate_bloom_sites("../results/bloom_valid_sites.npy")
