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


if __name__ == "__main__":
    seq1 = load_fasta_seq("../data/covid_spike_protein.fasta")
    g = acid_seq_to_genome(seq1, "../data/cov_spike_nt_genome.fasta", "Covid-19 Spike Protein")
