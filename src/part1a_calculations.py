import pandas as pd
import itertools
import numpy as np
from utils import load_fasta_seq

"""
@author Jack Ringer
Date: 3/2/2023
Description:
    File containing various calculations used for part 1a of the assignment.
"""


def averaged_non_synom_freqs():
    """
    Calculate averaged frequencies of non-synonymous changes and
    save result to a csv file.
    :return: Non
    """
    df = pd.read_csv("../data/codons.txt", sep='\t')
    df = df.rename({'    AminoAcid': "AminoAcid"}, axis="columns")
    non_synom_cnts = {}
    synom_cnts = {}
    nucleotides = ["A", "C", "G", "T"]
    for base, mutant in itertools.product(nucleotides, nucleotides):
        if base == mutant:
            continue
        synom = 0
        non_synom = 0
        for i in range(0, 3):
            df_sub = df[
                df["Codon"].map(lambda x: True if x[i] == base else False)]
            # df_sub = df[df["Codon"].str.startswith(base)]
            df_mut = pd.DataFrame({})
            df_mut["Codon"] = df_sub["Codon"].map(
                lambda x: x[:i] + x[i].replace(base, mutant) + x[i + 1:])
            # see where amino acid changed
            for j in range(len(df_mut["Codon"])):
                cdn = df_mut["Codon"].iloc[j]
                new_acid = df[df["Codon"] == cdn]["AminoAcid"].iloc[0]
                old_acid = df_sub["AminoAcid"].iloc[j]
                if new_acid == old_acid:
                    synom += 1
                else:
                    non_synom += 1
        non_synom_cnts[base] = non_synom_cnts.get(base, 0) + non_synom
        synom_cnts[base] = synom_cnts.get(base, 0) + synom
    freq_mp = {}
    for nt in nucleotides:
        ns_cnt = non_synom_cnts[nt]
        s_cnt = synom_cnts[nt]
        freq_mp[nt] = ns_cnt / (ns_cnt + s_cnt)
    df = pd.DataFrame(freq_mp.items(), columns=["Nucleotide", "Frequency"])
    df["Frequency"] = df["Frequency"].map(lambda x: np.round(x, decimals=4))
    df.to_csv("../data/batch_non_synom_freqs.csv", index=False)


def cov19_nt_freqs(fasta_pth: str):
    """
    Get frequency of nucleotides (nt) in original Covid-19 genome.
    :param fasta_pth: path to .fasta file containing Covid-19 genome.
    :return: dict, where dict[nt] = frequency of nt in genome
    """
    genome = load_fasta_seq(fasta_pth)
    genome_len = len(genome)
    freqs = {}
    total = 0
    for nt in ['A', 'C', 'T', 'G']:
        freqs[nt] = genome.count(nt) / genome_len
        total += freqs[nt]
    assert total == 1.0
    df = pd.DataFrame(freqs.items(), columns=["Nucleotide", "Frequency"])
    df.to_csv("../results/covid_freqs.csv", index=False)
    return freqs

def cov19_rbd_nt_freqs(fasta_pth: str):
    """
    Get frequency of nucleotides (nt) in original Covid-19 RBD. (same as above function with different path.)
    :param fasta_pth: path to .fasta file containing Covid-19 RBD.
    :return: dict, where dict[nt] = frequency of nt in genome
    """

    #first 1269 (0:1268) nt's are not part of it,
    #the following 669 are the RBD

    genome = load_fasta_seq(fasta_pth)
    RBD_genome = genome[1269:1938]
    print(genome_len)
    freqs = {}
    total = 0
    for nt in ['A', 'C', 'T', 'G']:
        freqs[nt] = RBD_genome.count(nt) / genome_len
        total += freqs[nt]
    assert total == 1.0
    df = pd.DataFrame(freqs.items(), columns=["Nucleotide", "Frequency"])
    df.to_csv("../results/covid_rbd_freqs.csv", index=False)
    return freqs


def q2_main():
    """
    Script used to answer question 2 of part 1a on the report
    :return:
    """
    n_muts = 89709
    covid_freq_df = pd.read_csv("../results/covid_freqs.csv")
    nonsynom_freq_df = pd.read_csv("../results/batch_non_synom_freqs.csv")
    probs_sum = np.dot(covid_freq_df["Frequency"],
                       nonsynom_freq_df["Frequency"])
    non_synom_muts = np.round(n_muts * probs_sum, decimals=0)
    synom_muts = n_muts - non_synom_muts
    print("# of non-synonymous mutations:", non_synom_muts)
    print("# of synonymous mutations:", synom_muts)


if __name__ == "__main__":
    averaged_non_synom_freqs()
    # calc_non_synom_freqs()
    freqs = cov19_nt_freqs("../results/covid_seq.fasta")
    freqs_rbd = cov19_rbd_nt_freqs("../results/cov_spike_nt_genome.fasta")

    rounded_freqs = {k: np.round(v, decimals=4) for k, v in freqs.items()}
    rounded_rbd_freqs = {k: np.round(v, decimals=4) for k, v in freqs_rbd.items()}
    print(rounded_freqs)
    print(rounded_rbd_freqs)
