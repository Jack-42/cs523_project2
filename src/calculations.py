import pandas as pd
import itertools
import numpy as np

"""
@author Jack Ringer
Date: 3/2/2023
Description:
    File containing various calculations used for the report.
"""


def calc_non_synom_freqs():
    """
    Calculate frequency of synonymous/non-synonymous changes and
    save result to a csv file.
    :return: Non
    """
    df = pd.read_csv("../data/codons.txt", sep='\t')
    df = df.rename({'    AminoAcid': "AminoAcid"}, axis="columns")
    mp = {}
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
        percent_non_synom = non_synom / (synom + non_synom)
        key = "%s:%s" % (base, mutant)
        mp[key] = percent_non_synom
    mp_df = pd.DataFrame(mp.items(), columns=["Base:Mut", "Freq"])
    mp_df["Freq"] = mp_df["Freq"].map(lambda x: np.round(x, decimals=4))
    mp_df.to_csv("../data/non_synom_freqs.csv", index=False)


def batched_non_synom_freqs():
    """
    Calculate frequency of synonymous/non-synonymous changes and
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
    with open(fasta_pth) as f:
        lines = f.readlines()
        lines = lines[1:]  # exclude header
    lines = list(map(lambda x: x.strip(), lines))  # remove \n chars
    genome = ''.join(lines)
    genome_len = len(genome)
    freqs = {}
    total = 0
    for nt in ['A', 'C', 'T', 'G']:
        freqs[nt] = genome.count(nt) / genome_len
        total += freqs[nt]
    assert total == 1.0
    return freqs


if __name__ == "__main__":
    batched_non_synom_freqs()
    """
    calc_non_synom_freqs()
    freqs = cov19_nt_freqs("../data/covid_seq.fasta")

    rounded_freqs = {k: np.round(v, decimals=4) for k, v in freqs.items()}
    print(rounded_freqs)
    """
