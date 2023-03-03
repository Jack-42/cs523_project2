import pandas as pd
import itertools
import numpy as np

"""
@author Jack Ringer
Date: 3/2/2023
Description:
    File containing various calculations used for the report.
"""


def main():
    # calculate number of synonymous/non-synonymous changes
    df = pd.read_csv("../data/codons.txt", sep='\t')
    df["Codon"] = df["Codon"].map(lambda x: x.replace("T", "U"))
    df = df.rename({'    AminoAcid': "AminoAcid"}, axis="columns")
    mp = {}
    nucleotides = ["A", "C", "G", "U"]
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
    mp_df.to_csv("data/non_synom_freqs.csv")


if __name__ == "__main__":
    main()
