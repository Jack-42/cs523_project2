import pandas as pd
import numpy as np

"""
@author Jack Ringer
Date: 3/4/2023
Description:
    Script used to answer question 2 on the report
"""


def main():
    n_muts = 89709
    covid_freq_df = pd.read_csv("../data/covid_freqs.csv")
    nonsynom_freq_df = pd.read_csv("../data/batch_non_synom_freqs.csv")
    probs_sum = np.dot(covid_freq_df["Frequency"], nonsynom_freq_df["Frequency"])
    non_synom_muts = np.round(n_muts * probs_sum, decimals=0)
    synom_muts = n_muts - non_synom_muts
    print("# of non-synonymous mutations:", non_synom_muts)
    print("# of synonymous mutations:", synom_muts)

if __name__ == "__main__":
    main()
