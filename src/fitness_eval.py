import pandas as pd
import numpy as np

"""
@author Jack Ringer
Date: 3/18/2023
Description:
    This file contains code used to evaluate the antigenic fitness of mutations
    1-jump away from our generated neutral network using the work of 
    J Bloom Lab:
    https://github.com/jbloomlab/SARS2-mut-fitness/blob/main/
"""


def get_fitness(fitness_df: pd.DataFrame, mut_code: str):
    aa = mut_code[-1]
    site = int(mut_code[1:-1])
    fitness_entry = fitness_df[
        (fitness_df["aa"] == aa) & (fitness_df["aa_site"] == site)]
    if len(fitness_entry) == 0:
        # fitness not recorded for this mutation
        return None
    # many mutations have multiple entries, entry with the largest expected
    # count is most accurate
    best_entry = fitness_entry.loc[fitness_entry["expected_count"].idxmax()]
    return best_entry["fitness"]


def main(save_pth: str):
    assert save_pth.endswith(".csv"), "Check save_pth"
    # from:
    # https://github.com/jbloomlab/SARS2-mut-fitness/blob/main/data/Neher_aa_fitness.csv
    fitness_df = pd.read_csv("../data/aa_fitness.csv")
    fitness_df = fitness_df[["aa_site", "aa", "fitness", "expected_count"]]
    fitness_df = fitness_df.dropna()
    edge_muts = np.load("../results/neutral_net_edgemuts.npy",
                        allow_pickle=True)
    mut_fn = lambda mut_code: get_fitness(fitness_df, mut_code)
    edge_fitnesses = list(map(mut_fn, edge_muts))
    mut_fit_mp = dict(zip(edge_muts, edge_fitnesses))
    mut_fit_df = pd.DataFrame(mut_fit_mp.items(), columns=["Mutation", "Fitness"])
    mut_fit_df = mut_fit_df.sort_values(by="Fitness")
    mut_fit_df = mut_fit_df.dropna()  # values not found by fitness calculator
    mut_fit_df.to_csv(save_pth, index=False)
    print("Saved fitness info to:", save_pth)


if __name__ == "__main__":
    main("../results/edge_fitness.csv")
