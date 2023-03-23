import numpy as np
import pandas as pd

from binding_calculator import BindingCalculator

"""
@author Jack Ringer
Date: 3/18/2023
Description:
    This file contains code used to evaluate the antigenic escape of mutations
    1-jump away from our generated neutral network using the work of 
    J Bloom Lab:
    https://github.com/jbloomlab/SARS2_RBD_Ab_escape_maps
"""


def calculate_difs(df: pd.DataFrame, viruses: list, edge_sites: list):
    rel_table = np.zeros((len(edge_sites), len(edge_sites)))
    for i in range(len(edge_sites)):
        site1 = edge_sites[i]
        for j in range(len(edge_sites)):
            site2 = edge_sites[j]
            diffs = np.abs(df.loc[site1] - df.loc[site2])  # arr of diffs
            rel_table[i, j] = diffs.mean()  # averaged difference
    return rel_table


def get_valid_sites(virus: str, csv_pth: str):
    cal = BindingCalculator(csv_or_url=csv_pth, eliciting_virus=virus,
                            mutation_escape_strength=10)
    return set(cal.sites)


def get_ebs(virus: str, csv_pth: str, sites: list):
    cal = BindingCalculator(csv_or_url=csv_pth, eliciting_virus=virus,
                            mutation_escape_strength=10)
    ebs = list(map(lambda site: cal.binding_retained([site]), sites))
    return ebs


def main(save_path: str):
    assert save_path.endswith(".csv")
    escape_data_csv = "../data/escape_calculator_data.csv"
    edge_muts = np.load("../results/neutral_net_edgemuts.npy",
                        allow_pickle=True)
    edge_sites = set(list(map(lambda mut_code: int(mut_code[1:-1]), edge_muts)))
    # antibodies = ['BA.2.75', 'BQ.1.1', 'BA.2', 'XBB', 'BA.1', 'BA.5', 'D614G']
    viruses = ['pre-Omicron SARS-CoV-2', 'SARS-CoV-1 then SARS-CoV-2',
               'pre-Omicron SARS-CoV-2 then BA.2',
               'pre-Omicron SARS-CoV-2 then BA.5', 'SARS-CoV-2',
               'pre-Omicron SARS-CoV-2 then BA.1']
    # use only sites with data
    for vir in viruses:
        valid_sites = get_valid_sites(vir, escape_data_csv)
        edge_sites = edge_sites & valid_sites
    edge_sites = list(edge_sites)
    # create 1 table for each site mutation
    print(edge_sites)
    table = {}
    for vir in viruses:
        table[vir] = get_ebs(vir, escape_data_csv, edge_sites)
    print(table)
    df = pd.DataFrame(table)
    df['antigen'] = edge_sites
    df = df.set_index('antigen')
    r_table = calculate_difs(df, viruses, edge_sites)
    print(r_table)
    return table


if __name__ == "__main__":
    save_pth = "../results/escape_distance_table.csv"
    main(save_pth)
