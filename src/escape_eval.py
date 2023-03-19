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


def get_valid_sites(antibody: str, csv_pth: str):
    cal = BindingCalculator(csv_or_url=csv_pth, eliciting_virus="SARS-CoV-2",
                            known_to_neutralize=antibody)
    return set(cal.sites)


def get_ebs(antibody: str, csv_pth: str, sites: list):
    cal = BindingCalculator(csv_or_url=csv_pth, eliciting_virus="SARS-CoV-2",
                            known_to_neutralize=antibody)
    ebs = list(map(lambda site: cal.binding_retained([site]), sites))
    return ebs


def main(save_path: str):
    assert save_path.endswith(".csv")
    escape_data_csv = "../data/escape_calculator_data.csv"
    edge_muts = np.load("../results/neutral_net_edgemuts.npy",
                        allow_pickle=True)
    edge_sites = set(list(map(lambda mut_code: int(mut_code[1:-1]), edge_muts)))
    # antibodies = ['BA.2.75', 'BQ.1.1', 'BA.2', 'XBB', 'BA.1', 'BA.5', 'D614G']
    antibodies = ['BA.2.75',
                  'BQ.1.1']  # = antiserum, sites are various antigens
    # use only sites with data
    for ab in antibodies:
        valid_sites = get_valid_sites(ab, escape_data_csv)
        edge_sites = edge_sites & valid_sites
    edge_sites = list(edge_sites)
    # create table
    print(edge_sites)
    table = {}
    for ab in antibodies:
        table[ab] = get_ebs(ab, escape_data_csv, edge_sites)
    print(table)
    df = pd.DataFrame(table)
    df['antigen'] = edge_sites
    df = df.set_index('antigen')
    print(df)
    return table


if __name__ == "__main__":
    main("../results/escape_distance_table.csv")
