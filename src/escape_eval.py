import numpy as np
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


def main():
    edge_muts = np.load("../results/neutral_net_edgemuts.npy",
                        allow_pickle=True)
    binding_cal = BindingCalculator()  # use default settings
    edge_sites = list(map(lambda mut_code: int(mut_code[1:-1]), edge_muts))
    edge_sites = set(edge_sites) & set(binding_cal.sites)
    edge_bindings = list(
        map(lambda site: binding_cal.binding_retained([site]), edge_sites))
    binding_mp = dict(zip(edge_sites, edge_bindings))
    print(binding_mp)


if __name__ == "__main__":
    main()
