"""Calculate residual antibody binding after some mutations.

This module can be downloaded from
`https://github.com/jbloomlab/SARS2_RBD_Ab_escape_maps/blob/main/bindingcalculator.py <https://github.com/jbloomlab/SARS2_RBD_Ab_escape_maps/blob/main/bindingcalculator.py>`_

The module defines :class:`BindingCalculator` which does the calculation.

Written by Jesse Bloom.

"""

__docformat__ = 'numpy'

import pandas as pd


class BindingCalculator:
    """Calculates residual polyclonal antibody binding after some mutations.

    The calculator is the one implemented interactively at
    `https://jbloomlab.github.io/SARS2_RBD_Ab_escape_maps/escape-calc/ <https://jbloomlab.github.io/SARS2_RBD_Ab_escape_maps/escape-calc/>`_

    Parameters
    ----------
    csv_or_url : str
        Path to CSV or URL of CSV containing the escape data. Should
        have columns 'condition', 'metric', and 'escape'.
    eliciting_virus : str
        Include antibodies elicited by these viruses.
    known_to_neutralize : str
        Include antibodies known to neutralize this virus.
    weight_by_log_IC50 : bool
        Weight antibodies by log IC50.
    mutation_escape_strength : float
        Scaling exponent :math:`s`; larger values mean stronger escape, see
        https://jbloomlab.github.io/SARS2_RBD_Ab_escape_maps/escape-calc/

    Attributes
    ----------
    escape_data : pandas.DataFrame
        The data frame read from `csv_or_url` after filtering for specified
        `eliciting_virus` and `known_to_neutralize`.
    sites : set
        All sites for which we have escape data. We can only calculate effects
        of mutations at these sites.
    weight_by_log_IC50 : bool
        Value of `weight_by_log_IC50` passed as parameter.
    """

    def __init__(self,
                 csv_or_url='https://raw.githubusercontent.com/jbloomlab/SARS2_RBD_Ab_escape_maps/main/processed_data/escape_calculator_data.csv',
                 *,
                 eliciting_virus='SARS-CoV-2',
                 known_to_neutralize="any",
                 weight_by_log_IC50=True,
                 mutation_escape_strength=2,
                 ):
        """See main class docstring."""
        # read escape data
        self.escape_data = (
            pd.read_csv(csv_or_url)
            .assign(
                eliciting_virus=lambda x: x["eliciting_virus"].str.split(";"),
                known_to_neutralize=lambda x: x[
                    "known_to_neutralize"].str.split(";"),
                neg_log_IC50=lambda x: x["neg_log_IC50"].map(
                    lambda s: tuple(
                        [pd.NA if si == "NA" else float(si) for si in
                         s.split(";")])
                )
            )
            .explode("eliciting_virus")
            .explode(["known_to_neutralize", "neg_log_IC50"])
        )
        assert all(self.escape_data["neg_log_IC50"] >= 0)

        # make sure escape data has expected columns
        if not set(self.escape_data.columns).issuperset({'condition',
                                                         'site',
                                                         'escape',
                                                         'eliciting_virus',
                                                         "known_to_neutralize",
                                                         "neg_log_IC50",
                                                         }):
            raise ValueError(
                f"{self.escape_data.columns} lacks expected columns")

        # filter by virus
        eliciting_viruses = set(self.escape_data["eliciting_virus"])
        if eliciting_virus not in eliciting_viruses:
            raise ValueError(f"{eliciting_virus} not in {eliciting_viruses}")
        self.escape_data = self.escape_data.query(
            'eliciting_virus == @eliciting_virus').drop(
            columns="eliciting_virus"
        )
        assert len(self.escape_data) == len(self.escape_data.drop_duplicates())

        # filter by known_to_neutralize
        if known_to_neutralize not in set(
                self.escape_data['known_to_neutralize']):
            raise ValueError(f"invalid {known_to_neutralize}")
        self.escape_data = (
            self.escape_data
            .query("known_to_neutralize == @known_to_neutralize")
            .drop(columns="known_to_neutralize")
        )
        assert len(self.escape_data) == len(self.escape_data.drop_duplicates())

        # get escape scaled by the max escape for that condition
        self.escape_data = (
            self.escape_data
            .assign(max_escape=lambda x: (x.groupby('condition')
                                          ['escape']
                                          .transform('max')
                                          ),
                    scale_escape=lambda x: x['escape'] / x['max_escape'],
                    )
        )

        # get all sites for which we have escape data
        self.sites = set(self.escape_data['site'])

        # set mutation escape strength
        self.mutation_escape_strength = mutation_escape_strength

        # do we weight by log IC50?
        self.weight_by_log_IC50 = weight_by_log_IC50

        assert (
                self.escape_data["condition"].nunique()
                == len(
            self.escape_data[["condition", "neg_log_IC50"]].drop_duplicates())
        )

        # number of conditions (antibodies), weighting by negative log IC50 if doing that
        if self.weight_by_log_IC50:
            self._n_conditions = (
                self.escape_data
                [["condition", "neg_log_IC50"]]
                .drop_duplicates()
                ["neg_log_IC50"]
                .sum()
            )
        else:
            self._n_conditions = self.escape_data['condition'].nunique()

    def escape_per_site(self, mutated_sites):
        """Escape at each site after mutating indicated sites.

        Parameters
        ----------
        mutated_sites : array-like of integers
            List of mutated sites, must all be in :attr:`BindingCalculator.sites`.

        Returns
        -------
        pandas.DataFrame
            For each site, gives the original escape and the escape
            retained after mutations.

        """
        mutated_sites = set(mutated_sites)
        if not mutated_sites.issubset(self.sites):
            raise ValueError(f"invalid sites: {mutated_sites - self.sites}")
        df = (
                 self.escape_data
                 .assign(
                     mutated=lambda x: x['site'].isin(mutated_sites).astype(
                         int),
                     site_bind_retain=lambda x: 1 - x['scale_escape'] * x[
                         'mutated']
                 )
                 .groupby(['condition', "neg_log_IC50"], as_index=False)
                 .aggregate(
                     cond_bind_retain=pd.NamedAgg('site_bind_retain', "prod"))
                 .assign(
                     cond_bind_retain=lambda x: x["cond_bind_retain"].pow(
                         self.mutation_escape_strength),
                     weight=lambda x: x[
                         "neg_log_IC50"] if self.weight_by_log_IC50 else 1,
                 )
                 [["condition", 'cond_bind_retain', "weight"]]
                 .merge(self.escape_data[['condition', 'site', 'escape']])
                 .assign(
                     escape=lambda x: x["escape"] * x["weight"],
                     retained_escape=lambda x: x['cond_bind_retain'] * x[
                         'escape'],
                 )
                 .groupby('site')
                 .aggregate(
                     original_escape=pd.NamedAgg('escape', 'sum'),
                     retained_escape=pd.NamedAgg('retained_escape', 'sum'),
                 )
             ) / self._n_conditions
        return df.reset_index()

    def binding_retained(self, mutated_sites):
        """Fraction binding retained after mutating indicated sites.

        Parameters
        ----------
        mutated_sites : array-like of integers
            List of mutated sites, must all be in :attr:`BindingCalculator.sites`.

        Returns
        -------
        float
            The fraction binding retained after these mutations.

        """
        mutated_sites = set(mutated_sites)
        if not mutated_sites.issubset(self.sites):
            raise ValueError(f"invalid sites: {mutated_sites - self.sites}")
        binding_retained = (
                               self.escape_data
                               .assign(
                                   mutated=lambda x: x['site'].isin(
                                       mutated_sites).astype(int),
                                   site_bind_retain=lambda x: 1 - x[
                                       'scale_escape'] * x['mutated'],
                               )
                               .groupby(['condition', 'neg_log_IC50'],
                                        as_index=False)
                               .aggregate(cond_bind_retain=pd.NamedAgg(
                                   'site_bind_retain', "prod"))
                               .assign(
                                   cond_bind_retain=lambda x: x[
                                       "cond_bind_retain"].pow(
                                       self.mutation_escape_strength),
                                   weight=lambda x: x[
                                       "neg_log_IC50"] if self.weight_by_log_IC50 else 1,
                                   weighted_cond_bind_retain=lambda x: x[
                                                                           "cond_bind_retain"] *
                                                                       x[
                                                                           "weight"],
                               )
                               ['weighted_cond_bind_retain']
                               .sum()
                           ) / self._n_conditions
        return binding_retained


if __name__ == '__main__':
    binding_cal = BindingCalculator()
    #print(binding_cal.binding_retained([420]))
    sites = []
    for site in binding_cal.escape_data['site']:
        sites.append(site)

    sites = [*set(sites)]
    print(sites)

    num = 0
    print(len(sites))
    for s in sites:
        escape = binding_cal.binding_retained([s])
        if (escape >= .94):
            num += 1

    print(num)