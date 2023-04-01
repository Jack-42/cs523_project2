# CS 523 Project 2

By James Morris, Carter Frost, Dominic Larranaga, and Jack Ringer

## Overview

The data, code, and results contained within this repo were used to generate
our report for project 2 of CS 523.

### data dir

The data directory contains data used as input for various parts of the project. Some of these files were generated by
our own code, whereas others were taken from external sources.
The filenames and sources for each are provided by the table below.

| Filename                       |                                                       Source                                                        | 
|:-------------------------------|:-------------------------------------------------------------------------------------------------------------------:|
| aa_fitness.csv                 |                     https://github.com/jbloomlab/SARS2-mut-fitness/tree/main/results/aa_fitness                     |
| codons.txt                     |                                https://github.com/zhanxw/anno/blob/master/codon.txt                                 |
| cov_spike_nt_genome.fasta      |                   generated from covid_spike_protein.fasta (see acid_seq_to_genome() in utils.py)                   |
| covid_rbd.fasta                |                                   https://www.ncbi.nlm.nih.gov/nuccore/MN908947.3                                   |
| covid_rbd_freqs.csv            |                                       generated (see part1a_calculations.py)                                        |
| covid_spike_protein.fasta      |                                   https://www.uniprot.org/uniprotkb/P0DTC2/entry                                    |
| epistatic_genome_frequency.csv |                                                      generated                                                      |
| escape_calculator_data.csv     | https://raw.githubusercontent.com/jbloomlab/SARS2_RBD_Ab_escape_maps/main/processed_data/escape_calculator_data.csv |

### src dir

The src directory contains scripts used to generate results for various parts of the project. Some of the code was taken
from
or heavily based on other sources. We document these cases below.

| Filename               |                                  Source/Based Upon                                   | 
|:-----------------------|:------------------------------------------------------------------------------------:|
| binding_calculator.py  | https://github.com/jbloomlab/SARS2_RBD_Ab_escape_maps/blob/main/bindingcalculator.py |
| create_antigenic_map.R |       https://acorg.github.io/Racmacs/articles/making-a-map-from-scratch.html        |

### results dir

The results directory contains output files (figures, tables, etc.) related to results discussed in
the project report.

## Setup

All necessary packages to run the code in this project can be installed with:

pip install -r requirements.txt

## Acknowledgement

The code for the binding calculator (src/binding_calculator.py) was taken
from the Bloom lab GitHub page:
https://github.com/jbloomlab/SARS2_RBD_Ab_escape_maps/blob/main/bindingcalculator.py

We also acknowledge that much of our data was taken from third parties, please see the first table shown above.