import queue

import networkx
import networkx as nx
import numpy as np
import pandas as pd
import plotly.graph_objects as go

from utils import load_fasta_seq, load_codon_whl

"""
@author Jack Ringer
Date: 3/5/2023
Description:
    Script used to generate and show a neutral network (part 2).
"""


def get_codon(seq: str, i: int):
    """
    Get the 3-letter codon from a given position
    :param seq: str, the full genome
    :param i: int, the position
    :return: str, the codon
    """
    codon_index = i - (i % 3)
    return seq[codon_index: codon_index + 3]


def generate_mut(seq: str, codon_df: pd.DataFrame, mut_type: str = "any",
                 exclude_pos=None, chosen_pos=None):
    """
    Generate a mutated sequence from the input.
    :param seq: str, the base sequence
    :param mut_type: str, one of "synonymous", "nonsynonymous", or "any"
    :param codon_df: pd.DataFrame, contains mapping of codons to amino acids
    :param exclude_pos: list(optional), positions to not apply
        any mutation to
    :param chosen_pos: list(optional), specific indices to pull mutations from.
    :return: str, the mutated sequence
    """
    if mut_type not in ["synonymous", "nonsynonymous", "any"]:
        raise ValueError(
            "argument 'mut_type' must be one of the following: 'synonymous', "
            "'nonsynonymous', 'any'")
    if exclude_pos is None:
        exclude_pos = []
    nts = ["A", "C", "G", "T"]
    if chosen_pos is None:
        r_idx = np.random.randint(0, len(seq))
    else:
        r_idx = np.random.choice(chosen_pos)
    mutated_seq = list(seq)
    n_tries = 0
    max_tries = 500
    while r_idx in exclude_pos and n_tries < max_tries:
        # this is a bad way to prevent the same position from being mutated
        # oh well
        if chosen_pos is None:
            r_idx = np.random.randint(0, len(seq))
        else:
            r_idx = np.random.choice(chosen_pos)
        n_tries += 1
    if n_tries >= max_tries:
        raise RuntimeError(
            "Exceeded max number of tries trying to find appropriate r_idx")
    base_nt = seq[r_idx]
    mut_nts = [x for x in nts if x != base_nt]
    mut_nt = np.random.choice(mut_nts)
    mutated_seq[r_idx] = mut_nt
    mutated_seq = ''.join(mutated_seq)  # back to string
    mut_code = get_mut_code(seq, mutated_seq, r_idx, codon_df)
    if mut_type == "any":
        return mutated_seq, r_idx, mut_code
    elif mut_type == "synonymous" and mut_code[0] == mut_code[-1]:
        return mutated_seq, r_idx, mut_code
    elif mut_type == "nonsynonymous" and mut_code[0] != mut_code[-1]:
        return mutated_seq, r_idx, mut_code
    # try again
    return generate_mut(seq, codon_df, mut_type, exclude_pos, chosen_pos)


def get_mut_code(base_seq: str, mut_seq: str, nt_idx: int,
                 codon_df: pd.DataFrame):
    """
    Generate mutation code (base aa, idx, mut aa) for a mutation.
    :param base_seq: str, the base sequence
    :param mut_seq: str, the mutated sequence
    :param nt_idx: int, the index where the mutation occurred
    :param codon_df: pd.DataFrame, mapping of codons to amino acids
    :return: str, the mutation code
    """
    base_cdn = get_codon(base_seq, nt_idx)
    base_acid = codon_df[codon_df["Codon"] == base_cdn]["Letter"].iloc[0]
    mut_cdn = get_codon(mut_seq, nt_idx)
    mut_acid = codon_df[codon_df["Codon"] == mut_cdn]["Letter"].iloc[0]
    acid_idx = (nt_idx - (nt_idx % 3)) / 3
    acid_idx += 1  # 1-indexing used by Bloom calculator
    mut_code = "%s%d%s" % (base_acid, acid_idx, mut_acid)
    return mut_code


def generate_neutral_net(base_seq: str, n_jumps: int, nodes_per_jump: int,
                         codon_df: pd.DataFrame):
    """
    Generate a neutral network with a given number of jumps
    :param base_seq: str, the base genome
    :param n_jumps: int, number of jumps
    :param nodes_per_jump: int, mutations for each jump
    :param codon_df: pd.DataFrame, contains mapping of codons to amino acids
    :return: nx.Graph object representing the neutral network
    """
    net = nx.Graph()
    base_node = (base_seq, 0)
    q = queue.Queue()
    q.put(base_node)
    excluded_pos = []
    while not (q.empty()):
        cur_node = q.get()
        seq = cur_node[0]
        cur_muts = cur_node[1]
        i = cur_muts + 1
        if i > n_jumps:
            break
        for _ in range(nodes_per_jump):
            # generate mutation 1-step away from cur_node
            mut_seq, idx, mut_code = generate_mut(seq, codon_df, "synonymous",
                                                  excluded_pos)
            mut_node = (mut_seq, i, mut_code)
            net.add_edge(cur_node, mut_node)
            q.put(mut_node)
            excluded_pos.append(idx)
    return net


def show_neutral_network(base_seq: str, n_jumps: int, nodes_per_jump: int,
                         codon_whl: pd.DataFrame):
    """
    Generate a neutral network from the given parameters using Plotly and
    networkx
    :param base_seq: str, the original genome sequence
    :param n_jumps: int, number of jumps
    :param nodes_per_jump: int, mutations for each jump
    :param codon_whl: pd.DataFrame, contains mapping of codons to amino acids
    :return: nx.Graph object representing the neutral network
    """
    g = generate_neutral_net(base_seq, n_jumps, nodes_per_jump, codon_whl)
    pos_ = nx.spring_layout(g)
    # For each edge, make an edge_trace, append to list
    edge_x = []
    edge_y = []
    for edge in g.edges():
        node1 = edge[0]
        node2 = edge[1]
        x0, y0 = pos_[node1]
        x1, y1 = pos_[node2]
        edge_x.append(x0)
        edge_x.append(x1)
        edge_x.append(None)
        edge_y.append(y0)
        edge_y.append(y1)
        edge_y.append(None)
    edge_trace = go.Scatter(
        x=edge_x, y=edge_y,
        line=dict(width=0.5, color='#888'),
        hoverinfo='none',
        mode='lines')

    # Make a node trace
    node_trace = go.Scatter(x=[],
                            y=[],
                            text=[],
                            textposition="top center",
                            mode='markers',
                            textfont_size=10,
                            marker=dict(
                                showscale=True,
                                # colorscale options
                                # 'Greys' | 'YlGnBu' | 'Greens' | 'YlOrRd' | 'Bluered' | 'RdBu' |
                                # 'Reds' | 'Blues' | 'Picnic' | 'Rainbow' | 'Portland' | 'Jet' |
                                # 'Hot' | 'Blackbody' | 'Earth' | 'Electric' | 'Viridis' |
                                colorscale='Rainbow',
                                reversescale=True,
                                color=[],
                                size=10,
                                colorbar=dict(
                                    thickness=15,
                                    title='Mutations',
                                    xanchor='left',
                                    titleside='right',
                                    dtick=1
                                ),
                                line_width=2))
    # For each node in g, get the position and size and add to the node_trace
    node_text = []
    node_colors = []
    for node in g.nodes():
        x, y = pos_[node]
        node_trace['x'] += tuple([x])
        node_trace['y'] += tuple([y])
        node_colors.append(node[1])
        node_text.append('# of mutations: ' + str(node[1]))

    node_trace.marker.color = node_colors
    node_trace.text = node_text

    fig = go.Figure(data=[edge_trace, node_trace],
                    layout=go.Layout(
                        title='Subspace of Neutral Network for COVID-19 Spike Protein',
                        titlefont_size=16,
                        showlegend=False,
                        hovermode='closest',
                        margin=dict(b=20, l=5, r=5, t=40),
                        xaxis=dict(showgrid=False, zeroline=False,
                                   showticklabels=False),
                        yaxis=dict(showgrid=False, zeroline=False,
                                   showticklabels=False))
                    )
    fig.show()
    return g


def get_edge_mutations(net: networkx.Graph, n_jumps: int,
                       mut_per_node: int, codon_df: pd.DataFrame,
                       idx_pool=None):
    """
    Get a list of mutations generated at the edge of the neutral network.
    :param net: networkx.Graph, neutral network
    :param n_jumps: int, the number of jumps we performed in the network
    :param mut_per_node: int, the number of mutated sequences to generate from
                        each edge node.
    :param codon_df: pd.DataFrame, dataframe containing codon wheel info
    :param idx_pool: list (optional), list of indices where mutations allowed
                    to occur. Useful for Bloom calculator since escape values
                    only available for a subset of positions.
    :return: list of edge mutations
    """
    # get edge nodes in neutral net
    edge_q = queue.Queue()
    for node in net.nodes():
        if node[1] == n_jumps:
            edge_q.put(node)
    # generate mutations
    edge_muts = []
    while not (edge_q.empty()):
        cur_node = edge_q.get()
        seq = cur_node[0]
        excluded_pos = []
        for _ in range(mut_per_node):
            # generate mutation 1-step away from cur_node
            mut_seq, idx, mut_code = generate_mut(seq, codon_df,
                                                  "nonsynonymous", excluded_pos,
                                                  idx_pool)
            edge_muts.append(mut_code)
            excluded_pos.append(idx)
    return edge_muts


def get_nt_sites(bloom_sites: np.array):
    sites = bloom_sites - 1
    nt_sites = []
    for aa_site in sites:
        nt_sites.append(aa_site * 3)
        nt_sites.append(aa_site * 3 + 1)
        nt_sites.append(aa_site * 3 + 2)
    return nt_sites


def main():
    # for reproducibility
    np.random.seed(43)

    # generate network for spike protein
    seq1 = load_fasta_seq("../data/cov_spike_nt_genome.fasta")
    codon_whl = load_codon_whl("../data/codons.txt")
    n_jumps1 = 4
    nodes_per_jump1 = 2
    neutral_net = show_neutral_network(seq1, n_jumps1, nodes_per_jump1,
                                       codon_whl)

    # generate mutations only where escape data is available
    escape_aa_sites = np.load("../results/bloom_valid_sites.npy",
                              allow_pickle=True)
    escape_nt_sites = get_nt_sites(escape_aa_sites)
    edge_muts = get_edge_mutations(neutral_net, n_jumps1,
                                   nodes_per_jump1,
                                   codon_whl, escape_nt_sites)

    # save data
    networkx.write_adjlist(neutral_net, "../results/neutral_net.adjlist")
    np.save("../results/neutral_net_edgemuts.npy", edge_muts)


if __name__ == "__main__":
    main()
    # for loading in data
    edge_muts1 = np.load("../results/neutral_net_edgemuts.npy",
                         allow_pickle=True)
    print(len(edge_muts1))
    print(edge_muts1)
    # neutral_net = networkx.read_adjlist("../results/neutral_net.adjlist")
