import networkx
import pandas as pd
import itertools
import plotly.graph_objects as go
import networkx as nx
import numpy as np
import queue

from utils import load_fasta_seq

"""
@author Jack Ringer
Date: 3/5/2023
Description:
    Script used to generate neutral network and to display a random
    walk from the network (part 2).
"""


def generate_synom_df() -> pd.DataFrame:
    """
    Generate Pandas Dataframe containing all synonymous mutations.
    Will have 3 columns, "BaseCodon", "Mutation", and "Position"
    :return: pd.DataFrame containing synom mutations (codon-wise)
    """
    df = pd.read_csv("../data/codons.txt", sep='\t')
    df = df.rename({'    AminoAcid': "AminoAcid"}, axis="columns")
    base_codons = []
    mut_nts = []
    positions = []
    nucleotides = ["A", "C", "G", "T"]
    for base, mutant in itertools.product(nucleotides, nucleotides):
        if base == mutant:
            continue
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
                    base_codon = cdn[:i] + base + cdn[i + 1:]
                    base_codons.append(base_codon)
                    mut_nts.append(mutant)
                    positions.append(i)
    mp_df = pd.DataFrame({"BaseCodon": base_codons, "Mutation": mut_nts, "Position": positions})
    # remove duplicates
    mp_df = mp_df.drop_duplicates()
    return mp_df


def get_codon(seq: str, i: int):
    """
    Get the 3-letter codon from a given position
    :param seq: str, the full genome
    :param i: int, the position
    :return: str, the codon
    """
    codon_index = i - (i % 3)
    return seq[codon_index: codon_index + 3]


def generate_mut(seq: str, exclude_pos=[]):
    """
    Generate a mutated sequence from the input.
    :param seq: str, the base sequence
    :param exclude_pos: list(optional), positions to not apply
        any mutation to
    :return: str, the mutated sequence
    """
    nts = ["A", "C", "G", "T"]
    r_idx = np.random.randint(0, len(seq))
    mutated_seq = list(seq)
    while r_idx in exclude_pos:
        # this is a bad way to prevent the same position from being mutated
        # watch out for infinite loops
        r_idx = np.random.randint(0, len(seq))
    base_nt = seq[r_idx]
    mut_nts = [x for x in nts if x != base_nt]
    mut_nt = np.random.choice(mut_nts)
    mutated_seq[r_idx] = mut_nt
    mutated_seq = ''.join(mutated_seq)  # back to string
    return mutated_seq, r_idx


def generate_synom_mut(seq: str, synom_df: pd.DataFrame, exclude_pos=None):
    """
    Generate a synonymous mutation for the given sequence.
    :param seq: str, genome
    :param synom_df: pd.DataFrame, output of generate_synom_df()
    :param exclude_pos: list (optional), positions mutation not allowed for
    :return: str, mutated (but synonymous) genome
    """
    # about 1/4 of mutations are synonymous
    if exclude_pos is None:
        exclude_pos = []
    mutated_seq = list(seq)
    r_idx = np.random.randint(0, len(seq))
    base_cdn = get_codon(seq, r_idx)
    if r_idx in exclude_pos:
        # this is a bad way to prevent the same position from being mutated
        # watch out for infinite loops
        return generate_synom_mut(seq, synom_df, exclude_pos)
    sub_df = synom_df.loc[(synom_df["BaseCodon"] == base_cdn) & (synom_df["Position"] == (r_idx % 3))]
    if len(sub_df) == 0:
        # watch out for infinite loops
        return generate_synom_mut(seq, synom_df, exclude_pos)
    mut_nt = sub_df["Mutation"].sample(n=1).iloc[0]
    mutated_seq[r_idx] = mut_nt
    mutated_seq = ''.join(mutated_seq)  # back to string
    mut_code = "%s%d%s" % (seq[r_idx], r_idx, mut_nt)
    return mutated_seq, r_idx, mut_code


def generate_neutral_net(base_seq: str, n_jumps: int, nodes_per_jump: int):
    """
    Generate a neutral network with a given number of jumps
    :param base_seq: str, the base genome
    :param n_jumps: int, number of jumps
    :param nodes_per_jump: int, mutations for each jump
    :return: nx.Graph object representing the neutral network
    """
    net = nx.Graph()
    synom_df = generate_synom_df()
    base_node = (base_seq, 0)
    q = queue.Queue()
    q.put(base_node)
    while not (q.empty()):
        cur_node = q.get()
        seq = cur_node[0]
        cur_muts = cur_node[1]
        i = cur_muts + 1
        if i > n_jumps:
            break
        excluded_pos = []
        for _ in range(nodes_per_jump):
            mut_seq, idx, mut_code = generate_synom_mut(seq, synom_df, excluded_pos)
            mut_node = (mut_seq, i, mut_code)
            net.add_edge(cur_node, mut_node)
            q.put(mut_node)
            excluded_pos.append(idx)
    return net


def show_neutral_network(base_seq: str, n_jumps: int, nodes_per_jump: int):
    """
    Generate a neutral network from the given parameters using Plotly and
    networkx
    :param base_seq: str, the original genome sequence
    :param n_jumps: int, number of jumps
    :param nodes_per_jump: int, mutations for each jump
    :return: None
    """
    g = generate_neutral_net(base_seq, n_jumps, nodes_per_jump)
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
                        title='Neutral Network',
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


def get_edge_mutations(net: networkx.Graph, n_jumps: int, mut_per_node: int):
    """
    Get a list of mutations generated at the edge of the neutral network.
    :param net: networkx.Graph, neutral network
    :param n_jumps: int, the number of jumps we performed in the network
    :param mut_per_node: int, the number of mutated sequences to generate from each edge node
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
        cur_muts = cur_node[1]
        i = cur_muts + 1
        excluded_pos = []
        for _ in range(mut_per_node):
            mut_seq, idx = generate_mut(seq, excluded_pos)
            mut_node = (mut_seq, i)
            edge_muts.append(mut_node)
            excluded_pos.append(idx)
    return edge_muts


if __name__ == "__main__":
    seq1 = load_fasta_seq("../data/covid_spike_protein.fasta")
    n_jumps1 = 4
    nodes_per_jump1 = 2
    neutral_net = show_neutral_network(seq1, n_jumps1, nodes_per_jump1)
