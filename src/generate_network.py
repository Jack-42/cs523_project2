import plotly.graph_objects as go
import networkx as nx
import numpy as np
import queue

"""
@author Jack Ringer
Date: 3/5/2023
Description:
    Script used to generate neutral network and to display a random
    walk from the network (part 2).
"""


def generate_mut(seq: str, k: int, exclude_pos=[]):
    """
    Generate a mutated sequence from the input.
    :param seq: str, the base sequence
    :param k: int, the number of mutations to apply
    :param exclude_pos: list(optional), positions to not apply
        any mutation to
    :return: str, the mutated sequence
    """
    nts = ["A", "C", "G", "T"]
    r_indices = np.random.randint(0, len(seq), size=k)
    mutated_seq = list(seq)
    positions_used = []
    for i in range(len(r_indices)):
        r_idx = r_indices[i]
        while r_idx in r_indices[:i] or r_idx in exclude_pos:
            r_idx = np.random.randint(0, len(seq))
        base_nt = seq[r_idx]
        mut_nts = [x for x in nts if x != base_nt]
        mut_nt = np.random.choice(mut_nts)
        mutated_seq[r_idx] = mut_nt
        positions_used.append(r_idx)
    mutated_seq = ''.join(mutated_seq)  # back to string
    return mutated_seq, positions_used


def generate_network(base_seq: str, n_muts_per: int, n_variants_away: int):
    """
    Generate a neutral network from the given parameters using Plotly and
    networkx
    :param base_seq: str, the original genome sequence
    :param n_muts_per: int, number of new nodes we attach to each node in
        perimeter
    :param n_variants_away: int,  max number of mutations from source before
        stopping
    :return: None
    """
    g = nx.Graph()
    # node[0] = seq, node[1] = n_muts away
    base_node = (base_seq, 0)
    q = queue.Queue()
    q.put(base_node)
    while not (q.empty()):
        cur_node = q.get()
        seq = cur_node[0]
        cur_muts = cur_node[1]
        i = cur_muts + 1
        if i > n_variants_away:
            break
        excluded_pos = []
        for j in range(n_muts_per):
            # 1 mutation away from current seq == i+1 away from base
            mut_seq, chosen_positions = generate_mut(seq, 1, excluded_pos)
            excluded_pos = excluded_pos + chosen_positions
            mut_node = (mut_seq, i)
            g.add_edge(cur_node, mut_node)
            q.put(mut_node)

    pos_ = nx.spring_layout(g)
    seqs = [x[0] for x in g.nodes]
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


if __name__ == "__main__":
    generate_network("ACGTACGGTAAG", 3, 4)
