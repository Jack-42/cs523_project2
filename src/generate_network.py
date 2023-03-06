import plotly.graph_objects as go
import networkx as nx
import numpy as np

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
    return mutated_seq


def main(base_seq: str):
    # in progress
    g = nx.Graph()
    n_muts_per = 4  # number of new nodes we attach to each node in perimeter
    n_variants_away = 4  # maximum number of mutations from source before stopping
    base_node = {"seq": base_seq, "n_muts": 0}
    g.add_node(base_node)
    cur_node = base_node
    excluded_pos = []
    for i in range(n_variants_away):
        for j in range(n_muts_per):
            mut_seq, chosen_positions = generate_mut(base_seq, i, excluded_pos)
            excluded_pos = excluded_pos + chosen_positions

def test():
    G = nx.random_geometric_graph(25, 0.125)
    edge_x = []
    edge_y = []
    for edge in G.edges():
        x0, y0 = G.nodes[edge[0]]['pos']
        x1, y1 = G.nodes[edge[1]]['pos']
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

    node_x = []
    node_y = []
    for node in G.nodes():
        x, y = G.nodes[node]['pos']
        node_x.append(x)
        node_y.append(y)

    node_trace = go.Scatter(
        x=node_x, y=node_y,
        mode='markers',
        hoverinfo='text',
        marker=dict(
            showscale=True,
            # colorscale options
            # 'Greys' | 'YlGnBu' | 'Greens' | 'YlOrRd' | 'Bluered' | 'RdBu' |
            # 'Reds' | 'Blues' | 'Picnic' | 'Rainbow' | 'Portland' | 'Jet' |
            # 'Hot' | 'Blackbody' | 'Earth' | 'Electric' | 'Viridis' |
            colorscale='YlGnBu',
            reversescale=True,
            color=[],
            size=10,
            colorbar=dict(
                thickness=15,
                title='Node Connections',
                xanchor='left',
                titleside='right'
            ),
            line_width=2))

    node_adjacencies = []
    node_text = []
    for node, adjacencies in enumerate(G.adjacency()):
        node_adjacencies.append(len(adjacencies[1]))
        node_text.append('# of connections: ' + str(len(adjacencies[1])))

    node_trace.marker.color = node_adjacencies
    node_trace.text = node_text

    fig = go.Figure(data=[edge_trace, node_trace],
                    layout=go.Layout(
                        title='<br>Network graph made with Python',
                        titlefont_size=16,
                        showlegend=False,
                        hovermode='closest',
                        margin=dict(b=20, l=5, r=5, t=40),
                        annotations=[dict(
                            text="Python code: <a href='https://plotly.com/ipython-notebooks/network-graphs/'> https://plotly.com/ipython-notebooks/network-graphs/</a>",
                            showarrow=False,
                            xref="paper", yref="paper",
                            x=0.005, y=-0.002)],
                        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False))
                    )
    fig.show()


if __name__ == "__main__":
    test()
