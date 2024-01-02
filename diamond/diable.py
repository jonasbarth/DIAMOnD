"""Module containing code for the DIAMOnD Background Local Expansion (DiaBLE) algorithm."""
from typing import Sequence

import networkx as nx
import pandas as pd

from .DIAMOnD import read_input, diamond_iteration_of_first_X_nodes, get_neighbors_and_degrees, compute_all_gamma_ln, \
    reduce_not_in_cluster_nodes, pvalue


def diable(network_file: str, seed_genes_file: str, num_genes_to_add: int):
    """Runs the DiaBLE algorithm on the provided network and seed genes.

    The DiaBLE algorithm is an iterative variant of DIAMOnD which considers a growing gene universe instead of
    always working with the entire network. At each DiaBLE iteration, the entire gene graph consists out of

    - current seed genes
    - candidate genes (genes with => 1 connection to current seed genes)
    - candidate gene neighbours

    Like in DIAMOnD, we then calculate the p-value of each gene in the set of candidate genes and their neighbours,
    adding the gene with the lowest p-value to the set of seed genes.

    Args:
        network_file (str): path to the network file, containing an edge list of nodes.
        seed_genes_file (str): path to the seed genes file, containing a list of seed genes.
        num_genes_to_add (int): the number of genes to add to the network.

    Return:
        None
    """
    G_original, seed_genes = read_input(network_file, seed_genes_file)

    # 1. throwing away the seed genes that are not in the network
    all_genes_in_network = set(G_original.nodes())
    seed_genes = set(seed_genes)
    disease_genes = seed_genes & all_genes_in_network

    if len(disease_genes) != len(seed_genes):
        print("DIAMOnD(): ignoring %s of %s seed genes that are not in the network" % (
            len(seed_genes - all_genes_in_network), len(seed_genes)))

    added_genes = []

    universe = create_diable_universe(G_original, disease_genes)

    while len(added_genes) < num_genes_to_add:
        new_gene, *_, p_value = diamond_iteration(universe, disease_genes, 1, 1)

        disease_genes.add(new_gene)
        added_genes.append(new_gene)

        universe = create_diable_universe(G_original, disease_genes)
        print(f"Universe size: {len(universe)}")

def find_nodes_with_links_to(graph: nx.Graph, nodes: Sequence[int]):
    """Finds all nodes in the graph that are directly linked to the given nodes.

    Args:
        graph (nx.Graph): the graph in which to find the nodes linking to the set of nodes.
        nodes (Sequence[int]): the nodes of which to find neighbours.

    Return:
        set: a list of nodes that are connected to the given nodes.
    """
    nodes = set(nodes)
    linked = set()
    for node in nodes:
        neighbours = set(nx.neighbors(graph, node))
        linked |= neighbours

    return linked - nodes


def create_diable_universe(network: nx.Graph, seed_genes: Sequence[int]):
    """Creates a diable universe from the given network and seed genes.

    The diable universe consists out of:

    - seed genes
    - candidate genes, which are seed gene neighbours
    - candidate gene neighbours.

    Args:
        network (nx.Graph): the entire gene network.
        seed_genes (Sequence[int]): a sequence of seed genes.

    Return:
        diable_universe (nx.Graph): a graph that is the diable universe.
    """
    G_universe = network.subgraph(seed_genes)

    # find the neighbours of the disease genes
    seed_neighbours = network.subgraph(find_nodes_with_links_to(network, seed_genes))

    # find the neighbours of the neighbours of the disease genes
    seed_neighbours_neighbours = network.subgraph(find_nodes_with_links_to(network, seed_neighbours.nodes()))

    return nx.compose(nx.compose(G_universe, seed_neighbours), seed_neighbours_neighbours)


def diamond_iteration(universe, seed_genes, X, alpha):
    """
    Parameters:
    ----------
    - G:     graph
    - S:     seeds
    - X:     the number of iterations, i.e only the first X gened will be
             pulled in
    - alpha: seeds weight
    Returns:
    --------

    - added_nodes: ordered list of nodes in the order by which they
      are agglomerated. Each entry has 4 info:
      * name : dito
      * k    : degree of the node
      * kb   : number of +1 neighbors
      * p    : p-value at agglomeration
    """
    size_universe = universe.number_of_nodes()

    # ------------------------------------------------------------------
    # Setting up initial set of nodes in cluster
    # ------------------------------------------------------------------

    cluster_nodes = set(seed_genes)
    s0 = len(cluster_nodes)

    s0 += (alpha - 1) * s0
    size_universe += (alpha - 1) * s0

    # ------------------------------------------------------------------
    # precompute the logarithmic gamma functions
    # ------------------------------------------------------------------
    gamma_ln = compute_all_gamma_ln(size_universe + 1)

    # ------------------------------------------------------------------
    # Setting initial set of nodes not in cluster
    # ------------------------------------------------------------------
    not_in_cluster = nx.Graph(universe.edges() - cluster_nodes).nodes()

    # ------------------------------------------------------------------
    #
    # M A I N     L O O P
    #
    # ------------------------------------------------------------------
    info = []
    for node in not_in_cluster:
        degree = nx.degree(universe, node)
        num_links_to_seed_genes = len(set(nx.neighbors(universe, node)) & set(seed_genes))
        p = pvalue(num_links_to_seed_genes, degree, size_universe, s0, gamma_ln)

        info.append([node, degree, num_links_to_seed_genes, p])

    candidate_genes = pd.DataFrame(info, columns=['node', 'degree', 'num_links_to_seed_genes', 'p_value'])
    candidate_genes.sort_values("p_value", ascending=True, inplace=True)
    return (*candidate_genes.head(1).iloc[0].tolist(),)
