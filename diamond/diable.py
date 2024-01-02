"""Module containing code for the DIAMOnD Background Local Expansion (DiaBLE) algorithm."""
from typing import Sequence

import networkx as nx

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
        new_node = diamond_iteration(universe, disease_genes, 1, 1)

        new_gene = new_node[0][0]
        disease_genes.add(new_gene)
        added_genes += new_node

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

    N = universe.number_of_nodes()

    added_nodes = []

    # ------------------------------------------------------------------
    # Setting up dictionaries with all neighbor lists
    # and all degrees
    # ------------------------------------------------------------------
    neighbors, all_degrees = get_neighbors_and_degrees(universe)

    # ------------------------------------------------------------------
    # Setting up initial set of nodes in cluster
    # ------------------------------------------------------------------

    cluster_nodes = set(seed_genes)
    not_in_cluster = set()
    s0 = len(cluster_nodes)

    s0 += (alpha - 1) * s0
    N += (alpha - 1) * s0

    # ------------------------------------------------------------------
    # precompute the logarithmic gamma functions
    # ------------------------------------------------------------------
    gamma_ln = compute_all_gamma_ln(N + 1)

    # ------------------------------------------------------------------
    # Setting initial set of nodes not in cluster
    # ------------------------------------------------------------------
    not_in_cluster = nx.Graph(universe.edges() - cluster_nodes).nodes()

    # ------------------------------------------------------------------
    #
    # M A I N     L O O P
    #
    # ------------------------------------------------------------------

    all_p = {}

    while len(added_nodes) < X:

        # ------------------------------------------------------------------
        #
        # Going through all nodes that are not in the cluster yet and
        # record k, kb and p
        #
        # ------------------------------------------------------------------

        info = {}

        pmin = 10
        next_node = 'nix'

        # for DiaBLE, I want to make the universe:
        #   seed genes + candidate genes (=> 1 link to seed genes) + candidate gene neighbours
        # the nodes to considered to be added to the seeed genes are just the candidate genes + their neighbours
        reduced_not_in_cluster = reduce_not_in_cluster_nodes(all_degrees,
                                                             neighbors, universe,
                                                             not_in_cluster,
                                                             cluster_nodes, alpha)

        # go through the entire universe, calculate p-values of each node not already in interesting genes
        # and finally add the node with the lowest p-value to the set of interesting genes
        for node, kbk in reduced_not_in_cluster.items():
            # Getting the p-value of this kb,k
            # combination and save it in all_p, so computing it only once!
            kb, k = kbk
            try:
                p = all_p[(k, kb, s0)]
            except KeyError:
                p = pvalue(kb, k, N, s0, gamma_ln)
                all_p[(k, kb, s0)] = p

            # recording the node with smallest p-value
            if p < pmin:
                pmin = p
                next_node = node

            info[node] = (k, kb, p)

        # ---------------------------------------------------------------------
        # Adding node with smallest p-value to the list of aaglomerated nodes
        # ---------------------------------------------------------------------
        added_nodes.append((next_node,
                            info[next_node][0],
                            info[next_node][1],
                            info[next_node][2]))

        # Updating the list of cluster nodes and s0
        cluster_nodes.add(next_node)
        s0 = len(cluster_nodes)
        not_in_cluster |= (neighbors[next_node] - cluster_nodes)
        not_in_cluster.remove(next_node)

    return added_nodes
