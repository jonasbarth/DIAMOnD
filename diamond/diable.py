"""Module containing code for the DIAMOnD Background Local Expansion (DiaBLE) algorithm."""
from typing import Sequence, Union, Dict, Iterable

import networkx as nx
import pandas as pd
from tqdm import tqdm
from scipy.stats import hypergeom

from .DIAMOnD import read_input, compute_all_gamma_ln, \
    pvalue


def diamond2(network_file: str, seed_genes_file: str, num_genes_to_add: int):
    G_original, seed_genes = read_input(network_file, seed_genes_file)

    # 1. throwing away the seed genes that are not in the network
    all_genes_in_network = set(G_original.nodes())
    seed_genes = set(seed_genes)
    disease_genes = seed_genes & all_genes_in_network

    if len(disease_genes) != len(seed_genes):
        print("DIAMOnD(): ignoring %s of %s seed genes that are not in the network" % (
            len(seed_genes - all_genes_in_network), len(seed_genes)))

    added_genes = []

    for i in range(num_genes_to_add):
        gene, degree, num_links_to_seed_genes, p_value = diamond_iteration(G_original, disease_genes)
        added_genes.append([gene, degree, num_links_to_seed_genes, p_value])

        disease_genes.add(gene)

    candidate_genes = pd.DataFrame(added_genes, columns=['node', 'degree', 'num_links_to_seed_genes', 'p_value'])
    return candidate_genes


def diable2(network_file: Union[str, pd.DataFrame], seed_genes_file: Union[str, pd.DataFrame], num_genes_to_add: int, ** kwargs):
    """Runs the DiaBLE algorithm on the provided network and seed genes.

    The DiaBLE algorithm is an iterative variant of DIAMOnD which considers a growing gene universe instead of
    always working with the entire network. At each DiaBLE iteration, the entire gene graph consists out of

    - current seed genes
    - candidate genes (genes with => 1 connection to current seed genes)
    - candidate gene neighbours

    Like in DIAMOnD, we then calculate the p-value of each gene in the set of candidate genes and their neighbours,
    adding the gene with the lowest p-value to the set of seed genes.

    Args:
        network_file (Union[str, pd.DataFrame]): path to the network file or a pandas DataFrame, containing an edge list
         of nodes.
        seed_genes_file (Union[str, pd.DataFrame]): path to the seed genes file or a pandas DataFrame, containing a list
        of seed genes.
        num_genes_to_add (int): the number of genes to add to the network.

    Return:
        (pd.DataFrame): a pandas dataframe with all the DiaBLE genes and their p-values.
    """
    if type(network_file) is not type(seed_genes_file):
        raise TypeError(
            f"Both the network and seed must be of the same type. Network: {type(network_file)}, Seed: {type(seed_genes_file)}.")

    if isinstance(network_file, pd.DataFrame):
        # ensure that genes are read as strings because when reading from .csv they are read as strings by default
        network_file.iloc[:, 0] = network_file.iloc[:, 0].astype(str)
        network_file.iloc[:, 1] = network_file.iloc[:, 1].astype(str)

        G_original = nx.Graph()
        G_original.add_edges_from(zip(network_file.iloc[:, 0], network_file.iloc[:, 1]))

    if isinstance(seed_genes_file, pd.DataFrame):
        seed_genes = set(map(str, seed_genes_file.iloc[:, 0]))

    else:
        G_original, seed_genes = read_input(network_file, seed_genes_file)

    # 1. throwing away the seed genes that are not in the network
    all_genes_in_network = set(G_original.nodes())
    seed_genes = set(seed_genes)
    disease_genes = seed_genes & all_genes_in_network

    if len(disease_genes) != len(seed_genes):
        print("DIAMOnD(): ignoring %s of %s seed genes that are not in the network" % (
            len(seed_genes - all_genes_in_network), len(seed_genes)))


def diable(network_file: Union[str, pd.DataFrame], seed_genes_file: Union[str, pd.DataFrame], num_genes_to_add: int, **kwargs):
    """Runs the DiaBLE algorithm on the provided network and seed genes.

    The DiaBLE algorithm is an iterative variant of DIAMOnD which considers a growing gene universe instead of
    always working with the entire network. At each DiaBLE iteration, the entire gene graph consists out of

    - current seed genes
    - candidate genes (genes with => 1 connection to current seed genes)
    - candidate gene neighbours

    Like in DIAMOnD, we then calculate the p-value of each gene in the set of candidate genes and their neighbours,
    adding the gene with the lowest p-value to the set of seed genes.

    Args:
        network_file (Union[str, pd.DataFrame]): path to the network file or a pandas DataFrame, containing an edge list
         of nodes.
        seed_genes_file (Union[str, pd.DataFrame]): path to the seed genes file or a pandas DataFrame, containing a list
        of seed genes.
        num_genes_to_add (int): the number of genes to add to the network.

    Return:
        (pd.DataFrame): a pandas dataframe with all the DiaBLE genes and their p-values.
    """
    if type(network_file) is not type(seed_genes_file):
        raise TypeError(f"Both the network and seed must be of the same type. Network: {type(network_file)}, Seed: {type(seed_genes_file)}.")

    if isinstance(network_file, pd.DataFrame):
        # ensure that genes are read as strings because when reading from .csv they are read as strings by default
        network_file.iloc[:, 0] = network_file.iloc[:, 0].astype(str)
        network_file.iloc[:, 1] = network_file.iloc[:, 1].astype(str)

        G_original = nx.Graph()
        G_original.add_edges_from(zip(network_file.iloc[:, 0], network_file.iloc[:, 1]))

    if isinstance(seed_genes_file, pd.DataFrame):
        seed_genes = set(map(str, seed_genes_file.iloc[:, 0]))

    else:
        G_original, seed_genes = read_input(network_file, seed_genes_file)

    # 1. throwing away the seed genes that are not in the network
    all_genes_in_network = set(G_original.nodes())
    seed_genes = set(seed_genes)
    disease_genes = seed_genes & all_genes_in_network

    if len(disease_genes) != len(seed_genes):
        print("DIAMOnD(): ignoring %s of %s seed genes that are not in the network" % (
            len(seed_genes - all_genes_in_network), len(seed_genes)))

    neighbours = {node: set(nx.neighbors(G_original, node)) for node in G_original}
    degrees = {node: nx.degree(G_original, node) for node in G_original}

    candidate_genes = set()
    for node in disease_genes:
        candidate_genes |= find_candidate_nodes(node, disease_genes, neighbours)

    added_genes = []
    for _ in tqdm(range(num_genes_to_add), disable=not kwargs.get("verbose", False)):
        new_gene, *_, p_value = diamond_iteration(candidate_genes, disease_genes, neighbours, degrees)

        added_genes.append([new_gene, p_value])

        disease_genes.add(new_gene)
        candidate_genes.remove(new_gene)
        candidate_genes |= find_candidate_nodes(new_gene, disease_genes, neighbours)

    added_genes = pd.DataFrame(added_genes, columns=['gene', "p_value"])
    added_genes.gene = added_genes.gene.astype(str)

    return added_genes


def find_candidate_nodes(gene: str, disease_genes: set, neighbours: dict):
    """Finds the candidate nodes related to the given gene.

    The candidate nodes are the neighbours of the given gene and the neighbours of the neighbours.

    Args:
        gene (str): The gene for which to find the candidate nodes.
        disease_genes (set): The disease genes.
        neighbours (dict): Map for all genes to their neighbours.

    Returns:
        set: The new candidate nodes.
    """
    candidate_nodes = neighbours[gene] - disease_genes
    for node in candidate_nodes.copy():
        candidate_nodes |= neighbours[node]

    candidate_nodes -= disease_genes

    return candidate_nodes


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
    seed_neighbours = find_nodes_with_links_to(network, seed_genes)

    # find the neighbours of the neighbours of the disease genes
    seed_neighbours_neighbours = find_nodes_with_links_to(network, seed_neighbours)

    return network.subgraph(list(G_universe.nodes()) +
                            list(seed_neighbours) +
                            list(seed_neighbours_neighbours))


def diamond_iteration(candidate_genes: Iterable, seed_genes: Sequence[int], neighbours: Dict, degrees: Dict, alpha=1):
    """Runs a single iteration of the DIAMOnD algorithm.

    In the DIAMOnD algorithm, we aim to find genes likely responsible for a disease. A gene that might be connected is
    one that has links to already known disease genes (seed genes). At each iteration in the algorithm, we find the most
    likely responsible gene by running a hypergeometric test for each gene in the universe that does not already exist
    in the seed genes. The new candidate gene that will be added to the seed genes after the iteration is the one with
    the lowest p-value.

    Args:
        universe (nx.Graph): the universe of genes and their interactions to consider. Should contain all the seed genes
        and other genes.
        seed_genes (Sequence[int]): the seed genes responsible for the disease.
        alpha (float): the edge weighting factor.

    Return:
        node, degree, num_links_to_seed_genes, p_value (str, int, int, float): the candidate node with the lowest p-value,
        its degree, number of links to seed genes, and p-value.
    """
    size_universe = len(candidate_genes) + len(seed_genes)

    # ------------------------------------------------------------------
    # Setting up initial set of nodes in cluster
    # ------------------------------------------------------------------
    cluster_nodes = set(seed_genes)
    s0 = len(cluster_nodes)

    s0 += (alpha - 1) * s0
    size_universe += (alpha - 1) * s0

    info = []

    for node in candidate_genes:
        degree = degrees[node]
        num_links_to_seed_genes = len(set(neighbours[node]) & set(seed_genes))

        if num_links_to_seed_genes == 0:
            p = 1.0
        else:
            p = 1 - hypergeom.cdf(num_links_to_seed_genes - 1, size_universe, s0, degree)

        info.append([node, degree, num_links_to_seed_genes, p])

    # find the most likely candidate gene
    candidate_genes = pd.DataFrame(info, columns=['node', 'degree', 'num_links_to_seed_genes', 'p_value'])
    candidate_genes.sort_values("p_value", ascending=True, inplace=True)

    return (*candidate_genes.head(1).iloc[0].tolist(),)
