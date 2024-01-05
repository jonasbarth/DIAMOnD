import os
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd
import pytest

from diamond.diable import find_nodes_with_links_to, create_diable_universe, diable, diamond_iteration, find_candidate_nodes


def test_that_find_nodes_with_links_to_single_node():
    graph = nx.Graph([(1, 2), (2, 3), (1, 3), (3, 4)])
    nodes = set([1])

    links = find_nodes_with_links_to(graph, nodes)

    assert len(links) == 2
    assert links == {2, 3}


def test_that_find_nodes_with_links_to_multiple_nodes():
    graph = nx.Graph([(1, 2), (2, 3), (1, 3), (3, 4)])
    nodes = set([1, 2])

    links = find_nodes_with_links_to(graph, nodes)

    assert len(links) == 1
    assert links == {3}


def test_that_create_diable_universe():
    graph = nx.Graph([(1, 2), (2, 3), (3, 4), (4, 5)])
    nodes = set([1])

    universe = create_diable_universe(graph, nodes)

    assert len(universe) == 3
    assert list(universe.nodes()) == [1, 2, 3]
    assert list(universe.edges()) == [(1, 2), (2, 3)]


def test_that_diable_output_correct():
    edge_list = [[1, 2], [1, 5], [2, 3], [2, 4], [2, 5], [3, 4], [4, 5]]
    tmp_dir = "tmp"
    Path(tmp_dir).mkdir(parents=True, exist_ok=True)
    network_file = os.path.join(tmp_dir, "test_edge_list.csv")
    seed_file = os.path.join(tmp_dir, "test_seed_genes.csv")
    pd.DataFrame(edge_list).to_csv(network_file, header=False, index=False)
    pd.DataFrame([1]).to_csv(seed_file, header=False, index=False)

    diable_result = diable(network_file, seed_file, 1)
    candidate_gene, p_value = (*diable_result.iloc[0],)

    assert len(diable_result) == 1
    assert candidate_gene == "5"
    assert np.isclose(p_value, 0.6)


def test_that_diable_raises_error_for_different_types():
    with pytest.raises(TypeError):
        diable("fake_path", pd.DataFrame(), 1)


def test_that_diable_output_correct_dataframe():
    edge_list = [[1, 2], [1, 5], [2, 3], [2, 4], [2, 5], [3, 4], [4, 5]]
    network_file = pd.DataFrame(edge_list)
    seed_file = pd.DataFrame([1])

    diable_result = diable(network_file, seed_file, 1)
    candidate_gene, p_value = (*diable_result.iloc[0],)

    assert len(diable_result) == 1
    assert candidate_gene == "5"
    assert np.isclose(p_value, 0.6)


def test_diable_with_slide_example():
    network_file = "slide_graph_universe.txt"
    seed_file = "slide_graph_seed_genes.txt"

    diable_result = diable(network_file, seed_file, 1)
    candidate_gene, p_value = (*diable_result.iloc[0],)

    assert len(diable_result) == 1
    assert candidate_gene == "7"
    assert np.isclose(p_value, 0.047, atol=0.001)


def test_diable_with_more_iterations():
    network_file = "slide_graph_universe.txt"
    seed_file = "slide_graph_seed_genes.txt"

    diable_result = diable(network_file, seed_file, 2)

    assert len(diable_result) == 2
    assert diable_result.gene.tolist() == ["7", "8"]
    assert np.isclose(diable_result.p_value.tolist(), [0.047, 0.20], atol=0.005).all()


def test_diamond_interation():
    network_file = "slide_graph_universe.txt"
    network_df = pd.read_csv(network_file, header=None)
    network = nx.Graph()
    network.add_edges_from(zip(network_df.iloc[:, 0], network_df.iloc[:, 1]))
    seed_file = "slide_graph_seed_genes.txt"
    seed_genes = pd.read_csv(seed_file, header=None).iloc[:, 0]
    candidate_genes = network.nodes() - seed_genes
    neighbours = {node: set(nx.neighbors(network, node)) for node in network}
    degrees = {node: nx.degree(network, node) for node in network}

    gene, degree, num_links_to_seed_genes, p_value = diamond_iteration(candidate_genes, seed_genes, neighbours, degrees)

    assert str(int(gene)) == "8"
    assert degree == 12
    assert num_links_to_seed_genes == 6
    assert np.isclose(p_value, 0.0025, atol=0.001)


def test_two_diamond_interation():
    network_file = "slide_graph_universe.txt"
    network_df = pd.read_csv(network_file, header=None)
    network = nx.Graph()
    network.add_edges_from(zip(network_df.iloc[:, 0], network_df.iloc[:, 1]))
    seed_file = "slide_graph_seed_genes.txt"
    seed_genes = pd.read_csv(seed_file, header=None).iloc[:, 0].tolist()
    candidate_genes = network.nodes() - seed_genes
    neighbours = {node: set(nx.neighbors(network, node)) for node in network}
    degrees = {node: nx.degree(network, node) for node in network}

    gene, *_ = diamond_iteration(candidate_genes, seed_genes, neighbours, degrees)
    seed_genes.append(str(int(gene)))
    candidate_genes |= neighbours[gene]
    for node in candidate_genes.copy():
        candidate_genes |= neighbours[node]
    candidate_genes -= set(seed_genes)
    gene, degree, num_links_to_seed_genes, p_value = diamond_iteration(candidate_genes, seed_genes, neighbours, degrees)

    assert str(int(gene)) == "7"
    assert degree == 5
    assert num_links_to_seed_genes == 4
    assert np.isclose(p_value, 0.0076, atol=0.001)


def test_that_find_candidate_genes():
    network_file = "slide_graph_universe.txt"
    network_df = pd.read_csv(network_file, header=None, dtype=str)
    network = nx.Graph()
    network.add_edges_from(zip(network_df.iloc[:, 0], network_df.iloc[:, 1]))
    seed_file = "slide_graph_seed_genes.txt"
    seed_genes = set(pd.read_csv(seed_file, header=None, dtype=str).iloc[:, 0].tolist())
    seed_genes.add("7")
    neighbours = {node: set(nx.neighbors(network, node)) for node in network}

    candidate_nodes = find_candidate_nodes("7", seed_genes, neighbours)

    assert candidate_nodes == {"9", "17", "18", "16"}