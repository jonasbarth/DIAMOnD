import os
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd

from diamond.diable import find_nodes_with_links_to, create_diable_universe, diable


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


def test_that_diable_p_value_accurate():
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
