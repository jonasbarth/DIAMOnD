import networkx as nx

from diamond.diable import find_nodes_with_links_to, create_diable_universe


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
