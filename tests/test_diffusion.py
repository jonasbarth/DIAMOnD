"""Diffusion tests."""
from diamond import diffusion


def test_diffusion():
    network_file = "slide_graph_universe.txt"
    seed_file = "slide_graph_seed_genes.txt"

    diffusion_result = diffusion(network_file, seed_file, 2)

    assert len(diffusion_result) == 2
