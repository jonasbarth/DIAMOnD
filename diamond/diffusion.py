"""Contains code for cytoscapes diffusion algorithm."""
from typing import Union

import ndex2
import networkx as nx
import pandas as pd
import requests

from .diable import read_input


def diffusion(network_file: Union[str, pd.DataFrame], seed_genes_file: Union[str, pd.DataFrame], num_genes_to_add: int,
              **kwargs):
    """Runs the cytoscape diffusion algorithm on the provided network and seed genes.

    The general idea behind diffusion analysis is to model the flow of information or influence from a set of source
    nodes to other nodes in the network. Each node is assigned a "heat" value, describing the amount of information that
    this node is able to spread throughout the network. A higher value means it can spread more information.

    Args:
        network_file (Union[str, pd.DataFrame]): path to the network file or a pandas DataFrame, containing an edge list
         of nodes.
        seed_genes_file (Union[str, pd.DataFrame]): path to the seed genes file or a pandas DataFrame, containing a list
        of seed genes.
        num_genes_to_add (int): the number of genes to add to the network.
        **kwargs: additional arguments
            time (float): the upper bound on the exponential multiplication performed by diffusion. Default value is 0.001.

    Return:
        (pd.DataFrame): a pandas dataframe with all the diffusion genes and their heat values.
    """
    if type(network_file) is not type(seed_genes_file):
        raise TypeError(f"Both the network and seed must be of the same type. Network: {type(network_file)}, "
                        f"Seed: {type(seed_genes_file)}.")

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

    cx_network = ndex2.create_nice_cx_from_networkx(G_original)

    for node_id, node in cx_network.get_nodes():
        if node["n"] in seed_genes:
            cx_network.set_node_attribute(type='double', node=node_id, attribute_name='diffusion_input', values='1.0')

    diffusion_response = _call_diffusion_service(cx_network, kwargs.get("time", 0.01))
    if diffusion_response.ok:
        diffusion_network = _add_diffusion_result_to_network(cx_network, diffusion_response)
        diffusion_df = _create_diffusion_result_df(diffusion_network)
        diffusion_df.sort_values("rank", inplace=True)
        diffusion_df.gene = diffusion_df.gene.astype(str)

        return diffusion_df.head(num_genes_to_add)[["gene", "heat"]]
    raise Exception(f"Error when running diffusion algorithm: {diffusion_response.content}")

def _call_diffusion_service(cx_network: ndex2.NiceCXNetwork, time: float):
    """Calls the cytoscape diffusion service.

    The cytoscape diffusion service is a webservice offered by cytoscape that runs the diffusion algorithm.
    """
    url = f'http://v3.heat-diffusion.cytoscape.io?time={time}'

    payload = cx_network.to_cx()
    for p in payload:
        k = list(p.keys())[0]
        if 'Attributes' in k:
            for i in range(len(p[k])):
                p[k][i]['v'] = str(p[k][i]['v'])

    response = requests.post(url, json=payload)
    return response


def _add_diffusion_result_to_network(network: ndex2.NiceCXNetwork, response: requests.Response):
    """Adds the diffusion HTTP response to the network.

    The diffusion HTTP response contains ranks and heat values for nodes. These will be assigned to an existing CX
    network.
    """
    for aspect in response.json()['data']:
        if 'nodeAttributes' not in aspect:
            continue
        for n_attr in aspect['nodeAttributes']:
            if n_attr['n'] == 'diffusion_output_rank' or n_attr['n'] == 'diffusion_output_heat':
                if n_attr['d'] == 'float':
                    n_type = 'double'
                else:
                    n_type = n_attr['d']
                network.add_node_attribute(property_of=int(n_attr['po']),
                                           name=n_attr['n'],
                                           values=n_attr['v'],
                                           type=n_type)

    return network


def _create_diffusion_result_df(network: ndex2.NiceCXNetwork):
    """Creates pandas dataframe of the diffusion result network.

    The dataframe will contain:
    - rank: rank of the node
    - gene: the id of the gene
    - heat: heat of the node
    """
    diffusion_nodes = []

    for node_id, node in network.get_nodes():
        rank = network.get_node_attribute_value(node_id, 'diffusion_output_rank')
        heat = network.get_node_attribute_value(node_id, 'diffusion_output_heat')

        diffusion_nodes.append([int(rank), node["n"], float(heat)])

    diffusion_pred = pd.DataFrame(diffusion_nodes, columns=["rank", "gene", "heat"])

    return diffusion_pred
