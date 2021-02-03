import sys

from . import vercov

import networkx as nx

from itertools import count
from time      import time
from math      import log10

def is_graph_empty(homology_graph):
    if not homology_graph.number_of_nodes():
        return True
    return False

def get_vercov_func(vercov_func, homology_graph):
    # TO DO: NEED TO USE DIFFERENT FUNCTIONS AT DIFFERENT EDGE COUNT SCALES!!
    if vercov_func == 'nrpG':
        return vercov.nrp_vercov_greedy, 'NRP Greedy'
    elif vercov_func == '2apx':
        # return vercov.std_vercov_approx, 'Standard 2-approximation'
        return vercov.nx_vercov, 'Standard 2-approximation'
    else:
        return vercov.nrp_vercov_approx, 'NRP 2-approximation'

def dump_homology_graph(homology_graph, graph_file):
    nx.write_adjlist(homology_graph, graph_file)

def load_homology_graph(graph_file):
    return nx.read_adjlist(graph_file, nodetype=int)

def get_recovered_non_homologs(homology_graph, graph_file, vercov_func=None, verbose=True):
    indiset_nodes   = set()
    # powertex_elimination(homology_graph, verbose)
    completex_nodes = [] #completex_elimination(homology_graph, verbose)
    # indiset_nodes.update(completex_nodes)

    possible_nodes  = set(homology_graph.nodes())

    if verbose:
        print('\n [+] Initial independent set = {}, computing vertex cover on remaining {} nodes.'.format(len(indiset_nodes), len(completex_nodes), len(possible_nodes)))

    if is_graph_empty(homology_graph):
        if verbose:
            print(' [X] Graph is empty, further independent set expansion not possible, terminating.')
    else:
        vercov_func, vercov_func_name = get_vercov_func(vercov_func, homology_graph)
        iteration = -1

        if verbose:
            print(' [+] Vertex Cover Function: {}'.format(vercov_func_name))
            sys.stdout.write(' [+] Dumping graph into: {}'.format(graph_file))

        t0 = time()
        dump_homology_graph(homology_graph, graph_file)
        if verbose:
            sys.stdout.write(' in {} seconds\n'.format(time()-t0))

        while True:
            iteration += 1

            if verbose:
                print('\n----------------------')
                print('Now running iteration: {}'.format(iteration))
                print('----------------------')

            t0 = time()

            if iteration > 0:
                homology_graph = nx.Graph(homology_graph.subgraph(possible_nodes))

            vercov_nodes = vercov_func(homology_graph, verbose)

            if verbose:
                print('\n [+] Computed vertex cover of size: {} (in {:.4} seconds)'.format(len(vercov_nodes), time()-t0))
                print(' [+] Loading graph from: {}'.format(graph_file))

            homology_graph = load_homology_graph(graph_file)

            new_indiset_nodes = possible_nodes - vercov_nodes
            indiset_nodes |= new_indiset_nodes

            possible_nodes = vercov_nodes
            prev_possibility_count = len(possible_nodes)
            for indi_node in new_indiset_nodes:
                possible_nodes.difference_update(homology_graph[indi_node])
            curr_possibility_count = len(possible_nodes)

            if verbose:
                print(' [+] Current independent set size:  {}'.format(len(indiset_nodes)))
                print(' [+] Potential nodes for expansion: {} (projected independent set size: {})'.format(len(possible_nodes), len(indiset_nodes)+len(possible_nodes)))
            if len(possible_nodes) == 0 or prev_possibility_count == curr_possibility_count:
                if verbose:
                    print(' [X] Cannot expand independent set, terminating.')
                break

    return indiset_nodes