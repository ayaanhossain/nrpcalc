import hgraph
import recovery
import utils
import kmerSetDB
import projector

import os
import uuid
import atexit
import shutil

import time

def nrp_finder(
    seq_list,
    homology,
    background=None,
    internal_repeats=False,
    vercov_func=None,
    verbose=True):

    # Program Verbage
    if verbose:
        print '\n[Non-Repetitive Parts Calculator - Finder Mode]\n'

    # Background Check
    if not background is None:
        valid_background = True
        if not isinstance(background, kmerSetDB.kmerSetDB):
            print '\n [ERROR]    Background Object is not kmerSetDB'
            print ' [SOLUTION] Please Instantiate Background via nrpcalc.background(...)'
            valid_background = False            
        elif background.K != homology:
            print '\n [ERROR]    Background Lmax is {}, but Part Lmax is {}'.format(
                        background.K,
                        homology)
            print ' [SOLUTION] Please Use Same Lmax Values'
            valid_background = False        

    # Vercov Function Check
    valid_vercov = True
    if not vercov_func in ['2apx', 'nrpG', 'nrp2']:
        print ' [ERROR]    vercov_func must be \'2apx\', \'nrpG\', or \'nrp2\' not \'{}\''.format(
            vercov_func)
        print ' [SOLUTION] Please use One of the Correct Options'
        valid_vercov = False


    if not valid_background or not valid_vercov:
        print '\nPlease Eliminate [ERROR]s Found Above.'
        return {} # Empty Dict ... no parts were subselected

    # Setup Project
    current_uuid = str(uuid.uuid4())
    projector.setup_proj_dir(current_uuid)
    seq_file   = './{}/seq_list.txt'.format(current_uuid)
    graph_file = './{}/repeat_graph.txt'.format(current_uuid)

    # Build Repeat Graph
    homology_graph = hgraph.get_homology_graph(
        list(seq_list),
        background,
        homology,
        internal_repeats,
        seq_file,
        verbose)

    # Find Non-Repetitive Subset
    indi_seqs_indx = recovery.get_recovered_non_homologs(
        homology_graph,
        graph_file,
        vercov_func,
        verbose)

    # Prepare Output
    parts_dict = {}
    with open(seq_file, 'r') as infile:
        for line in infile:
            index,seq = line.strip().split(',')
            index = int(index)
            if index in indi_seqs_indx:
                parts_dict[index] = seq

    # Cleanups and Return
    projector.remove_proj_dir()
    if verbose:
        print '\nNon-Repetitive Toolbox Size: {}'.format(
            len(parts_dict))
    return parts_dict

def main():
    t0 = time.time()

    # Setup Parameters
    homology = 16
    fasta_filename = 'riboz.fa' #'input.fa.bk104'
    seq_list = utils.get_fasta_seq_list(fasta_filename)
    internal_repeats = False

    # Setup Background
    background = kmerSetDB.kmerSetDB(
        path='./testDB',
        homology=homology,
        verbose=True)
    background.multiadd(utils.get_fasta_seq_list(
        fasta_filename='input.fa'))

    # Find Non-Repetititive Subset
    non_homologs = nrp_finder(
        seq_list,
        homology,
        background,
        internal_repeats,
        verbose=True)

    # Drop Background
    background.drop()

    # Report Time Elapsed
    print '\nWall Time {} sec'.format(time.time()-t0)


if __name__ == '__main__':
    main()
