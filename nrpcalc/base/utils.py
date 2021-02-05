import os
import sys
from typing import Tuple

import pkg_resources

from Bio import SeqIO

import RNA

import numpy as np

complement_table = str.maketrans('ATGCU', 'TACGA')

def stream_fasta_seq_list(fasta_filename):
    with open(fasta_filename, "rU") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            yield str(record.seq)

def get_fasta_seq_list(fasta_filename):
    return list(stream_fasta_seq_list(fasta_filename))

def stream_txt_seq_list(text_filename):
    with open(text_filename) as infile:
        for line in infile:
            yield line.strip()

def get_txt_seq_list(text_filename):
    return list(stream_txt_seq_list(text_filename))

def uniquify_background_list(background_list):
    uniq_background_set = set()
    while background_list:
        uniq_background_set.add(background_list.pop())
    background_list = []
    while uniq_background_set:
        background_list.append(uniq_background_set.pop())
    return background_list

def stream_kmers(seq, k):
    if k >= len(seq):
        return [seq]
    return (seq[i:i+k] for i in range(len(seq)-k+1))

def get_comp(seq):
    return seq.translate(complement_table)

def get_revcomp(seq):
    return get_comp(seq)[::-1]

def stream_min_kmers(seq, k):
    for kmer in stream_kmers(seq, k):
        yield min(kmer, get_revcomp(kmer))

class Fold(object):

    def __init__(
        self,
        temp=37.0,
        dangles=2,
        part_type='RNA'):

        if not part_type in ['RNA', 'DNA']:
            part_type = 'RNA'

        if part_type == 'DNA':
            RNA.cvar.noGU        = True
            RNA.cvar.noGUclosure = True

        self.parameter_directory = os.path.dirname(
            os.path.abspath(__file__))#"/usr/local/share/ViennaRNA/"
        # Temperature in Celsius;
        # default=37.0 (float)
        RNA.cvar.temperature = temp
        # Dangling end energies (0,1,2);
        # see RNAlib documentation;
        # default=2 (int)
        RNA.cvar.dangles = dangles
        self.settings    = RNA.md()

        self.part_type = part_type

        parameter_file = pkg_resources.resource_filename(
            'nrpcalc', 'base/{}.par'.format(
                self.part_type))
        RNA.read_parameter_file(parameter_file)

        if part_type == 'DNA':
            self.clear_warning()

        self.adjust = self.adjust_dG(temp)

    def adjust_dG(self, temp):
        # Adjustment according to Dirks et al.

        kB = 0.00198717 # Boltzmann constant in kcal/mol/K
        T = temp
        a = [-3.983035, 301.797, 522528.9, 69.34881, 999.974950]

        # Calculate the number of moles of water per liter (molarity) at temperature (T in deg C)
        # Density of water calculated using data from
        # Tanaka M., Girard, G., Davis, R., Peuto A., Bignell, N.
        # Recommended table for the density of water..., Metrologia, 2001, 38, 301-309
        pH2O = a[4] * (
            1 - (T+a[0])**2 * (T+a[1]) / \
                (a[2]) /                 \
                (T+a[3])) /              \
            18.0152

        return -kB * (T + 273.15) * np.log(pH2O)

    def clear_warning(self):
        clrlen = len('WARNING: stacking enthalpies not symmetric')
        sys.stdout.write('\033[F\033[F\033[F\033[F')
        sys.stdout.write(' '*clrlen+'\n')
        sys.stdout.write(' '*clrlen+'\n')
        sys.stdout.write(' '*clrlen+'\n')
        sys.stdout.write(' '*clrlen+'\n')
        sys.stdout.write('\033[F\033[F\033[F\033[F')
        sys.stdout.flush()

    def evaluate_mfe(self, seq, dg=False):
        # MFE Structure Only
        fc_obj = RNA.fold_compound(seq, self.settings)
        struct,energy = fc_obj.mfe()
        if not dg:
            return struct
        else:
            return struct, energy

    def evaluate_centroid(self, seq, dg=False):
        # Centroid Structure Only
        fc_obj = RNA.fold_compound(seq, self.settings)
        fc_obj.pf()
        struct,energy = fc_obj.centroid()
        if not dg:
            return struct
        else:
            return struct, energy

    def design(self, seq, struct):
        # Closest MFE Structure Sequence
        inv = RNA.inverse_fold(seq, struct)[0]
        if self.part_type == 'DNA':
            inv = inv.replace('U', 'T').replace('u', 't')
        return inv

    def evaluate_mfe_dimer(self, seq1, seq2):
        # MFE Dimer Structure and Energy
        fc_obj = RNA.fold_compound(seq1+'&'+seq2, self.settings)
        struct,energy = fc_obj.mfe_dimer()
        struct1 = struct[:len(seq1)]
        struct2 = struct[len(seq1):]
        energy += self.adjust
        return (struct1, struct2, energy)

if __name__ == '__main__':
    pass