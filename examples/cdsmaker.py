import nrpcalc

iupac_space = {
    'A': {'A'},
    'C': {'C'},
    'G': {'G'},
    'T': {'T'},
    'R': {'A', 'G'},
    'Y': {'C', 'T'},
    'S': {'G', 'C'},
    'W': {'A', 'T'},
    'K': {'G', 'T'},
    'M': {'A', 'C'},
    'B': {'C', 'G', 'T'},
    'V': {'A', 'C', 'G'},
    'D': {'A', 'G', 'T'},
    'H': {'A', 'C', 'T'},
    'N': {'A', 'T', 'G', 'C'}}

inv_iupac_space = {
    tuple(sorted(v)):k for k,v in iupac_space.items()
}

codon_table = {
    'A': {'GCT', 'GCC', 'GCA', 'GCG'},
    'R': {'CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'},
    'N': {'AAT', 'AAC'},
    'D': {'GAT', 'GAC'},
    'C': {'TGT', 'TGC'},
    'Q': {'CAA', 'CAG'},
    'E': {'GAA', 'GAG'},
    'G': {'GGT', 'GGC', 'GGA', 'GGG'},
    'H': {'CAT', 'CAC'},
    'I': {'ATT', 'ATC', 'ATA'},
    'L': {'TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'},
    'K': {'AAA', 'AAG'},
    'M': {'ATG'},
    'F': {'TTT', 'TTC'},
    'P': {'CCT', 'CCC', 'CCA', 'CCG'},
    'S': {'TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'},
    'T': {'ACT', 'ACC', 'ACA', 'ACG'},
    'W': {'TGG'},
    'Y': {'TAT', 'TAC'},
    'V': {'GTT', 'GTC', 'GTA', 'GTG'},
    '*': {'TAA', 'TGA'},
    'O': {'TAG'}}

inv_codon_table = {u:k for k,v in codon_table.items() for u in v}

tracedcodons = set(['R', 'L', 'S'])


def get_metaN(codons):
    '''
    Return meta nucleotide composition without
    caring for codon compatibility.
    '''
    return ''.join(
        inv_iupac_space[tuple(
            sorted(set(s)))] for s in zip(*codons))

def get_AA2DNA(aa_seq):
    '''
    Convert AA Sequence Constraint to Semi-Correct
    Nucleotide Sequence Constraint.
    '''
    return ''.join(
        get_metaN(
            codon_table[aa]) for aa in aa_seq)

def get_DNA2AA(seq):
    '''
    Convert DNA Sequence to AA Sequence.
    '''

    return ''.join(inv_codon_table[
        seq[i:i+3]] for i in range(0, len(seq), 3))

def is_DNA2AA(seq, aa_seq):
    '''
    Local Model Function to preferentially correct
    a partially explored sequence path.
    '''
    if len(seq) % 3 == 0:
        idx = (len(seq) // 3) - 1
        if aa_seq[idx] in tracedcodons:
            if not seq[len(seq)-3:] in codon_table[aa_seq[idx]]:
                assert aa_seq[idx] in tracedcodons
                # print('Violated', seq, seq[len(seq)-3:], aa_seq[idx], codon_table[aa_seq[idx]])
                return (False, len(seq)-1)
    return (True, None)

def augmented_local_model_fn(
    seq,
    aa_seq,
    local_model_fn):
    '''
    Apply local model function specified only if
    the AA Sequence Constraint is satisfied.
    '''

    # Is codon-config correct?
    state,traceloc = is_DNA2AA(
        seq=seq,
        aa_seq=aa_seq)

    # Nope
    if not state:
        return (state, traceloc)

    # How about the remaining?
    else:
        if not local_model_fn is None:
            return local_model_fn(seq)

    # Everything OK!
    return (True, None)


def augmented_global_model_fn(
    seq,
    aa_seq,
    global_model_fn):
    '''
    Apply global model function specified only if
    the AA Sequence Constraint is satisfied.
    '''

    # Is codon-config correct?
    state = get_DNA2AA(seq=seq) == aa_seq

    # Nope
    if not state:
        return False

    # How about the remaining?
    else:
        if not global_model_fn is None:
            return global_model_fn(seq)

    # Everything OK!
    return True

def cdsmaker(
    aa_seq_constr,
    Lmax,
    target_size,
    internal_repeats=False,
    background=None,
    seed=None,
    synth_opt=False,
    local_model_fn=None,
    global_model_fn=None,
    jump_count=10,
    fail_count=1000,
    output_file=None,
    verbose=True):
    '''
    NRP Calculator CDS Maker Mode for building non-repetitive
    coding sequence toolboxes from user defined amino acid
    sequence constraints and based on custom local and/or
    global model functions. All shared repeats longer than
    Lmax are eliminated, and internal repeats longer than
    Lmax are preserved if desired. Parts are optimized for
    DNA synthesis if desired. Error tolerance is adaptive
    and auto-adjusted based on recorded failures. Designed
    toolbox is returned as a dictionary of parts indexed by
    their order of design, and optionally written to a FASTA
    output file.
    '''

    # Basic Book-keeping
    aa_seq     = str(aa_seq_constr)
    seq_constr = get_AA2DNA(
        aa_seq=aa_seq)

    # Augment Local Model Function
    almfn = lambda seq: augmented_local_model_fn(
        seq=seq,
        aa_seq=aa_seq,
        local_model_fn=local_model_fn)

    # Augment Global Model Function
    agmfn = lambda seq: augmented_global_model_fn(
        seq=seq,
        aa_seq=aa_seq,
        global_model_fn=global_model_fn)

    # Compute Parts
    return nrpcalc.maker(
        seq_constr=seq_constr,
        struct_constr='.'*len(seq_constr),
        part_type='DNA',
        Lmax=Lmax,
        target_size=target_size,
        internal_repeats=internal_repeats,
        background=background,
        struct_type='mfe',
        seed=seed,
        synth_opt=synth_opt,
        local_model_fn=almfn,
        global_model_fn=agmfn,
        jump_count=jump_count,
        fail_count=fail_count,
        output_file=output_file,
        verbose=verbose)


def main():

    # Just an example
    def local_prev_cutsites(seq):
        if len(seq) >= 6:
            if seq[-6:] == 'TCTAGA':
                return False, len(seq) - 1

        return True, None

    # Just an example
    def global_prev_cutsites(seq):
        return not 'TCTAGA' in seq

    # Glucagon-like peptide
    aa_seq_constr_1 = 'MGAHGEGTFTSDVSSYLEEQAAKEFIAWLVKGRGAHGEGTFTSDVSSYLEEQAAKEFIAWLVKGRGAHGEGTFTSDVSSYLEEQAAKEFIAWLVKGRGAHGEGTFTSDVSSYLEEQAAKEFIAWLVKGRGAHGEGTFTSDVSSYLEEQAAKEFIAWLVKGRGAHGEGTFTSDVSSYLEEQAAKEFIAWLVKGR*'

    # Designed Ankyrin Repeat Protein
    aa_seq_constr_2 = 'MDLGKKLLEAARAGQDDEVRILMANGADVNADDTWGWTPLHLAAYQGHLEIVEVLLKNGADVNYDYIGWTPLHLAADGHLEIVEVLLKNGADVNASDYIGDTPLHLAAHNGHLEIVEVLLKHGADVNAQDKFGKTAFDISIDNGNEDLAEILQ*'

    # GVGVP Repeats
    aa_seq_constr_3 = 'GVGVP'*100

    for aa_seq_constr in (aa_seq_constr_3,):#, aa_seq_constr_2, aa_seq_constr_3):
        for seq in cdsmaker(
            aa_seq_constr=aa_seq_constr,
            Lmax=19,
            target_size=30,
            internal_repeats=False,
            background=None,
            seed=None,
            synth_opt=False,
            local_model_fn=local_prev_cutsites,
            global_model_fn=global_prev_cutsites,
            jump_count=100,
            fail_count=1000,
            output_file=None,
            verbose=True).values():

            assert get_DNA2AA(seq=seq) == aa_seq_constr
            assert not 'TCTAGA' in seq

            print(seq)


if __name__ == '__main__':
    main()