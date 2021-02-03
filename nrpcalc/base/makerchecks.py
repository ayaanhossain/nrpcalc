import random
import functools
import collections

iupac_count = {
    'A': 1, 'C': 1, 'G': 1, 'T': 1,
    'R': 2, 'Y': 2, 'S': 2, 'W': 2,
    'K': 2, 'M': 2, 'B': 3, 'D': 3,
    'H': 3, 'V': 3, 'N': 4, 'U': 1
}

iupac_space_rna = {
    'A': {'A'}, 'C': {'C'}, 'G': {'G'}, 'U': {'U'},
    'R': {'A', 'G'}, 'Y': {'C', 'U'},
    'S': {'G', 'C'}, 'W': {'A', 'U'},
    'K': {'G', 'U'}, 'M': {'A', 'C'},
    'B': {'C', 'G', 'U'}, 'V': {'A', 'C', 'G'},
    'D': {'A', 'G', 'U'}, 'H': {'A', 'C', 'U'},
    'N': {'A', 'U', 'G', 'C'}
}

iupac_compl_rna = {
    'A': 'U', 'C': 'G', 'G': 'Y', 'U': 'R',
    'R': 'Y', 'Y': 'R', 'S': 'B', 'W': 'D',
    'K': 'N', 'M': 'K', 'B': 'N', 'V': 'D',
    'D': 'N', 'H': 'N', 'N': 'N'
}

iupac_space_dna = {
    'A': {'A'}, 'C': {'C'}, 'G': {'G'}, 'T': {'T'},
    'R': {'A', 'G'}, 'Y': {'C', 'T'},
    'S': {'G', 'C'}, 'W': {'A', 'T'},
    'K': {'G', 'T'}, 'M': {'A', 'C'},
    'B': {'C', 'G', 'T'}, 'V': {'A', 'C', 'G'},
    'D': {'A', 'G', 'T'}, 'H': {'A', 'C', 'T'},
    'N': {'A', 'T', 'G', 'C'}
}

iupac_compl_dna = {
    'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
    'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W',
    'K': 'M', 'M': 'K', 'B': 'V', 'V': 'B',
    'D': 'H', 'H': 'D', 'N': 'N'
}

def is_homolog_legal(seq, homology):
    '''
    Check is homology is satisfiable
    '''
    if homology > len(seq):
        return False
    return True

def is_seq_constr_legal(seq, part_type):
    '''
    Check if all characters used is standard IUPAC code.
    On failure returns (False, a list of illegal chars used).
    '''
    global iupac_space_dna
    global iupac_space_rna
    if part_type == 'DNA':
        iupac_space = iupac_space_dna
    else:
        iupac_space = iupac_space_rna
    chars    = set(seq)
    alphabet = set(iupac_space.keys())
    if chars <= alphabet:
        return (True, None)
    else:
        return (False, sorted(chars-alphabet))

def get_computable_form(struct):
    '''
    Get computable components from an RNA secondary structure.
    '''
    opened  = []
    closed  = []
    pairs   = []
    invalid = []
    for pos in range(len(struct)):
        pairing = struct[pos]
        if pairing == '(':
            opened.append(pos)
        elif pairing == ')':
            closed.append(pos)
            try:
                pairs.append((opened.pop(), closed.pop()))
            except:
                pass
        elif pairing not in ['x', '.']:
            invalid.append(pos)
    return pairs, opened, closed, invalid

def get_k_mer_count(seq_list, start, end, rev_dict):
    '''
    Get the possible number of candidate sequences.
    '''
    global iupac_count
    product = 1
    for j in range(start, end):
        nt = seq_list[j]
        if j in rev_dict:
            i = rev_dict[j]
            if i >= start:
                continue
        product *= iupac_count[nt]
    return product

def compress_locs(locs, homology):
    '''
    Compress contigs in locs.
    '''
    if not locs:
        return None
    else:
        compressed = []
        x = locs[0]
        y = x+homology
        for pos in range(1, len(locs)):
            if locs[pos]-locs[pos-1] == 1:
                y = locs[pos]+homology
            else:
                compressed.append((x, y))
                x = locs[pos]
                y = x+homology
        compressed.append((x, y))
    return compressed

def is_seq_constr_sufficient(seq, struct, homology, toolbox_size):
    '''
    Check if a desired size toolbox can be generated.
    On failure returns (False, a list of start:end tuples with constricted motifs).
    '''
    seq_list = list(seq)
    sufficiency_status = True
    insufficiency_locs = []
    rev_dict = {j:i for i,j in get_computable_form(struct)[0]}
    for i in range(len(seq)-homology+1):
        start = i
        end = start + homology
        kmer_count = get_k_mer_count(
            seq_list=seq_list,
            start=start,
            end=end,
            rev_dict=rev_dict)
        if kmer_count < toolbox_size:
            sufficiency_status = False
            insufficiency_locs.append(i)
    return (sufficiency_status, compress_locs(locs=insufficiency_locs, homology=homology))

def is_structure_valid(struct):
    '''
    Check if the given RNA structure is legal and balanced.
    On failure returns (False, three lists of indices:
    unclosed, unopened and invalid charas/parens)
    '''
    pairs, opened, closed, invalid = get_computable_form(struct)
    if not opened:
        opened  = None
    if not closed:
        closed  = None
    if not invalid:
        invalid = None
    if opened or closed or invalid:
        return (False, opened, closed, invalid)
    else:
        return (True, None, None, None)

def hamming(a, b):
    '''
    Get hamming distance between two strings.
    '''
    hdist = 0
    for i in range(len(a)):
        if a[i] != b[i]:
            hdist += 1
    return hdist

def is_contiguous(item1, item2):
    '''
    Check if two pairs are contiguous.
    '''
    if (item1[0]-item2[0] <= 1) and (item2[1]-item1[1] <= 1):
        return True
    return False

def is_structure_not_conflict(struct, homology):
    '''
    Q. Does the structure not enforce internal repeats?
    A. Stack evalation of stack of pairs.
    On failure returns (False, a list of locations with stems >= homology)
    '''
    pairs, opened, closed, invalid = get_computable_form(struct)
    if pairs:
        bad_contigs = []
        stretch     = 1
        pairs       = pairs[::-1]
        start       = pairs.pop()
        end         = None
        item1       = start
        item2       = end
        while pairs:
            end   = pairs.pop()
            item2 = end
            if is_contiguous(item1, item2):
                item1 = item2
                item2 = None
                stretch += 1
            else:
                if stretch >= homology:
                    bad_contigs.append(((item1[0], start[0]), (start[1], item1[1])))
                stretch = 1
                start = item2
                item1 = start
        if stretch >= homology:
            bad_contigs.append(((end[0], start[0]), (start[1], end[1])))
        if bad_contigs:
            return (False, bad_contigs)
    return (True, None)

def is_pairing_compatible(seq, struct, part_type):
    global iupac_space_dna
    global iupac_compl_dna
    global iupac_space_rna
    global iupac_compl_rna
    pairs, opened, closed, invalid = get_computable_form(struct)
    incompat_locs = []
    reduced_locs = []
    if pairs:
        if part_type == 'DNA':
            iupac_space = iupac_space_dna
            iupac_compl = iupac_compl_dna
        else:
            iupac_space = iupac_space_rna
            iupac_compl = iupac_compl_rna
        while pairs:
            i,j = pairs.pop()
            nti = seq[i]
            ntj = seq[j]
            if len(iupac_space[nti]) > len(iupac_space[ntj]):
                nti,ntj = ntj,nti
            if not iupac_space[nti] <= iupac_space[iupac_compl[ntj]]:
                incompat_locs.append(tuple(sorted([i,j])))
            else:
                if iupac_space[iupac_compl[ntj]] > iupac_space[nti]:
                    reduced_locs.append(tuple(sorted([i,j])))
        if len(incompat_locs) > 0 or \
           len(reduced_locs) > 0:
            return (False, incompat_locs, reduced_locs)
    return (True, incompat_locs, reduced_locs)

def main():
    print('Testing Sequence Constraint Validation Checkers ... ',)
    legal_alphabet  = set('ACGTRYSWKMBDHVN')
    seq_legal = ''.join(random.sample(legal_alphabet, 1)[0] for _ in range(100))
    assert is_seq_constr_legal(seq=seq_legal, part_type='DNA')[0] == True

    illegal_alphabet = set('ABCDEFGHIJKLMNOPQRSTUVWXYZ')
    seq_illegal = ''.join(random.sample(illegal_alphabet, 1)[0] for _ in range(100))
    assert is_seq_constr_legal(seq=seq_illegal, part_type='DNA')[0] == False
    print('OK')

    print('Testing Sequence Constraint Sufficiency Checkers ... ',)
    seq      = 'N'*20 + 'TTGACA' + 'N'*17 + 'TATAAT' + 'N'*6 + 'CCN' + 'N'*20
    struct   = '.'*len(seq)
    homology = 11
    assert is_seq_constr_sufficient(seq, struct, homology, toolbox_size=1000)[0] == True
    assert is_seq_constr_sufficient(seq, struct, homology, toolbox_size=2000)[0] == False
    struct = '......((((((((((xxxxxxxx))))))))))......'
    seq    = 'NNNNNNNNNNNNNNNNAAAAAAAANNNNNNNNNNNNNNNN'
    assert is_seq_constr_sufficient(seq, struct, homology, toolbox_size=1000)[0] == False
    assert is_seq_constr_sufficient(seq, struct, 17, toolbox_size=1000)[0] == True
    print('OK')

    print('Testing Struture Constraint Validation Checkers ...',)
    struct = '....(((.xx.)))...()...(())...'
    assert is_structure_valid(struct)[0] == True
    struct = '....(((.xx.))).).()...(())...'
    assert is_structure_valid(struct)[0] == False
    struct = '.(..(((.xx.)))...()...(())...'
    assert is_structure_valid(struct)[0] == False
    struct = '....(((.xx.))).X.()...(())...'
    assert is_structure_valid(struct)[0] == False
    struct = '.B..(((.xx.))).).().(.(())...'
    assert is_structure_valid(struct)[0] == False
    print('OK')

    print('Testing Struture Constraint Homology Conflict Checkers ...',)
    struct = '(((((((((((......)))))))))))...(((((((((((......)))))))))))'
    homology = 10
    assert is_structure_not_conflict(struct, homology)[0] == False
    struct = '(((((((((((......))))).))))))...(((((((((((......)))))))))))'
    homology = 10
    assert is_structure_not_conflict(struct, homology)[0] == False
    struct = '(((((((((((......)))).)))))))...(((((.((((((......)))))))))))'
    homology = 10
    assert is_structure_not_conflict(struct, homology)[0] == True
    print('OK')

    print('Testing Sequence Constraint Base Pairing Conflict Checkers ...',)
    seq    = 'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN'
    struct = '.....(((..(((((....)))))..))).....'
    assert is_pairing_compatible(seq, struct, part_type='DNA') == (True, [], [])
    seq    = 'NNNNNNNNNNSNNNNNNNNNNNNSNNNNNNNNNN'
    struct = '.....(((..(((((....)))))..))).....'
    assert is_pairing_compatible(seq, struct, part_type='DNA') == (True, [], [])
    seq    = 'NNNNNNNNNNSNNNNNNNNNNNNWNNNNNNNNNN'
    struct = '.....(((..(((((....)))))..))).....'
    assert is_pairing_compatible(seq, struct, part_type='DNA') == (False, [(10, 23)], [])
    seq    = 'NNNNNNNNNNNNNNNNNNNNNNNANNNNNNNNNN'
    struct = '.....(((..(((((....)))))..))).....'
    assert is_pairing_compatible(seq, struct, part_type='RNA') == (False, [], [(10, 23)])
    seq    = 'NNNNNANNNNNNNNNNNNNNNNNNNNNNNNNNNN'
    struct = '.....(((..(((((....)))))..))).....'
    assert is_pairing_compatible(seq, struct, part_type='RNA') == (False, [], [(5, 28)])
    seq    = 'NNNNNMNNNNNNNNNNNNNNNNNNNNNNVNNNNN'
    struct = '.....(((..(((((....)))))..))).....'
    assert is_pairing_compatible(seq, struct, part_type='RNA') == (False, [(5, 28)], [])
    seq    = 'NNNNNHNNNNNNNNNNNNNNNNNNNNNNCNNNNN'
    struct = '.....(((..(((((....)))))..))).....'
    assert is_pairing_compatible(seq, struct, part_type='RNA') == (False, [], [(5, 28)])
    seq    = 'NNNNNGNNNNNNNNNNNNNNNNNNNNNNTNNNNN'
    struct = '.....(((..(((((....)))))..))).....'
    assert is_pairing_compatible(seq, struct, part_type='DNA') == (False, [(5, 28)], [])
    seq    = 'NNNNNGNNNNNNNNNNNNNNNNNNNNNNUNNNNN'
    struct = '.....(((..(((((....)))))..))).....'
    assert is_pairing_compatible(seq, struct, part_type='RNA') == (False, [], [(5, 28)])
    print('OK')

if __name__ == '__main__':
    main()

