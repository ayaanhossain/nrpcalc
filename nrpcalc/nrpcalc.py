from .base import maker     as nrpmaker
from .base import finder    as nrpfinder
from .base import kmerSetDB


__version__ = '1.6.3'

__authors__ = '''
Ayaan Hossain <auh57@psu.edu>
Howard Salis  <salis@psu.edu>
'''

__doc__ = '''
Non-Repetitive Parts Calculator

Automated design and discovery of non-repetitive genetic
parts for engineering stable systems.

Version: {}

Authors: {}
         {}

The Non-Repetitive Parts Calculator offers two modes of operation:

- Finder Mode: Discover toolboxes of non-repetitive parts
               from a list of candidate parts

-  Maker Mode: Design toolboxes of non-repetitive parts
               based on sequence, structure and model
               constraints

Additionally, a 'background' object is available that stores
background sequences. When the 'background' object is used, 
designed genetic parts will also be non-repetitive with respect
to these sequences.

You can learn more about the two modes and background via
  print(nrpcalc.background.__doc__)
  print(nrpcalc.finder.__doc__)
  print(nrpcalc.maker.__doc__)
'''.format(
    __version__,
    *__authors__.strip().split('\n'))


def background(
    path,
    Lmax,
    verbose=True):
    '''
    NRP Calculator kmerSetDB background object for on-disk
    storage of background sequence k-mers (where k=Lmax+1).
    When a sequence is added to background, k-mers from the
    sequence are added instead of the actual sequence itself.
    A sequence queried for existence in the given background
    is evaluated to be True if any k-mer from the sequence
    exists in the background. This object is useful when
    chaining multiple Maker Mode and Finder Mode jobs as
    well as persisting any backgrounds such as small genomes
    or other part toolboxes from earlier design rounds.

    :: path
       type - string
       desc - ./a/path/to/store/background/object for later
              reuse and part verification
    :: Lmax
       type - integer
       desc - maximum allowed shared repeat length between
              all sequences in a given toolbox
    :: verbose
       type - boolean
       desc - if True displays progress
              (default=True)

    Returns: A kmerSetDB object.

    Note: If the path provided points to an existing background
          object, then that background is opened for reading, and
          and the new Lmax is ignored, otherwise, a new background
          is instantiated at the given path.

    background / kmerSetDB API Examples
    
    >>> from pprint import pprint
    >>>
    >>> import nrpcalc
    >>>
    >>> my_background_list = [
        'ATGAGATCGTAGCAACC',
        'GACGATTACGTCAGGTA',
        'ACAGTAGAGACGAGTAA',
        'CCAGTACGAAAAGGCCC',
        'TTAGCTTGATAGTTTTA']
    >>> bkg = nrpcalc.background(
            path='./prj_bkg/',
            Lmax=15)
    >>> bkg
    kmerSetDB stored at ./prj_bkg/ with 0 16-mers
    >>>    

    background / kmerSetDB offers the following methods

    (1) add(seq) - adds an IUPAC string 'seq' to background

    >>> bkg.add('ATGCTTAGTGCCATACC')
    
    (2) multiadd(seq_list) - adds multiple sequences from
                             the list to background

    >>> bkg.multiadd(my_background_list)
    
    [Background Processing]
      Adding Seq 4: TTAGCTTGAT...
    
    (3) __contains__(seq) - checks if all k-mers from seq
                            is present in background

    >>> 'ATGCTTAGTGCCATACC' in bkg
    True

    (4) multicheck(seq_list) - checks if all k-mers from given
                               seq_list present in background

    >>> assert all(bkg.multicheck(my_background_list))

    (5) __iter__() - iterate over all k-mers in background

    >>> pprint(list(bkg))
    ['AAAACTATCAAGCTAA',
     'ACAGTAGAGACGAGTA',
     'ACCTGACGTAATCGTC',
     'ACGATTACGTCAGGTA',
     'ATGAGATCGTAGCAAC',
     'ATGCTTAGTGCCATAC',
     'CAGTACGAAAAGGCCC',
     'CAGTAGAGACGAGTAA',
     'CCAGTACGAAAAGGCC',
     'GGTATGGCACTAAGCA',
     'GGTTGCTACGATCTCA',
     'TAAAACTATCAAGCTA']

    (6) __len__() - returns the number of k-mers in background

    >>> len(bkg)
    12

    (7) remove(seq) - removes all k-mers in seq from the
                      background, freeing them up for use

    >>> 'ATGCTTAGTGCCATACC' in bkg
    True
    >>> bkg.remove('ATGCTTAGTGCCATACC')
    >>> 'ATGCTTAGTGCCATACC' in bkg
    False
    >>> len(bkg)
    10

    (8) multiremove(seq_list) - removes all k-mers from given
                                seq_list from background

    >>> bkg.multiremove(my_background_list)

    [Background Processing]
      Removing Seq 4: TTAGCTTGAT
    >>> len(bkg)
    0

    (9) clear(Lmax=None) - removes all k-mers stored in background,
                           and optionally resets background Lmax

    >>> bkg.add('ATGCTTAGTGCCATACC')
    >>> len(bkg)
    2
    >>> bkg.clear()
    >>> len(bkg)
    0
    >>> bkg
    kmerSetDB stored at ./prj_bkg/ with 0 16-mers
    >>> bkg.clear(Lmax=17)
    >>> bkg
    kmerSetDB stored at ./prj_bkg/ with 0 18-mers
    >>> bkg.clear(Lmax=15)
    >>> bkg
    kmerSetDB stored at ./prj_bkg/ with 0 16-mers

    (10) close() - closes background instance; once closed
                   operations on the instance raises error
    
    >>> bkg.close()
    True
    >>> bkg.close()
    False
    >>> bkg.add('ATGCTTAGTGCCATACC')
    Traceback (most recent call last):
      File "<stdin>", line 1, in <module>
      File "nrpcalc/base/kmerSetDB.py", line 97, in wrapper
        raise RuntimeError('kmerSetDB was dropped')
    RuntimeError: kmerSetDB was closed or dropped

    (11) drop() - deletes unclosed `background` from disk;
                  once dropped operations on the instance
                  raises error

    >>> bkg.drop()
    False
    >>> bkg = nrpcalc.background(
            path='./prj_bkg/',
            Lmax=15)
    >>> bkg
    kmerSetDB stored at ./prj_bkg/ with 0 16-mers
    >>> bkg.drop()
    True
    >>> bkg.add('ATGCTTAGTGCCATACC')
    Traceback (most recent call last):
      File "<stdin>", line 1, in <module>
      File "nrpcalc/base/kmerSetDB.py", line 97, in wrapper
        raise RuntimeError('kmerSetDB was closed or dropped')
    RuntimeError: kmerSetDB was closed or dropped
    '''
    return kmerSetDB.kmerSetDB(
        path=path,
        homology=Lmax+1,
        verbose=verbose)

def finder(
    seq_list,
    Lmax,
    internal_repeats=False,
    background=None,
    vercov='nrp2',
    output_file=None,
    verbose=True):
    '''
    NRP Calculator Finder Mode for discovering non-repetitive
    subset of parts from a given list. All parts sharing any
    repeat longer than Lmax are eliminated from seq_list, and
    the approximately largest subset of non-repetitive parts
    is returned in a dictionary indexed by their position in
    seq_list. If internal_repeats is set to True, then parts
    with internal repeats are preserved, otherwise such parts
    are eliminated from seq_list. Optionally, the discovered
    subset of parts is written to an output FASTA file.

    :: seq_list
       type - list
       desc - a list of IUPAC strings representing a genetic
              part toolbox
    :: Lmax
       type - integer
       desc - maximum allowed shared repeat length between
              all sequences in a given toolbox
    :: internal_repeats
       type - boolean
       desc - if False then parts containing internal repeats
              longer than Lmax are eliminated; shared repeats
              are eliminated regardless
              (default=False)
    :: background
       type - kmerSetDB / None
       desc - the background object containing k-mers (k=Lmax+1)
              which must be absent in discovered non-repetitive
              subset of parts
              (default=None)
    :: vercov
       type - string
       desc - must be either '2apx', 'nrpG', or 'nrp2'
              '2apx' - use standard 2-approximation Vertex
                       Cover Elimination algorithm
              'nrpG' - use Greedy Vertex Cover Elimination
                       algorithm
              'nrp2' - use Finder Mode 2-approximation Vertex
                       Cover Elimination algorithm
              (default='nrp2')
    :: output_file
       type - string / None
       desc - filename to store discovered non-repetitive parts
              indexed by their position in seq_list; sequences
              are written in FASTA format
              (default=None)
    :: verbose
       type - boolean
       desc - if True displays progress
              (default=True)

    Returns: A dictionary of DNA or RNA strings with integer keys.

    Finder Mode API Example

    >>> import nrpcalc
    >>> 
    >>> # define background corpus
    >>> my_chromosomes = [
        'ATGAGATCGTAGCAACC',
        'GACGATTACGTCAGGTA',
        'ACAGTAGAGACGAGTAA',
        'CCAGTACGAAAAGGCCC',
        'AAAAAAAAAAAAAAAAA']
    >>> 
    >>> # initialize background
    >>> genomic_kmers = nrpcalc.background(
        path='./my_genome/',
        Lmax=15)
    >>> genomic_kmers.multiadd(
        my_chromosomes)

    [Background Processing]
      Adding Seq 4: AAAAAAAAAA...
    >>> 
    >>> # fetch part toolbox
    >>> my_toolbox = [
        'AGAGCTATGACTGACGT',
        'GCAGATAGGGGGTAGTA',
        'TAAAAAAAAAAAAAAAA', # Repeats with last chromosome
        'GAGCTATGACTGACGTC'] # Repeats with first part
    >>> 
    >>> # find non-repetitive subset
    >>> # with respect to background
    >>> nrpcalc.finder(
        seq_list=my_toolbox,
        Lmax=15,
        background=genomic_kmers)

    [Non-Repetitive Parts Calculator - Finder Mode]

    [Checking Constraints]
     Sequence List   : 4 parts
              Lmax   : 15 bp
     Internal Repeats: False

     Check Status: PASS

    [Checking Background]
     Background: kmerSetDB stored at ./my_genome/ with 10 16-mers

     Check Status: PASS

    [Checking Arguments]
       Vertex Cover: nrp2
       Output  File: None

     Check Status: PASS

    Extracted 4 unique sequences out of 4 sequences in 3.91e-05 seconds

    Written 4 unique sequences out to ./30c21235-e3f6-47f5-bce1-f99f47053e0b/seq_list.txt in 0.0002167 seconds

     [Sequence processing remaining] = 1 
     [Cliques inserted] = 2 

    Built homology graph in 0.0006673 seconds. [Edges = 1] [Nodes = 3]
     [Intital Nodes = 4] - [Repetitive Nodes = 1] = [Final Nodes = 3]

     [+] Initial independent set = 0, computing vertex cover on remaining 0 nodes.
     [+] Vertex Cover Function: NRP 2-approximation
     [+] Dumping graph into: ./30c21235-e3f6-47f5-bce1-f99f47053e0b/repeat_graph.txt in 0.00023055076599121094 seconds

    ----------------------
    Now running iteration: 0
    ----------------------

     Pendant checking is in progress...
      [+] 3 Pendants found

     Pendant elimination initiated...
      [x] Isolated node 1 eliminated
      [x] Pendant node 0 eliminated
      [+] Node 3 covered


     [+] Computed vertex cover of size: 1 (in 0.0002351 seconds)
     [+] Loading graph from: ./30c21235-e3f6-47f5-bce1-f99f47053e0b/repeat_graph.txt
     [+] Current independent set size:  2
     [+] Potential nodes for expansion: 0 (projected independent set size: 2)
     [X] Cannot expand independent set, terminating.

    Non-Repetitive Toolbox Size: 2
    {1: 'GCAGATAGGGGGTAGTA', 0: 'AGAGCTATGACTGACGT'}
    >>>
    >>> genomic_kmers.drop() # we're done with this background
    True
    '''
    return nrpfinder.nrp_finder(
        seq_list=seq_list,
        homology=Lmax+1,
        background=background,        
        allow_internal_repeat=internal_repeats,
        vercov_func=vercov,
        output_file=output_file,
        verbose=verbose,
        check_constraints=True)

def maker(
    seq_constr,
    struct_constr,
    part_type,
    Lmax,
    target_size,
    internal_repeats=False,
    background=None,
    struct_type='mfe',
    seed=None,
    synth_opt=False,
    local_model_fn=None,
    global_model_fn=None,
    jump_count=10,
    fail_count=1000,
    output_file=None,
    verbose=True):
    '''
    NRP Calculator Maker Mode for designing non-repetitive
    genetic part toolboxes from user defined sequence and
    structure constraints and based on custom local and/or
    global model functions. All shared repeats longer than
    Lmax are eliminated, and internal repeats longer than
    Lmax are preserved if desired. Parts are optimized for
    DNA synthesis if desired. Error tolerance is adaptive
    and auto-adjusted based on recorded failures. Designed
    toolbox is returned as a dictionary of parts indexed by
    their order of design, and optionally written to a FASTA
    output file.

    :: seq_constr
       type - string
       desc - a string in IUPAC degenerate code describing
              all valid nucleotide choices at each position
              e.g. 'NNNNWWWWSSSSTTTT' implies that the first
                   four bases can be either 'A'/'T'/'G'/'C',
                   the next four bases can be either 'A'/'T',
                   followed by either 'G'/'C' for the next
                   four basses, and finally ending with 'T's
    :: struct_constr
       type - string / None
       desc - a string in dot-parenthesis-x notation that
              describe the secondary base pairing across all
              nucleotide positions; a string of only dots implies
              an absence of any structure constraint
              e.g. '..((xx))..' implies that the first, second,
                   and the last two bases are free to either
                   base pair or not ('.'), the third and fourth
                   bases are paired with the eighth and the
                   seventh bases respectively ('(' and ')'),
                   while the fifth and the sixth base must
                   not take part in any base pairing ('x') at all
    :: part_type
       type - string
       desc - must be either 'RNA' or 'DNA' depending on the
              type of genetic part being designed; ensures
              that correct folding free-energy parameters are
              used during structure evaluation
    :: Lmax
       type - integer
       desc - maximum allowed shared repeat length between
              all sequences in designed toolbox
    :: target_size
       type - integer
       desc - maximum number of genetic parts to be designed
              for the generated toolbox; target_size may not
              be reached if the constraints are too strict,
              for example, due to low degeneracy in the given
              sequence constraint, or a low Lmax
    :: internal_repeats
       type - boolean
       desc - if True then internal repeats in designed parts
              are not eliminated; useful when designing parts
              such as rho-independent terminators with structure
              constraints that necessitate internal repeats;
              shared repeats are eliminated regardless
              (default=False)
    :: background
       type - kmerSetDB / None
       desc - a background object containing k-mers (k=Lmax+1)
              which must be absent in the designed toolbox
              (default=None)
    :: struct_type
       type - string
       desc - must be either 'mfe', 'centroid', or 'both'
              'mfe' - use minimum free energy structure evaluation
              'centroid' - use centroid structure evaluation
              'both' - use both 'mfe' and 'centroid' evaluation
              (default='mfe')
    :: seed
       type - integer / None
       desc - integer used to seed random number generations;
              two Maker runs with same constraints and seed
              value will generate the exact same toolbox;
              if None then a random seed value is used
              (default=None)
    :: synth_opt
       type - boolean
       desc - if True then designed parts containing features
              that complicate DNA synthesis are eliminated
              (default=False)
    :: local_model_fn
       type - function / None
       desc - a function with signature 'fn_name(seq)' that
              takes in a partial genetic part sequence, and
              returns a tuple (state, index), where state is
              either True or False and index is a traceback
              index / location or None depending on whether
              a custom design objective was met or not;
              useful for providing concurrent feedback to the
              path-finding process and steering nucleotide
              selection choices
              e.g. prevent_cutsites(seq) maybe be a local
                   function that takes in a partial sequence
                   as it is built and returns (True, None)
                   if the last six bases of the partial
                   seq[-6:] is not the same as any of the
                   cutsites used for cloning the part, else
                   returns a tuple (False, len(seq)-6) as
                   the traceback location to reselect the
                   last six bases; this naturally ensures
                   the final part is devoid of any cutsites
                   used experimentally throughout the part
              (default=None)
    :: global_model_fn
       type - function / None
       desc - a function with signature 'fn_name(seq)' that
              takes in a complete genetic part sequence, and
              returns either True or False depending on whether
              a custom design objective was met; useful for
              design criteria that can only be evaluated when
              the complete genetic part is available; parts
              that are evaluated to be False are rejected and
              a new part generation is started; global model
              functions are evaluated only after the last base
              has been added to a genetic part under design
              e.g. gc_content(seq) may be a global function
                   that takes in a complete sequence and only
                   accepts parts with GC content greater than
                   threshold percentage (although, technically
                   one can enforce this condition via a well
                   planned local model function)
              (default=None)
    :: jump_count
       type - integer
       desc - maximum number of restarts in path finding due
              to failure in finding suitable k-mers, meeting
              local_model_fn, or being stuck in local optima
              (auto-adjusted with each iteration)
              (default=10)
    :: fail_count
       type - integer
       desc - maximum number of consecutive failures tolerated
              when structure constraints, global model functions
              and synthesis objectives are not met
              (auto-adjusted with each iteration)
              (default=1000)
    :: output_file
       type - string / None
       desc - filename to store designed non-repetitive parts
              as they are generated consecutively; sequences
              are written in FASTA format
              (default=None)
    :: verbose
       type - boolean
       desc - if True displays progress
              (default=True)

    Returns: A dictionary of DNA or RNA strings with integer keys.
    
    Maker Mode API Example
    
    >>> import nrpcalc
    >>>
    >>> Lmax = 15  # for global use
    >>>
    >>> # initialize background
    >>> bkg = nrpcalc.background(
        path='./my_toolbox_kmers',
        Lmax=Lmax)
    >>> bkg
    kmerSetDB stored at ./my_toolbox_kmers/ with 0 16-mers
    >>>
    >>> # a local model function
    >>> def prevent_cutsites(seq):
            # part is long enough to evaluate
            if len(seq) >= 6:
                BamHI = 'GGATCC' # cutsite 1
                XbaI  = 'TCTAGA' # cutsite 2
                # evaluate part concurrently
                if seq[-6:] in [BamHI, XbaI]:
                    return (False, len(seq)-6)
                else:
                    return (True, None)
            # part is short enough .. pass
            return (True, None)
    >>>
    >>> # a global model function
    >>> def optimal_gc_content(seq):
            # compute gc count
            gcount = seq.count('G') * 1.
            ccount = seq.count('C') * 1.
            # evaluate on complete part
            if 0.4 <= (gcount + ccount) / len(seq) <= 0.6:
                return True
            else:
                return False
    >>>
    >>> # final toolbox storage
    >>> final_toolbox = []
    >>>
    >>> # build a toolbox of non-repetitive sigma70 promoters
    >>> promoters_strong = nrpcalc.maker(
        seq_constr='N'*20+'TTGACA'+'N'*17+'TATAAT'+'NNNNNN',
        struct_constr='.'*55,
        part_type='DNA',
        Lmax=Lmax,
        target_size=500,
        internal_repeats=False,
        background=None,
        local_model_fn=prevent_cutsites,
        global_model_fn=optimal_gc_content)
    WARNING: stacking enthalpies not symmetric
    WARNING: stacking enthalpies not symmetric
    WARNING: stacking enthalpies not symmetric
    WARNING: stacking enthalpies not symmetric

    [Non-Repetitive Parts Calculator - Maker Mode]

    [Checking Constraints]
      Sequence Constraint: NNNNNNNNNNNNNNNNNNNNTTGACANNNNNNNNNNNNNNNNNTATAATNNNNNN
     Structure Constraint: .......................................................
          Part Type      : DNA
               Lmax      : 15 bp
        Target Size      : 500 parts
      Internal Repeats   : False

     Check Status: PASS

    [Checking Arguments]
     Struct Type : mfe
      Synth Opt  : False
       Jump Count: 10
       Fail Count: 1000
     Output File : None

     Check Status: PASS

    Constructing Toolbox:

     [part] 1, [16-mers] 40, [iter 1] 0.00s, [avg] 0.00s, [total time] 0.00h
     [part] 2, [16-mers] 80, [iter 2] 0.01s, [avg] 0.01s, [total time] 0.00h
     ...
     [part] 499, [16-mers] 19960, [iter 499] 0.00s, [avg] 0.00s, [total time] 0.00h
     [part] 500, [16-mers] 20000, [iter 500] 0.00s, [avg] 0.00s, [total time] 0.00h

    Construction Complete.

    Non-Repetitive Toolbox Size: 500
    >>>
    >>> # add promoter toolbox to background
    >>> bkg.multiadd(promoters_strong.values())

    [Background Processing]
      Adding Seq 499: AACATCGATG... 
    >>>
    >>> # add promoters to final_toolbox
    >>> final_toolbox.extend(promoters_strong.values())
    >>>
    >>> # build another toolbox of sigma70 promoters
    >>> # non-repetitive to the previous toolbox
    >>> # (this is called toolbox chaining)
    >>> promoters_variable = nrpcalc.maker(
        seq_constr='N'*20+'TTGACA'+'N'*16+'WWWWWWW'+'NNNNN',
        struct_constr='.'*55,
        part_type='DNA',
        Lmax=Lmax,
        target_size=500,
        internal_repeats=False,
        background=bkg,
        local_model_fn=prevent_cutsites,
        global_model_fn=optimal_gc_content)
    WARNING: stacking enthalpies not symmetric
    WARNING: stacking enthalpies not symmetric
    WARNING: stacking enthalpies not symmetric
    WARNING: stacking enthalpies not symmetric

    [Non-Repetitive Parts Calculator - Maker Mode]

    [Checking Constraints]
      Sequence Constraint: NNNNNNNNNNNNNNNNNNNNTTGACANNNNNNNNNNNNNNNNWWWWWWWNNNNN
     Structure Constraint: .......................................................
        Target Size      : 500 parts
               Lmax      : 15 bp
      Internal Repeats   : False

     Check Status: PASS

    [Checking Background]
     Background: kmerSetDB stored at ./my_toolbox_kmers/ with 20000 16-mers

     Check Status: PASS

    [Checking Arguments]
       Part Type : DNA
     Struct Type : mfe
      Synth Opt  : False
       Jump Count: 10
       Fail Count: 1000
     Output File : None

     Check Status: PASS

    Constructing Toolbox:

     [part] 1, [16-mers] 39, [iter 1] 0.01s, [avg] 0.01s, [total time] 0.00h
     [part] 2, [16-mers] 78, [iter 2] 0.00s, [avg] 0.01s, [total time] 0.00h
     ...
     [part] 499, [16-mers] 19461, [iter 499] 0.00s, [avg] 0.00s, [total time] 0.00h
     [part] 500, [16-mers] 19500, [iter 500] 0.00s, [avg] 0.00s, [total time] 0.00h

    Construction Complete.

    Non-Repetitive Toolbox Size: 500
    >>>
    >>> # remove background .. it's served its purpose
    >>> bkg.drop()
    True
    >>>
    >>> # add new promoters to final_toolbox
    >>> final_toolbox.extend(promoters_variable.values())
    >>>
    >>> # verify all promoters are non-repetitive
    >>> assert len(nrpcalc.finder(
        seq_list=final_toolbox,
        Lmax=Lmax,
        verbose=False)) == 1000
    '''
    # Initialize object
    _maker = nrpmaker.NRPMaker(
        part_type=part_type,
        seed=seed)
    
    # Compute parts
    parts = _maker.nrp_maker(
        seq_constr=seq_constr,
        struct_constr=struct_constr,
        struct_type=struct_type,
        target_size=target_size,
        homology=Lmax+1,
        allow_internal_repeat=internal_repeats,
        background=background,
        synth_opt=synth_opt,
        local_model_fn=local_model_fn,
        global_model_fn=global_model_fn,
        jump_count=jump_count,
        fail_count=fail_count,
        output_file=output_file,
        verbose=verbose,
        abortion=True,
        check_constraints=True)

    # Cleanups and return
    del _maker
    return parts