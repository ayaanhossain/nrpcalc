from .base import maker     as nrpmaker
from .base import finder    as nrpfinder
from .base import kmerSetDB


__version__ = '1.1.0'

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

NRP Calculator offers two modes of operation:

- Finder Mode: Discover toolboxes of non-repetitive parts
               from a list of candidate parts

-  Maker Mode: Design toolboxes of non-repetitive parts
               based on sequence, structure and model
               constraints

Additionally, a 'background' object is available which can
be used to store background sequences against which parts
discovered or designed are ensured to be non-repetitive.

You can learn more about the two modes and background via
  print nrpcalc.background.__doc__
  print nrpcalc.finder.__doc__
  print nrpcalc.maker.__doc__
'''.format(
    __version__,
    *__authors__.strip().split('\n'))


def background(
    path,
    Lmax,
    verbose=True):
    '''
    NRP Calculator kmerSetDB 'background' object for on-disk
    storage of background sequence k-mers. When a sequence is
    added to background, k-mers from the sequence is added in
    instead of the actual sequence itself. A sequence queried
    for existence in background, is evaluated to be True if
    any k-mer from the sequence exists in the background.
    Useful in chaining multiple Maker and Finder jobs as well
    as persisting any backgrounds such as small genomes or
    other part toolboxes from earlier design rounds.

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

    Note that if the path provided points to an existing
    background object, then that background is opened for
    reading, otherwise, a new background is instantiated
    at the given path.

    Background / kmerSetDB API Examples
    
    >>> import nrpcalc
    >>>
    >>> my_background_list = [...]
    >>> bkg = nrpcalc.background(
           path='./my_backgound/',
           Lmax=15)

    (1) add(seq) - adds an IUPAC string 'seq' to background

    >>> bkg.add('ATGCTAGGCCAACC')
    
    (2) multiadd(seq_list) - adds multiple sequences in
                             the list to background

    >>> bkg.multiadd(my_background_list)
    
    (3) __contains__(seq) - check if all k-mers from 'seq'
                            is present in background

    >>> 'ATGCTAGGCCAACC' in bkg

    (4) multicheck(seq_list) - check if all k-mers from given
                               seq_list present in background

    >>> assert all(bkg.multicheck(my_background_list))

    (5) __iter__() - iterate over all k-mers in background

    >>> e.g. kmers = list(iter(bkg))

    (6) remove(seq) - removes all k-mers in 'seq' from the
                      background, freeing them up for use

    >>> e.g. bkg.remove('ATGCTAGGCCAACC')

    (7) multiremove(seq_list) - removes all k-mers in the given
                                seq_list from background

    >>> e.g. bkg.multiremove(my_background_list)

    (8) clear() - removes all k-mers stored in background

    >>> bkg.clear()

    (9) close() - closes background instance
    
    >>> bkg.close()

    (10) drop() - deletes background instance from disk

    >>> bkg.drop()
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
    are eliminated from seq_list.

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
              longer than Lmax are eliminated
              (default=False)
    :: background
       type - kmerSetDB / None
       desc - the background object containg k-mers (k=Lmax+1)
              which must be absent in returned non-repetitive
              subset of parts
              (default=None)
    :: vercov
       type - string
       desc - must be either '2apx', 'nrpG', or 'nrp2';
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

    Finder Mode API Example

    >>> import nrpcalc
    >>>
    >>> my_background_list = [...]
    >>> bkg = nrpcalc.background(...)
    >>> bkg.multiadd(my_background_list)
    >>>
    >>> my_parts = [...]
    >>> nrpset = nrpcalc.finder(
    ...     seq_list=my_parts,
    ...     Lmax=15,
    ...     internal_repeats=False,
    ...     background=bkg,
    ...     vercov='nrp2')
    >>>
    '''
    return nrpfinder.nrp_finder(
        seq_list=seq_list,
        homology=Lmax+1,
        background=background,        
        allow_internal_repeat=internal_repeats,
        vercov_func=vercov,
        verbose=verbose)

def maker(
    seq_constr,
    struct_constr,
    target_size,
    Lmax,
    internal_repeats=False,
    background=None,
    part_type='RNA',
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
    genetic part toolboxes from user defined constraints
    such as sequence and structure constraints and ensured
    to satisfy specified local and global model functions.
    The designed toolbox is returned as a dictionary of
    parts indexed by their order of design, and optionally
    written to a FASTA output file.

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
       type - string
       desc - a string in dot-parenthesis-x notation that
              describe the secondary base pairing across all
              nucleotide positions
              e.g. '..((xx))..' implies that the first, second,
                   and the last two bases are free to either
                   base pair or not (dot), the third and fourth
                   bases are paired with the eighth and the
                   seventh bases respectively (parenthesis),
                   while the fifth and the sixth base must
                   not  take part in any base pairing (x)
    :: target_size
       type - integer
       desc - maximum number of genetic parts to be designed
              for the generated toolbox; target_size may not
              be reached if the constraints are too strict,
              for example, due to low degeneracy in the given
              sequence constraint
    :: Lmax
       type - integer
       desc - maximum allowed shared repeat length between
              all sequences in a given toolbox
    :: internal_repeats
       type - boolean
       desc - if True then designed parts are not eliminated
              due to internal repeats; useful when designing
              parts such as rho-independent terminators with
              structure constraints that necessitate internal
              repeats; shared repeats are always eliminated
              (default=False)
    :: background
       type - kmerSetDB / None
       desc - the background object containg k-mers (k=Lmax+1)
              which must be absent in the designed toolbox
              (default=None)
    :: part_type
       type - string
       desc - must be either 'RNA' or 'DNA' depending on the
              type of genetic part being designed; ensures
              that correct free-energy folding parameters are
              used during structure evaluation
              (default='RNA')
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
              that complicate synthesis are eliminated
              (default=False)
    :: local_model_fn
       type - function / None
       desc - a function with signature 'fn_name(seq)' that
              takes in a partial genetic part sequence, and
              returns a tuple (state, index), where state is
              either True or False and index is a traceback
              index / location or None depending on whether
              a custom design objective was met or not;
              useful for providing direct feedback to the
              path-finding process by steering nucleotide
              selection choices
              e.g. prevent_cutsites(seq) maybe be a local
                   function that takes in a partial sequence
                   as it is built and returns (True, None)
                   if the last six bases of the partial
                   part[-6:] is not the same as any of the
                   cutsites used for cloning the part, else
                   returns a tuple (False, len(seq)-6) as
                   the traceback location to reselect the
                   last six bases; this naturally ensures
                   the final part is devoid of all cutsites
                   throughout the part
              (default=None)
    :: global_model_fn
       type - function / None
       desc - a function with signature 'fn_name(seq)' that
              takes in a complete genetic part sequence, and
              returns either True or False depending on whether
              a custom design objective was met; useful for
              design criteria that can only be evlauated when
              the complete genetic part is available; parts
              that are evaluated to be False are rejected and
              a new part generation is started; global model
              functions are evaluated after the last base has
              been added to the genetic part being evaluated
              e.g. gc_content(seq) may be a global function
                   that takes in a complete sequence and only
                   accepts parts with GC content greater than
                   threshold percentage
              (default=None)
    :: jump_count
       type - integer
       desc - maximum number of restarts in path finding due
              to failure in meeting local_model_fn or getting
              stuck in a local optima
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
    
    Maker Mode API Example
    
    >>> import nrpcalc
    >>>
    >>> my_background_list = [...]
    >>> bkg = nrpcalc.background(...)
    >>> bkg.multiadd(my_background_list)
    >>>
    >>> def prevent_cutsites(seq):
    ...     BamHI = 'GGATCC'
    ...     XbaI  = 'TCTAGA'
    ...     if seq[-6:] in [BamHI, XbaI]:
    ...         return (False, len(seq)-6)
    ...     else:
    ...         return (True, None)
    >>>
    >>> def optimal_gc_content(seq):
    ...     gcount = seq.count('G') * 1.
    ...     ccount = seq.count('C') * 1.
    ...     if 0.4 <= (gcount + ccount) / len(seq) <= 0.6:
    ...         return True
    ...     else:
    ...         return False 
    >>>
    >>> promoters_1 = nrpcalc.maker(
    ...     seq_constr='N'*20+'TTGACA'+'N'*17+'TATAAT'+'NNNNNN',
    ...     struct_constr='.'*55,
    ...     target_size=500,
    ...     Lmax=15,
    ...     internal_repeats=False,
    ...     background=bkg,
    ...     part_type='DNA',
    ...     local_model_fn=prevent_cutsites,
    ...     global_model_fn=optimal_gc_content)
    >>>
    '''
    _maker = nrpmaker.NRPMaker(
        part_type=part_type,
        seed=seed)
    
    return _maker.nrp_maker(
        seq_constr=seq_constr,
        part_type=part_type,
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
        abortion=True)