from .base import maker
from .base import finder
from .base import kmerSetDB


__version__ = '1.0.0'
__authors__ = '''Ayaan Hossain <auh57@psu.edu>
\r\t  Howard Salis <salis@psu.edu>'''

__doc__ = '''
Non-Repetitive Parts Calculator

Automated design and discovery of non-repetitive 
genetic parts for engineering stable genetic systems

Version = {}
Authors = {}

NRP Calculator offers two modes of operations
- Finder Mode: Discover toolboxes of non-repetitive parts
               from a list of candidate sequences
-  Maker Mode: Automated design of non-repetitive parts
               based on constraints

Additionally, a 'background' object is present which can
be used to store background sequences against which parts
discovered to designed are ensured to be non-repetitive

You can learn more about the modes and background via
  print nrpcalc.background.__doc__
  print nrpcalc.finder.__doc__
  print nrpcalc.maker.__doc__
'''.format(
    __version__,
    __authors__)


def background(
    path,
    Lmax,
    verbose=True):
    '''
    NRP Calculator 'background' object for on-disk storage of background
    sequences, as well as any sequences built using Maker. Useful when
    chaining several Maker and Finder jobs as well as persisting background
    developed from previous runs.

    :: path
       type - string
       desc - ./a/path/to/store/background/object
    :: Lmax
       type - integer
       desc - maximum allowed shared repeat length
    :: verbose
       type - boolean
       desc - if True displays progress
              (default=True)
    '''
    return kmerSetDB(
        path=path,
        homology=Lmax+1,
        verbose=verbose)

def finder(
    seq_list,
    Lmax,
    internal_repeats=False,
    background=None,
    vercov_func=None,
    verbose=True):
    '''
    Holder function for Finder Mode. Given a list of sequences, and a Lmax length,
    finds the largest possible set of sequences that are non-repetitive, i.e. they
    do not share a subsequences greater than or equal to Lmax length.

    seq_list        : a list of strings in the alphabet {A, T, C, G}
    Lmax            : an integer, such that found sequences are non-repetitive over that length
    internal_repeats: optional boolean, if true allows internal repeats in sequences, else selects sequences without internal repeats
    background      : a kmerSetDB object storing strings against which found sequences must be non-repetitive
    vercov_func     : optional string, can be '2apx' (classical 2-approximation) / 'dosG' (greedy DOS Vercov wtih Recovery) / None (default, 2-Approx DOS Vercov with recovery)
    verbose         : optional boolean, enables/disables printing of progress
    '''
    return finder.nrp_finder(
        seq_list=seq_list,
        homology=Lmax+1,
        background=background,        
        internal_repeats=internal_repeats,
        vercov_func=vercov_func,
        verbose=verbose)

def maker(
    seq_list,
    struct_list,
    target_list,
    Lmax,
    internal_repeats=False,
    background=None,
    struct_type=None,
    seed=None,
    synth_opt=True,
    local_model_fn=None,
    global_model_fn=None,
    jump_count=10,
    fail_count=1000,
    output_file=None,
    verbose=False,
    abortion=True):
    '''
    Holder function for Maker Mode. Given a list of sequence and structure constraints, their 
    target counts, and a Lmax length, attempts to design the largest possible set of 
    sequences that are non-repetitive, while optionally optimizing for synthesizability 
    and user specified model functions, within some specified jump and failure counts.

    seq_list        : a list of strings in IUPAC degenerate code alphabet
    struct_list     : a list of strings in the alphabet {(, ), ., x}, can be None or a zero length string if no structure desired
    target_list     : a list of integers specifying the target number of parts
    Lmax            : an integer, such that made sequences are non-repetitive over that length
    internal_repeats: optional boolean, if true allows internal repeats in sequences, else builds sequences without internal repeats
    background      : a kmerSetDB object storing strings against which made sequences must be non-repetitive
    struct_type     : optional string, 'mfe', 'centroid', 'both' or None (default None, uses mfe if structure given)
    seed            : optional integer, helps the maker reproduce results by seeding the random number generator (default None)
    synth_opt       : optional boolean, to enable/disable synthesis optimization (default True)
    local_model_fn  : optional user supplied function, which must take a partial sequence as input, and produce a boolean tuple (boolean, index), enables local tracebacks
    global_model_fn : optional user supplied function, which must take a sequence as input and produce a boolean as output, enables global model optimization
    jump_count      : optional integer, a maximum count of jumps allowed to the candidate sequence generator, before following abortion parameter
    fail_count      : optional integer, a maximum count of consecutive failures in structure, synthesis, and model function satisfiability to be tolerated before termination
    output_file     : optional string, all generated strings are written to this output file in FASTA format
    verbose         : optional boolean, enables/disables printing of progress
    abortion        : optional boolean to disable traceback jumps, not recommeded to set it to False (default True)
    '''
    _maker = maker.NRPMaker(seed)
    return _maker.nrp_maker(
        seq_list=seq_list,
        struct_list=struct_list,
        struct_type=struct_type,
        target_list=target_list,
        homology=Lmax+1,
        allow_internal_repeat=internal_repeats,
        background_list=background_list,        
        synth_opt=synth_opt,
        local_model_fn=local_model_fn,
        global_model_fn=global_model_fn,
        jump_count=jump_count,
        fail_count=fail_count,
        output_file=output_file,
        verbose=verbose,
        abortion=abortion)
