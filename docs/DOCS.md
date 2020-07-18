

<h1 align="center">
    <a href="https://github.com/ayaanhossain/nrpcalc/">
        <svg>
            <img src="https://raw.githubusercontent.com/ayaanhossain/nrpcalc/master/img/logo.svg?sanitize=true"  alt="Non-Repetitive Parts Calculator" width="418" class="center"/>
        </svg>
    </a>
</h1>

<p align="center">
  <a href="#Background">Background</a> •
  <a href="#Finder-Mode">Finder Mode</a> •
  <a href="#Maker-Mode">Maker Mode</a> •
  <a href="../README.md">README</a> •
  <a href="../README.md">EXAMPLES</a>
</p>

## Background

`NRP Calculator` `kmerSetDB` `background` object for on-disk storage of background sequence _k_-mers (where _k_=`Lmax`+1). When a sequence is added to `background`, _k_-mers from the sequence are added instead of the actual sequence itself. A sequence queried for existence in the given `background` is evaluated to be `True` if any _k_-mer from the sequence exists in the background. This object is useful when chaining multiple `Maker Mode` and `Finder Mode` jobs as well as persisting any backgrounds such as small genomes or other part toolboxes from earlier design rounds.

**nrpcalc.background(path, Lmax, verbose=True)**

| argument | type | description | default |
|--|--|--|--|
| `path` | `string` | ./a/path/to/store/background/object for later reuse and part verification | -- |
| `Lmax` | `integer` | maximum allowed shared repeat length between all sequences in a given toolbox | -- |
| `verbose` | `boolean` | if `True` displays progress | `True` |

**_Returns_**: A `kmerSetDB` object.

**_Note:_** If the path provided points to an existing `background` object, then that `background` is opened for reading, and the new `Lmax` is ignored, otherwise, a new `background` is instantiated at the given path.

`background` / `kmerSetDB` **API Examples**

```python
>>> from pprint import pprint
>>> import nrpcalc
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
```

`background` / `kmerSetDB` offers the following methods

(1) **add(seq)** - adds an IUPAC string `seq` to `background`
```python
>>> bkg.add('ATGCTTAGTGCCATACC')
```
(2) **multiadd(seq_list)** - adds multiple sequences from the list to `background`
```python
>>> bkg.multiadd(my_background_list)

[Background Processing]
  Adding Seq 4: TTAGCTTGAT...
```
(3) **\_\_contains__(seq)** - checks if all _k_-mers from `seq`  is present in `background`
```python
>>> 'ATGCTTAGTGCCATACC' in bkg
True
```
(4) **multicheck(seq_list)** - checks if all _k_-mers from given `seq_list` present in `background`
```python
>>> assert all(bkg.multicheck(my_background_list))
```
(5) **\_\_iter__()** - iterates over all _k_-mers in `background`
```python
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
```
(6) **\_\_len()__** - returns the number of _k_-mers in `background`
```python
>>> len(bkg)
12
```
(7) **remove(seq)** - removes all _k_-mers in seq from the `background`, freeing them up for use
```python
>>> 'ATGCTTAGTGCCATACC' in bkg
True
>>> bkg.remove('ATGCTTAGTGCCATACC')
>>> 'ATGCTTAGTGCCATACC' in bkg
False
>>> len(bkg)
10
```
(8) **multiremove(seq_list)** - removes all _k_-mers from given seq_list from `background`
```python
>>> bkg.multiremove(my_background_list)

[Background Processing]
  Removing Seq 4: TTAGCTTGAT
>>> len(bkg)
0
```
(9) **clear()** - removes all _k_-mers stored in `background`
```python
>>> bkg.add('ATGCTTAGTGCCATACC')
>>> len(bkg)
2
>>> bkg.clear()

[Background Processing]
  Removing Seq 1: GGTATGGCAC...
>>> len(bkg)
0
```
(10) **close()** - closes `background` instance; once closed operations on the instance raises error
```python
>>> bkg.close()
True
>>> bkg.close()
False
>>> bkg.add('ATGCTTAGTGCCATACC')
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
  File "nrpcalc/base/kmerSetDB.py", line 96, in wrapper
    raise RuntimeError('kmerSetDB was dropped')
RuntimeError: kmerSetDB was closed or dropped
```
(11) **drop()** - deletes unclosed `background` from disk; once dropped operations on the instance raises error
```python
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
  File "nrpcalc/base/kmerSetDB.py", line 96, in wrapper
    raise RuntimeError('kmerSetDB was closed or dropped')
RuntimeError: kmerSetDB was closed or dropped
```

## Finder Mode

`NRP Calculator` `Finder Mode` for discovering non-repetitive subset of parts from a given list. All parts sharing any repeat longer than `Lmax` are eliminated from `seq_list`, and the approximately largest subset of non-repetitive parts is returned in a dictionary indexed by their position in `seq_list`. If `internal_repeats` is set to `True`, then parts with internal repeats are preserved, otherwise such parts are eliminated from `seq_list`. Optionally, the discovered subset of parts is written to an output `FASTA` file.

<h3 align="center">
    <a href="https://github.com/ayaanhossain/nrpcalc/img/Fig2.svg">
        <svg>
            <img src="https://raw.githubusercontent.com/ayaanhossain/nrpcalc/master/img/Fig2.svg?sanitize=true"  alt="NRP Calculator Finder Mode Algorithm" width="800" class="center"/>
        </svg>
    </a>
</h3>

**nrpcalc.finder(seq_list, Lmax, internal_repeats=False, background=None, vercov='nrp2', output_file=None, verbose=True)**

| argument | type | description | default |
|--|--|--|--|
| `seq_list` | `list` | a list of IUPAC strings representing a genetic part toolbox | -- |
| `Lmax` | `integer` | maximum allowed shared repeat length between all sequences in a given toolbox | -- |
| `internal_repeats` | `boolean` | if `False` then parts containing internal repeats longer than `Lmax` are eliminated; shared repeats are eliminated regardless | `False` |
| `background` | `kmerSetDB`/`None` | the `background` object containing _k_-mers (_k_=`Lmax`+1) which must be absent in discovered non-repetitive subset of parts | `None` |
| `vercov` | `string` | must be either `'2apx'`, `'nrpG'`, or `'nrp2'` <br> `'2apx'` - use standard 2-approximation Vertex Cover Elimination algorithm <br> `'nrpG'` - use Greedy Vertex Cover Elimination algorithm <br> `'nrp2'` - user `Finder Mode` 2-approximation Vertex Cover Elimination Algorithm | `'nrp2'` |
| `output_file` | `string`/`None` | filename to store discovered non-repetitive parts indexed by their position in `seq_list`; sequences are written in `FASTA` format | `None` |
| `verbose` | `boolean` | if `True` displays progress | `True` |

**_Returns_**: A `dictionary` of IUPAC strings with integer keys.

`Finder Mode` **API Examples**

```python
>>> import nrpcalc
>>> 
>>> my_chromosomes = [
    'ATGAGATCGTAGCAACC',
    'GACGATTACGTCAGGTA',
    'ACAGTAGAGACGAGTAA',
    'CCAGTACGAAAAGGCCC',
    'AAAAAAAAAAAAAAAAA']
>>> 
>>> genomic_kmers = nrpcalc.background(
    path='./my_genome/',
    Lmax=15)
>>> 
>>> genomic_kmers.multiadd(
    my_chromosomes)

[Background Processing]
  Adding Seq 4: AAAAAAAAAA...
>>> 
>>> my_toolbox = [
    'AGAGCTATGACTGACGT',
    'GCAGATAGGGGGTAGTA',
    'TAAAAAAAAAAAAAAAA', # Repeats with last chromosome
    'CAGATGATGCTAGGACT']
>>> 
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

[Checking Background]:
 Background: kmerSetDB stored at ./my_genome/ with 10 16-mers

 Check Status: PASS

[Checking Arguments]
   Vertex Cover: nrp2
   Output  File: None

 Check Status: PASS

Extracted 4 unique sequences out of 4 sequences in 1.693e-05 seconds

Written 4 unique sequences out to ./5ebb5779-5314-41ce-8114-2d858ef41e2e/seq_list.txt in 9.68e-05 seconds

 [Sequence processing remaining] = 1 
 [Cliques inserted] = 3 

Built homology graph in 0.000263 seconds. [Edges = 0] [Nodes = 3]
 [Intital Nodes = 4] - [Repetitive Nodes = 1] = [Final Nodes = 3]

 [+] Initial independent set = 0 [0 completex], computing vertex cover on remaining 3 nodes.
 [+] Vertex Cover Function: NRP 2-approximation
 [+] Dumping graph into: ./5ebb5779-5314-41ce-8114-2d858ef41e2e/repeat_graph.txt in 0.000402927398682 seconds

----------------------
Now running iteration: 0
----------------------

 Pendant checking is in progress...
  [+] 3 Pendants found

 Pendant elimination initiated...
  [x] Isolated node 0 eliminated
  [x] Isolated node 1 eliminated
  [x] Isolated node 3 eliminated


 [+] Computed vertex cover of size: 0 (in 0.000123 seconds)
 [+] Loading graph from: ./5ebb5779-5314-41ce-8114-2d858ef41e2e/repeat_graph.txt
 [+] Current independent set size:  3
 [+] Potential nodes for expansion: 0 (projected independent set size: 3)
 [X] Cannot expand independent set, terminating.

Non-Repetitive Toolbox Size: 3
{0: 'AGAGCTATGACTGACGT', 1: 'GCAGATAGGGGGTAGTA', 3: 'CAGATGATGCTAGGACT'}
```

## Maker Mode

`NRP Calculator` `Maker Mode` for designing non-repetitive genetic part toolboxes from user defined sequence and structure constraints and based on custom local and/or global model functions. All shared repeats longer than `Lmax` are eliminated, and internal repeats longer than `Lmax` are preserved if desired. Parts are optimized for `DNA` synthesis if desired. Error tolerance is adaptive and auto-adjusted based on recorded failures. Designed toolbox is returned as a `dictionary` of parts indexed by their order of design, and optionally written to a `FASTA` output file.

<h3 align="center">
    <a href="https://github.com/ayaanhossain/nrpcalc/img/Fig3.svg">
        <svg>
            <img src="https://raw.githubusercontent.com/ayaanhossain/nrpcalc/master/img/Fig3.svg?sanitize=true"  alt="NRP Calculator Finder Mode Algorithm" width="800" class="center"/>
        </svg>
    </a>
</h3>

**nrpcalc.maker(seq_constr, struct_constr, target_size, Lmax, internal_repeats=False, background=None, part_type='RNA', struct_type='mfe', seed=None, synth_opt=False, local_model_fn=None, global_model_fn=None, jump_count=10, fail_count=1000, output_file=None, verbose=True)**

| argument | type | description | default |
|--|--|--|--|
| `seq_constr` | `string` | a string in IUPAC degenerate code describing all valid nucleotide choices at each position <br> **e.g.** `'NNNNWWWWSSSSTTTT'` implies that the first four bases can be either `'A'`/`'T'`/`'G'`/`'C'`, the next four bases can be either `'A'`/`'T'`, followed by either `'G'`/`'C'` for the next four basses, and finally ending with `'T'`s | -- |
| `struct_constr` | `string` | a string in `dot-parenthesis-x` notation that describe the secondary base pairing across all nucleotide positions <br> **e.g.** `'..((xx))..'` implies that the first, second, and the last two bases are free to either base pair or not (`dot`), the third and fourth bases are paired with the eighth and the seventh bases respectively (`parenthesis`), while the fifth and the sixth base must not take part in any base pairing (`x`) at all | -- |
| `target_size` | `integer` | maximum number of genetic parts to be designed for the generated toolbox; `target_size` may not be reached if the constraints are too strict, for example, due to low degeneracy in the given sequence constraint, or a low `Lmax` | -- |
| `Lmax` | `integer` | maximum allowed shared repeat length between all sequences in designed toolbox | -- |
| `internal_repeats` | `boolean` | if `True` then internal repeats in designed parts are not eliminated; useful when designing parts such as rho-independent terminators with structure constraints that necessitate internal repeats; shared repeats are always eliminated | `False` |
| `background` | `kmerSetDB`/`None` | a `background` object containing _k_-mers (_k_=`Lmax`+1) which must be absent in the designed toolbox | `None` |
| `part_type` | `string` | must be either `'RNA'` or `'DNA'` depending on the type of genetic part being designed; ensures that correct folding free-energy parameters are used during structure evaluation | `'RNA'` |
| `struct_type` | `string` | must be either `'mfe'`, `'centroid'`, or `'both'` <br> `'mfe'` - use minimum free energy structure evaluation <br> `'centroid'` - use centroid structure evaluation <br> `'both'` - use both `'mfe'` and `'centroid'` evaluation | `'mfe'` |
| `seed` | `integer`/`None` | integer used to seed random number generations; two Maker runs with same constraints and seed value will generate the exact same toolbox; if `None` then a random seed value is used | `None` |
| `synth_opt` | `boolean` | if `True` then designed parts containing features that complicate DNA synthesis are eliminated | `False` |
| `local_model_fn` | `function`/`None` | a function with signature `'fn_name(seq)'` that takes in a partial genetic part sequence, and returns a tuple `(state, index)`, where `state` is either `True` or `False` and `index` is a traceback index / location or `None` depending on whether a custom design objective was met or not; useful for providing concurrent feedback to the path-finding process and steering nucleotide selection choices <br> **e.g.** `prevent_cutsites(seq)` maybe be a local function that takes in a partial sequence as it is built and returns `(True, None)` if the last six bases of the partial `seq[-6:]` is not the same as any of the cutsites used for cloning the part, else returns a tuple `(False, len(seq)-6)` as the traceback location to reselect the last six bases; this naturally ensures the final part is devoid of any cutsites used experimentally throughout the part | `None` |
| `global_model_fn` | `function`/`None` | a function with signature `'fn_name(seq)'` that takes in a complete genetic part sequence, and returns either `True` or `False` depending on whether a custom design objective was met; useful for design criteria that can only be evaluated when the complete genetic part is available; parts that are evaluated to be `False` are rejected and a new part generation is started; global model functions are evaluated only after the last base has been added to a genetic part under design <br> **e.g.** `gc_content(seq)` may be a global function that takes in a complete sequence and only accepts parts with GC content greater than threshold percentage (although, technically one can enforce this condition via a well planned local model function) | `None` |
| `jump_count` | `integer` | maximum number of restarts in path finding due to failure in finding suitable _k_-mers, meeting `local_model_fn`, or being stuck in local optima <br> (auto-adjusted with each iteration) | `10` |
| `fail_count` | `integer` | maximum number of consecutive failures tolerated when structure constraints, global model functions and synthesis objectives are not met <br> (auto-adjusted with each iteration) | `1000` |
| `output_file` | `string`/`None` | filename to store designed non-repetitive parts as they are generated consecutively; sequences are written in `FASTA` format | `None` |
| `verbose` | `boolean` | if `True` displays progress | `True` |

**_Returns_**: A `dictionary` of IUPAC strings with integer keys.

`Maker Mode` **API Examples**