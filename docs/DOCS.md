<h1 align="center">
    <a href="https://github.com/ayaanhossain/nrpcalc/">
        <img src="https://raw.githubusercontent.com/ayaanhossain/nrpcalc/master/img/logo.svg?sanitize=true"  alt="Non-Repetitive Parts Calculator" width="418" class="center"/>
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

**background(path, Lmax, verbose=True)**

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
RuntimeError: kmerSetDB was dropped
```

## Finder Mode