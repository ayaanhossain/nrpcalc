

<h1 align="center">
    <a href="https://github.com/ayaanhossain/nrpcalc/">
        <img src="https://raw.githubusercontent.com/ayaanhossain/nrpcalc/master/img/logo.svg?sanitize=true"  alt="Non-Repetitive Parts Calculator" width="418" class="center"/>
    </a>
</h1>

<p align="center">
  <a href="#Overview">Overview</a> •
  <a href="#Installation">Installation</a> •
  <a href="#License">License</a> •
  <a href="#Citation">Citation</a> •
  <a href="#Documentation">Documentation</a> 
</p>

## Overview

Engineering large genetic circuits, pathways and whole genomes require large toolboxes of well characterized genetic parts. Existing toolboxes for composing such genetic systems from scratch are either small (forcing designers to re-use a particular part at multiple locations) and/or frequently share many repeats. Systems composed of repetitive elements are challenging to synthesize using commercial `DNA` synthesis approaches, often resulting in incomplete or mixed products. Moreover, when a repetitive genetic system introduced in a host organism causes biochemical stress, the cell may respond with homologous recombination – sections of the introduced system are deleted, effectively breaking it.  In _E. coli_, a `20`-bp repeat is enough  to trigger detectable levels of recombination.

To solve this _repeat challenge_ in synthetic biology, we developed and experimentally validated a suite of algorithms that automatically discovers and designs large toolboxes of highly non-repetitive genetic parts, collectively called the `Non-Repetitive Parts Calculator`. In `Finder Mode`, the algorithm discovers subsets of non-repetitive genetic parts from pre-existing toolboxes of characterized parts, while in `Maker Mode`, the algorithm designs large toolboxes of non-repetitive genetic parts for downstream characterization, based on specified design constraints such as a degenerate `DNA` or `RNA` sequence, an `RNA` secondary structure, and/or arbitrary model-based criteria.

Non-repetitiveness is a global property of the entire genetic part toolbox, and it is quantified by a repeat distribution (the number of repeats of length `L`) or more simply by the maximum shared repeat length, `Lmax`.  Both algorithms generate toolboxes according to a user-specified  `Lmax`  to control the desired level of non-repetitiveness. For example, when a toolbox of genetic parts is designed using an  `Lmax` of `10` base pairs, every genetic part in the toolbox can be simultaneously utilized in the same genetic system without introducing a repetitive sequence longer than `10` base pairs, which is necessary to ensure successful `DNA` synthesis and assembly.

<h3 align="center">
    <a href="https://github.com/ayaanhossain/nrpcalc/img/Fig1.svg">
        <img src="https://raw.githubusercontent.com/ayaanhossain/nrpcalc/master/img/Fig1.svg?sanitize=true"  alt="NRP Calculator Maker Mode" width="1000" class="center"/>
    </a>
</h3>

## Installation

`Non-Repetitive Parts Calculator` is a `Linux/MacOS`-tested software, and built with `Python 2.7`. The software is not `Python 3.x` compatible at the moment, but will be pretty soon.

The best way to install `Non-Repetitive Parts Calculator` is via `conda`. If you have either `anaconda` or `miniconda` installed on your system, you are good to proceed. Otherwise, you may first need to install [miniconda](https://docs.conda.io/en/latest/miniconda.html).

Assuming you have `conda` available on your system, you may first create a new environment (here we are naming it `nrpcalc` but you can name it as you like)
```bash
$ conda create -n nrpcalc python=2.7
```
Once your new environment is ready, deactivate your current environment
```bash
$ conda deactivate
```
and activate the newly created `nrpcalc` environment
```bash
$ conda activate nrpcalc
```

**Note** You don't necessarily need to make a new environment, as long as you know which environment you want `NRP Calculator` installed in, and have it activated during installation. For example, you might have a different project environment inside which you want to install `NRP Calculator`.

The first thing we will need to install in a new environment is `ViennaRNA` which is an external dependency for `NRP Calculator` and not available in [PyPI](https://pypi.org). `ViennaRNA` is easily installed with

```bash
$ conda install -c bioconda viennarna
```
If you are following the instructions as is, this will install the latest copy of `ViennaRNA` inside the new `nrpcalc` environment.

If you do not want to install `conda` on your system, you will need to install `ViennaRNA` [manually](https://www.tbi.univie.ac.at/RNA/ViennaRNA/doc/html/install.html), ensuring proper setup of `PYTHONPATH` etc.

**Note** If you are using `conda` then the following step will make `NRP Calculator` available only inside the environment in which it is being installed. If you want to use software, you will always need to activate that specific environment, for example, here the environment is named `nrpcalc`.

Once `ViennaRNA` is installed, you can easily install `NRP Calculator` via
```bash
$ pip install --upgrade nrpcalc
```
which will additionally install all PyPI dependencies.

If everything went well, `NRP Calculator` is now available in your environment under the `nrpcalc` package. You may verify it like so:
```bash
$ python
Python 2.7.15 | packaged by conda-forge | (default, Mar  5 2020, 14:56:06) 
[GCC 7.3.0] on linux2
Type "help", "copyright", "credits" or "license" for more information.
>>> 
>>> import nrpcalc
>>> 
>>> print nrpcalc.__doc__

Non-Repetitive Parts Calculator

Automated design and discovery of non-repetitive genetic
parts for engineering stable systems.

Version: 1.1.3

Authors: Ayaan Hossain <auh57@psu.edu>
         Howard Salis  <salis@psu.edu>

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

>>> 
```

**Uninstallation** of `Non-Repetitive Parts Calculator` is as easy as
```bash
$ pip uninstall nrpcalc
```

## License

`Non-Repetitive Parts Calculator` (c) 2020 Ayaan Hossain.

`Non-Repetitive Parts Calculator` is an **open-source software** under [MIT](https://opensource.org/licenses/MIT) License.

See [LICENSE](./LICENSE) file for more details.

## Citation

If you use `Non-Repetitive Parts Calculator` or toolboxes designed or discovered using the algorithm in your research publication, please cite

```
Hossain A., Lopez E., Halper S.M., Cetnar D.P., Reis A.C., Strickland D., Klavins E., and Salis H.M.
Automated design of thousands of nonrepetitive parts for engineering stable genetic systems
Nature Biotechnology, doi:10.1038/s41587-020-0584-2
```

You can read the complete article online at [Nature Biotechnology](https://www.nature.com/articles/s41587-020-0584-2).

**Abstract** Engineered genetic systems are prone to failure when their genetic parts contain repetitive sequences. Designing many non-repetitive  genetic  parts  with  desired  functionalities  remains  a  significant  challenge  with  high  computational  complexity.  To  overcome this challenge, we developed the nonrepetitive parts calculator to rapidly generate thousands of highly nonrepetitive genetic  parts  from  specified  design  constraints,  including  promoters,  ribosome-binding  sites  and  terminators.  As  a  demonstration, we designed and experimentally characterized 4,350 nonrepetitive bacterial promoters with transcription rates that varied across a 820,000-fold range, and 1,722 highly nonrepetitive yeast promoters with transcription rates that varied across a  25,000-fold  range.  We  applied  machine  learning  to  explain  how  specific  interactions  controlled  the  promoters’  transcription rates. We also show that using nonrepetitive genetic parts substantially reduces homologous recombination, resulting in greater genetic stability.

**Acknowledgements** This project was supported by funds from the Air Force Office of Scientific Research (grant no. FA9550-14-1-0089), the Defense Advanced Research Projects Agency (grant nos. FA8750-17-C-0254 and HR001117C0095), the Department of Energy (grant no. DE-SC0019090), and a Graduate Research Innovation award from the Huck Institutes of the Life Sciences.

**Contributions** A.H. and H.M.S. conceived the study. A.H., E.L., D.P.C., S.M.H., A.C.R. and D.S. designed and carried out the experiments. A.H., A.C.R. and H.M.S. developed the algorithms and performed the data analysis. A.H., D.S., E.K. and H.M.S. wrote the manuscript.

**Maintenance** `Non-Repetitive Parts Calculator` is currently maintained by
- Ayaan Hossain | [github.com/ayaanhossain](https://github.com/ayaanhossain) | [@bioalgorithmist](https://twitter.com/bioalgorithmist)

## Documentation
A write-up explaining the `nrpcalc` API as well as a toy example demonstrating `NRP Calculator` is in the works, and will be released by Friday, June 17, 2020. For the moment, please access the API documentation from within the `Python` REPL via `print nrpcalc.finder.__doc__`, `print nrpcalc.maker.__doc__` and `print nrpcalc.background.__doc__`. Meanwhile, please enjoy the paper!