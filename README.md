<h1 align="center">
    <a href="https://github.com/ayaanhossain/nrpcalc/">
        <img src="https://raw.githubusercontent.com/ayaanhossain/nrpcalc/master/img/logo.svg?sanitize=true"  alt="Non-Repetitive Parts Calculator" width="460" class="center"/>
    </a>
</h1>

<p align="center">
  <a href="#Overview">Overview</a> •
  <a href="#Installation">Installation</a> •
  <a href="#License">License</a> •
  <a href="#Citation">Citation</a> •
  <a href="https://github.com/ayaanhossain/nrpcalc/blob/master/docs/DOCS.md">DOCS</a> •
  <a href="#NRP-Calculator-in-Action">NRP Calculator in Action!</a>
</p>

## Overview

Engineering large genetic circuits, pathways and whole genomes require large toolboxes of well characterized genetic parts. Existing toolboxes for composing such genetic systems from scratch are either small (forcing designers to re-use a particular part at multiple locations) and/or frequently share many repeats. Systems composed of repetitive elements are challenging to synthesize using commercial DNA synthesis approaches, often resulting in incomplete or mixed products. Moreover, when a repetitive genetic system introduced in a host organism causes biochemical stress, the cell may respond with homologous recombination – sections of the introduced system are deleted, effectively breaking it.  In _E. coli_, a 20-bp repeat is enough  to trigger detectable levels of recombination.

To solve this _repeat challenge_ in synthetic biology, we developed and experimentally validated a suite of algorithms that automatically discovers and designs large toolboxes of highly non-repetitive genetic parts, collectively called the `Non-Repetitive Parts Calculator`. In `Finder Mode`, the algorithm discovers subsets of non-repetitive genetic parts from pre-existing toolboxes of characterized parts, while in `Maker Mode`, the algorithm designs large toolboxes of non-repetitive genetic parts for downstream characterization, based on specified design constraints such as a degenerate DNA or RNA sequence, an RNA secondary structure, and/or arbitrary model-based criteria.

Non-repetitiveness is a global property of the entire genetic part toolbox, and it is quantified by a repeat distribution (the number of repeats of length `L`) or more simply by the maximum shared repeat length, `Lmax`.  Both algorithms generate toolboxes according to a user-specified  `Lmax`  to control the desired level of non-repetitiveness. For example, when a toolbox of genetic parts is designed using an  `Lmax` of 10 base pairs, every genetic part in the toolbox can be simultaneously utilized in the same genetic system without introducing a repetitive sequence longer than 10 base pairs, which is necessary to ensure successful DNA synthesis and assembly.

<h3 align="center">
    <a href="https://github.com/ayaanhossain/nrpcalc/img/Fig1.svg">
        <img src="https://raw.githubusercontent.com/ayaanhossain/nrpcalc/master/img/Fig1.svg?sanitize=true"  alt="NRP Calculator Maker Mode" width="800" class="center"/>
    </a>
</h3>

## Installation

`Non-Repetitive Parts Calculator` is a `Linux`/`MacOS`-tested software, and was originally built with `Python2.7`. The software is now `Python3` exclusive, and compatible only with `Python3.6` and above. `Python2.7` is no longer supported.


**Setting up the Environment**

The best way to install `NRP Calculator` is via `conda`. If you have either `anaconda` or `miniconda` installed on your system, you are good to proceed. Otherwise, you may first need to install [miniconda](https://docs.conda.io/en/latest/miniconda.html).

Assuming you have `conda` available on your system, you may first create a new environment. Here we are naming it `nrpenv` but you can name it as you like. Also, we are setting `nrpenv` up with `Python3.6` but feel free to install `Python3.7`, or one of the newer ones.
```bash
$ conda create -n nrpenv python=3.6
```
Once your new environment is ready, deactivate your current environment
```bash
$ conda deactivate
```
and activate the newly created `nrpenv` environment
```bash
$ conda activate nrpenv
```

> **Note** You don't necessarily need to make a new environment, as long as you know which environment you want `NRP Calculator` installed in, and have it activated during installation. For example, you might have a different project environment inside which you want to install `NRP Calculator`.


**Approach 1: Install from PyPI**

You can install `NRP Calculator` from PyPI, where it is published as the `nrpcalc` package. This is as easy as
```bash
$ pip install --upgrade nrpcalc --no-cache-dir
```
which will also install all dependencies from PyPI.


**Approach 2: Install from GitHub**

Alternatively, you can install `NRP Calculator` from GitHub. To do so, first clone the repository with
```bash
$ git clone https://github.com/ayaanhossain/nrpcalc.git
```
Once downloaded, navigate into the `nrpcalc` directory with
```bash
$ cd nrpcalc
```
and install `NRP Calculator` using the included `setup.py`
```bash
$ python setup.py install
```

**Verifying Installation**

If everything went well, `NRP Calculator` is now available in your environment under the `nrpcalc` package. You may verify it like so:
```python
$ python
Python 3.6.10 | packaged by conda-forge | (default, Apr 24 2020, 16:44:11)
[GCC 7.3.0] on linux
Type "help", "copyright", "credits" or "license" for more information.
>>>
>>> import nrpcalc
>>>
>>> print(nrpcalc.__doc__)

Non-Repetitive Parts Calculator

Automated design and discovery of non-repetitive genetic
parts for engineering stable systems.

Version: 1.6.3

Authors: Ayaan Hossain <auh57@psu.edu>
         Howard Salis  <salis@psu.edu>

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

>>>
```

> **Note** Remember to activate the specific environment in which `NRP Calculator` is installed in order to use it. For example, if you followed the instructions above exactly, then the environment you created and installed `NRP Calculator` in is named `nrpenv`.

**Uninstalling `NRP Calculator`**

You can easily remove `Non-Repetitive Parts Calculator` with
```bash
$ pip uninstall nrpcalc
```

**Reporting Installation Issues**

If you encounter any problems during installation, please feel free to [open an issue](https://github.com/ayaanhossain/nrpcalc/issues) describing your problem along with your OS details, and any console output that shows the error.

## License

`Non-Repetitive Parts Calculator` (c) 2020 Ayaan Hossain.

`Non-Repetitive Parts Calculator` is an **open-source software** under [MIT](https://opensource.org/licenses/MIT) License.

See [LICENSE](https://github.com/ayaanhossain/nrpcalc/blob/master/LICENSE) file for more details.

## Citation

If you use `Non-Repetitive Parts Calculator` or toolboxes designed or discovered using the algorithm in your research publication, please cite

```
Hossain A., Lopez E., Halper S.M., Cetnar D.P., Reis A.C., Strickland D., Klavins E., and Salis H.M.
Automated design of thousands of nonrepetitive parts for engineering stable genetic systems
Nature Biotechnology, doi:10.1038/s41587-020-0584-2
```

You can read the complete article online at [Nature Biotechnology](https://www.nature.com/articles/s41587-020-0584-2). A free PDF copy of the paper is accessible [here](https://rdcu.be/b5Aw0).

**Abstract** Engineered genetic systems are prone to failure when their genetic parts contain repetitive sequences. Designing many non-repetitive  genetic  parts  with  desired  functionalities  remains  a  significant  challenge  with  high  computational  complexity.  To  overcome this challenge, we developed the nonrepetitive parts calculator to rapidly generate thousands of highly nonrepetitive genetic  parts  from  specified  design  constraints,  including  promoters,  ribosome-binding  sites  and  terminators.  As  a  demonstration, we designed and experimentally characterized 4,350 nonrepetitive bacterial promoters with transcription rates that varied across a 820,000-fold range, and 1,722 highly nonrepetitive yeast promoters with transcription rates that varied across a  25,000-fold  range.  We  applied  machine  learning  to  explain  how  specific  interactions  controlled  the  promoters’  transcription rates. We also show that using nonrepetitive genetic parts substantially reduces homologous recombination, resulting in greater genetic stability.

**Acknowledgements** This project was supported by funds from the Air Force Office of Scientific Research (grant no. FA9550-14-1-0089), the Defense Advanced Research Projects Agency (grant nos. FA8750-17-C-0254 and HR001117C0095), the Department of Energy (grant no. DE-SC0019090), and a Graduate Research Innovation award from the Huck Institutes of the Life Sciences.

**Contributions** A.H. and H.M.S. conceived the study. A.H., E.L., D.P.C., S.M.H., A.C.R. and D.S. designed and carried out the experiments. A.H., A.C.R. and H.M.S. developed the algorithms and performed the data analysis. A.H., D.S., E.K. and H.M.S. wrote the manuscript.

**Maintenance** `Non-Repetitive Parts Calculator` is currently maintained by
- Ayaan Hossain | [github.com/ayaanhossain](https://github.com/ayaanhossain) | [@bioalgorithmist](https://twitter.com/bioalgorithmist)
- Howard Salis | [github.com/hsalis](https://github.com/hsalis) | [@hsalis](https://twitter.com/hsalis)

## API Documentation

A detailed description of `nrpcalc` `Maker Mode`, `Finder Mode`, and `Background` API, along with example use cases can be found in [DOCS.md](https://github.com/ayaanhossain/nrpcalc/blob/master/docs/DOCS.md). You may also check the API documentation from within the `Python` REPL via
* `print(nrpcalc.background.__doc__)`,
* `print(nrpcalc.finder.__doc__)`, and
* `print(nrpcalc.maker.__doc__)`.

If you enconter any problem using the API, please feel free to [open an issue](https://github.com/ayaanhossain/nrpcalc/issues) describing your use case, along with minimal code snippets reproducing the problem and any console output that shows the problem or error.

## `NRP Calculator` in Action!

A `jupyter` [notebook](https://github.com/ayaanhossain/nrpcalc/blob/master/examples/NRPCalcInAction.ipynb) describing the use of `NRP Calculator` in designing thousands of commonly used genetic parts from constitutive promoters to intrinsic terminators for _E. coli_ systems engineering is available inside `/examples/` directory.

If you have [installed](#Installation) `nrpcalc` then you have everything you need to open and execute the notebook.

Simply clone this repository
```bash
$ git clone https://github.com/ayaanhossain/nrpcalc.git
```

and navigate to `nrpcalc/examples/` directory
```bash
$ cd nrpcalc/examples/
```

and fire up `jupyter`.
```bash
$ jupyter notebook
```

This should now open a browser tab showing a single `jupyter` notebook named `NRPCalcInAction.ipynb`. Click on the notebook to open it, and follow the content of the notebook.

If any part of the notebook is unclear or confusing to you, please feel free to reach any of the authors via Email or Twitter, who would be more than happy to clarify anything in the notebook.

We hope you enjoy the paper and the notebook!
