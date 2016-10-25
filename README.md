![THAPBI-pycits TravisCI build status](https://api.travis-ci.org/widdowquinn/THAPBI-pycits.svg?branch=master) ![codecov.io coverage status](https://img.shields.io/codecov/c/github/widdowquinn/THAPBI-pycits.svg)

* [Travis-CI page](https://travis-ci.org/widdowquinn/THAPBI-pycits/branches)
* [`codecov`.io page](https://codecov.io/gh/widdowquinn/THAPBI-pycits)

# README.py - `THAPBI-pycits`
This repository is for development of ITS1-based diagnostic/profiling tools for the THAPBI Phyto-Threats project, funded by BBSRC.

## Using a virtual environment

In the root directory of the repository:

```
$ virtualenv -p python3.5 venv-THAPBI-pycits
$ source venv-THAPBI-pycits/bin/activate
<activity>
$ deactivate
```

## Dependencies

### Python packages

* `biom-format`: available through `PyPI` at [https://pypi.python.org/pypi/biom-format](https://pypi.python.org/pypi/biom-format) (tested with v2.1.5)

### Third-party tools

Required for both the Python and `Makefile` conversions.

These tools should either be in the `${PATH}`, or their executables should be specified in the script options.

* `FastQC`: tested with v0.11.3
* `seq_crumbs`: tested with v0.1.9 [https://github.com/JoseBlanca/seq_crumbs](https://github.com/JoseBlanca/seq_crumbs)
* `BLAST+`: tested with v2.2.31+
* `QIIME`: tested with v1.9.1
* `MUSCLE`: tested with v3.8.31


## Python conversion

The initial Python pipeline is a `Python` transliteration of Santi's script, with supporting modules that wrap the packages to be called by the pipeline, and provide helper functions. The script and modules can be installed to the system with:

```
$ python setup.py install
```

which makes the script available at the command-line:

```
$ thapbi_santi_otus.py -h
usage: thapbi_santi_otus.py [-h] [-p PREFIX] [-i INDIRNAME] [-o OUTDIRNAME]
                            [-r REFERENCE_FASTA] [-v] [-l LOGFILE]
                            [-t THREADS] [--fastqc FASTQC]
                            [--trim_quality TRIM_QUALITY]
                            [--join_paired_ends JOIN_PAIRED_ENDS]
                            [--convert_format CONVERT_FORMAT]
                            [--blastclust BLASTCLUST] [--muscle MUSCLE]
                            [--pick_otus PICK_OTUS]
                            [--pick_closed_reference_otus PICK_CLOSED_REFERENCE_OTUS]
optional arguments:
  -h, --help            show this help message and exit
  -p PREFIX, --prefix PREFIX
                        Paired readfiles prefix
  -i INDIRNAME, --indir INDIRNAME
                        Path to directory containing input reads
  -o OUTDIRNAME, --outdir OUTDIRNAME
                        Path to directory to write output
  -r REFERENCE_FASTA, --reference REFERENCE_FASTA
                        Path to reference sequence FASTA file
  -v, --verbose         Report verbose output
  -l LOGFILE, --logfile LOGFILE
                        Logfile location
  -t THREADS, --threads THREADS
                        Number of threads to use (default: all)
  --fastqc FASTQC       Path to FastQC executable
  --trim_quality TRIM_QUALITY
                        Path to seq_crumbs trim_quality script
  --join_paired_ends JOIN_PAIRED_ENDS
                        Path to ea-utils join_paired_ends.py script
  --convert_format CONVERT_FORMAT
                        Path to seq_crumbs convert_format script
  --blastclust BLASTCLUST
                        Path to blastclust
  --muscle MUSCLE       Path to MUSCLE
  --pick_otus PICK_OTUS
                        Path to QIIME pick_otus.py script
  --pick_closed_reference_otus PICK_CLOSED_REFERENCE_OTUS
                        Path to QIIME pick_closed_reference_otus.py script
```

Modules are found under the `thapbi_santi` directory:

```
$ tree thapbi_santi
thapbi_santi
├── __init__.py
├── blast.py
├── ea_utils.py
├── fastqc.py
├── muscle.py
├── qiime.py
├── seq_crumbs.py
└── tools.py
```

The advantages of the Python script are that it can be installed once on the system, and is available everywhere; a comprehensive logfile is produced, which facilitates reproducibility of all commands that were run; input file prefix and location, and output directory can be specified directly, facilitating distribution of the script across several nodes; and output from each package is segregated into its own directory.

## Python style conventions

In this project, we're trying to keep to the Python [PEP8 style convention](https://www.python.org/dev/peps/pep-0008/), the [PEP257 docstring conventions](https://www.python.org/dev/peps/pep-0257/), and the [Zen of Python](https://www.python.org/dev/peps/pep-0020/). To help in this, a pre-commit hook script is provided in the `git_hooks` subdirectory that, if deployed in the `Git` repository, checks Python code for PEP8 correctness before permitting a `git commit` command to go to completion.

If the `pep8` module is not already present, it can be installed using `pip install pep8`

### Installing the `git hook`

To install the pre-commit hook:

1. clone the repository with `git clone https://github.com/widdowquinn/THAPBI` (you may already have done this)
2. change directory to the root of the repository with `cd THAPBI-pycits`
3. copy the pre-commit script to the `.git/hooks` directory with `cp git_hooks/pre-commit .git/hooks/`

### More information

* Git hooks (`git`): [https://git-scm.com/book/en/v2/Customizing-Git-Git-Hooks](https://git-scm.com/book/en/v2/Customizing-Git-Git-Hooks)
* Git hooks (tutorial): [http://githooks.com/](http://githooks.com/)
* PEP8: [https://www.python.org/dev/peps/pep-0008/](https://www.python.org/dev/peps/pep-0008/)
* PEP257: [https://www.python.org/dev/peps/pep-0257/](https://www.python.org/dev/peps/pep-0257/)
* Zen of Python (PEP20): [https://www.python.org/dev/peps/pep-0020/](https://www.python.org/dev/peps/pep-0020/)


