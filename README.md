![THAPBI-pycits TravisCI build status](https://api.travis-ci.org/widdowquinn/THAPBI-pycits.svg?branch=master) ![codecov.io coverage status](https://img.shields.io/codecov/c/github/widdowquinn/THAPBI-pycits.svg) [![Code Health](https://landscape.io/github/widdowquinn/THAPBI-pycits/master/landscape.svg?style=flat)](https://landscape.io/github/widdowquinn/THAPBI-pycits/master)

# README.py - `THAPBI-pycits`
This repository is for development of ITS1-based diagnostic/profiling tools for the THAPBI Phyto-Threats project, funded by BBSRC.

# DEVELOPER NOTES

## Python style conventions

In this repository, we're trying to keep to the Python [PEP8 style convention](https://www.python.org/dev/peps/pep-0008/), the [PEP257 docstring conventions](https://www.python.org/dev/peps/pep-0257/), and the [Zen of Python](https://www.python.org/dev/peps/pep-0020/). To help in this, a pre-commit hook script is provided in the `git_hooks` subdirectory that, if deployed in the `Git` repository, checks Python code for PEP8 correctness before permitting a `git commit` command to go to completion.

If the `pep8` module is not already present, it can be installed using `pip install pep8`

Whether you choose to use this or not, the `THAPBI-pycits` repository is registered with `landscape.io`, and the "health" of the code is assessed and reported for every repository push.

* [`landscape.io` page](https://landscape.io/github/widdowquinn/THAPBI-pycits)

### Installing the `git hook`

To install the pre-commit hook:

1. clone the repository with `git clone https://github.com/widdowquinn/THAPBI` (you may already have done this)
2. change directory to the root of the repository with `cd THAPBI-pycits`
3. copy the pre-commit script to the `.git/hooks` directory with `cp git_hooks/pre-commit .git/hooks/`

## Using a virtual environment with the repository

In the root directory of the repository:

```
$ virtualenv -p python3.5 venv-THAPBI-pycits
$ source venv-THAPBI-pycits/bin/activate
<activity>
$ deactivate
```

# INSTALLATION

## Dependencies: Python modules

All Python module dependencies are described in `requirements.txt` and can be installed using

```
pip install -r requirements.txt
```

There may be issues with `biom-format` and `biopython` installations due to ordering of module installation. If this is the case for you, then it might be solved by installing `numpy` at the command-line first, with:

```
pip install numpy
pip install -r requirements.txt
```

## Dependencies: Third-party applications


### `pear`

* [home page](http://sco.h-its.org/exelixis/web/software/pear/)

`pear` is a paired-end read merger, used by the pipeline to merge ITS paired-end reads into a single ITS sequence. It is available from the [`pear` home page](http://sco.h-its.org/exelixis/web/software/pear/) as a precompiled executable that can be placed in your `$PATH`, and it can be installed on the Mac with [Homebrew](http://brew.sh/) and [homebrew-science](https://github.com/Homebrew/homebrew-science), using:

```
brew install pear
```


### `Trimmomatic`

* [home page](http://www.usadellab.org/cms/?page=trimmomatic)

`Trimmomatic` is used to trim and quality-control the input reads. `pycits` expects `Trimmomatic` to be available at the command-line as `trimmomatic`. You can check if the tool is installed this way with the command:

```
which trimmomatic
```

To obtain `Trimmomatic` with this installation type on Linux systems, you can use:

```
apt-get install trimmomatic
```

and on the Mac (with [Homebrew](http://brew.sh/) and [homebrew-science](https://github.com/Homebrew/homebrew-science)):

```
brew install trimmomatic
```

If you have downloaded the Java `.jar.` file from [`trimmomatic`'s home page](http://www.usadellab.org/cms/?page=trimmomatic), you can wrap the `.jar` file with a Bash script called `trimmomatic` in your `$PATH`, such as

```
#!/bin/bash
exec java -jar $TRIMMOMATIC "$@"
```

where `$TRIMMOMATIC` is the path to your `trimmomatic .jar` file.




### More information

* Git hooks (`git`): [https://git-scm.com/book/en/v2/Customizing-Git-Git-Hooks](https://git-scm.com/book/en/v2/Customizing-Git-Git-Hooks)
* Git hooks (tutorial): [http://githooks.com/](http://githooks.com/)
* PEP8: [https://www.python.org/dev/peps/pep-0008/](https://www.python.org/dev/peps/pep-0008/)
* PEP257: [https://www.python.org/dev/peps/pep-0257/](https://www.python.org/dev/peps/pep-0257/)
* Zen of Python (PEP20): [https://www.python.org/dev/peps/pep-0020/](https://www.python.org/dev/peps/pep-0020/)


