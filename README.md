# README.py - `THAPBI-pycits`
This repository is for development of ITS1-based diagnostic/profiling tools for the THAPBI Phyto-Threats project, funded by BBSRC.

## Dependencies

### Python packages

Only required for the Python conversion, not the `Makefile`.

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