# README.py  - `santi_script`

This directory contains conversions of Santiago Garcia's OTU-classification script from the original Word document reference to Bash, Makefile and Python.

## Original pipeline

The original pipeline information and Santi's Python scripts are in the subdirectory `santi_docs`:

```
$ tree santi_docs/
santi_docs/
├── Santi pipeline.docx
├── blastclust_lst2fasta.py
├── database.fasta
└── trim_longitudes.py
```

* `Santi pipeline.docx`: The source material for reconstructing the pipeline.
* `database.fasta`: A dataset of reference *Phytophthora* ITS sequences
* `blastclust_lst2fasta.py`: Python script to convert `BLASTCLUST` output to a set of `FASTA` files, one file per cluster.
* `trim_longitudes.py`: Python script to trim joined reads by a fixed number of bases, from each end.

## Shell script conversion

The `bash` shell script conversion is a transliteration of the pipeline given in the file `santi_docs/Santi pipeline.docx`, with comments to outline the actions at each step. This script is not very flexible, and does not segregate input from output from intermediate files very clearly. The script must be run in the same directory as the input read files.

* `santi_script.sh`: Transliteration of the original pipeline into `bash` shell script.

**Usage:**

```
$ santi_script.sh
```

## Makefile conversion

The `Make` conversion introduces the ability to pass the input file prefix and sequence directory as arguments, so no longer requires to be run in the same directory as the input file. Pipeline dependencies are recognised, which saves time on re-running if the pipeline fails. However, the `Makefile` remains relatively inflexible, in that all commands/packages must be in the `${PATH}`, and the `Makefile` must still be copied to, and run in, the intended output directory (though a small modification to the script could avoid this).

* `Makefile`: Conversion of original pipeline to `Makefile`

**Usage:**

```
$ make otus -e PREFIX=<prefix of sequence data> SEQDIR=<path to sequence directory>
```

## Python conversion

The pipeline is also converted to a `Python` script, with supporting modules that wrap the packages to be called by the pipeline, and provide helper functions. The script and modules can be installed to the system with:

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