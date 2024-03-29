[![DOI](https://zenodo.org/badge/598950812.svg)](https://zenodo.org/badge/latestdoi/598950812)

# pebbles

Pebbles is a package for calling variants from single reads into HGVS format.

In many use cases we are interested in calling variants from a single read rather than the more
typical variant calling case where evidence from multiple reads is accumulated by genome position.
Examples of such uses include saturation genome editing and other multiplex assays of variant effect.

Pebbles is a lightweight converter of BAM or SAM files to genomic HGVS. This means that the resulting
calls will be in the format g.REFNAME:REFPOSITION and will not be corrected to be the most 3' variant 
(due to the way most aligners work they will be 5' shifted on the + reference strand).
For most uses you will want to process the output of pebbles with a HGVS validator and usually project
them into a different coordinate space such as the c. coordinates of the gene/transcript of interest.
This can be done with packages such as https://hgvs.readthedocs.io or websites like mutalyzer.nl

Currently pebbles ignores per base quality scores and assumes all reads are a single unpaired unique observation.
In most uses it will be desirable to align overlaps for paired end sequencing and trim reads before calling.


# Installation
To install from PyPI using pip:

```shell
pip install pebbles
```

To install from github using pip:

```shell
pip install git+https://github.com/genomematt/pebbles
```

Pebbles requires pysam and has only been tested with versions >= 0.20.0 (htslib 1.6)
Hatch is used as the build system, and will be required for source installs (and you need an up to date pip).

# Usage
For input pebbles requires a SAM or BAM file of alignment segments with MD tags. When using `minimap2` you
will need to map with the `--MD` argument.

Pebbles is in early and active development. Features and usage is likely to change as it is integrated into
other tools.

Pebbles can be used to call per read or to count occurrences of variants

To call per read
```shell
pebbles call myalignedsequences.bam > output.tsv
```

This will produce a tab seperated file with the read name, and a list of variants
identified in the read as a python list. For the test data of reads named with the expected variant in the tammar wallaby opsin
gene this output looks like:

```text
readname        call
WT      None
16_18delGAC     ['AY286018:g.16_18delGAC']
18_19insATG     ['AY286018:g.18_19insATG']
19_20delinsAG   ['AY286018:g.19_20delinsAG']
19_20delinsAG   ['AY286018:g.19_20delinsAG']
19_21delinsATG  ['AY286018:g.19_21delinsATG']
59A>T   ['AY286018:g.59A>T']
59A>T   ['AY286018:g.59A>T']
```

To generate counts
```shell
pebbles count myalignedsequences.bam > output.tsv
```

This will produce a tab seperated file with a column of variants and a column of counts.
For the test data this output looks like:

```text
variant      count
AY286018:g.16_18delGAC  1
AY286018:g.18_19insATG  1
AY286018:g.19_20delinsAG    2
AY286018:g.19_21delinsATG   1
AY286018:g.59A>T    2
```

For more detailed usage information see 
```shell
pebbles --help
```

# Usage as a CountESS plugin

The CountESS project is a graphical workflow manager for analysing count based datasets, in particular Deep Mutational
Scanning (DMS) and other Multiplex Assays of Variant Effects. CountESS is built with an entrypoint and inheritance
based plugin system, that pebbles implements.

To use pebbles in a CountESS workflow both Pebbles and CountESS need to be installed in the same python environment.
Once correctly installed the BAM and SAM parsing workflow steps should automatically be detected by CountESS and made
available in the user interface.

For further information on CountESS see 
[https://github.com/CountESS-Project/CountESS](https://github.com/CountESS-Project/CountESS)

# Contributing to Pebbles
Pebbles is licensed under the BSD-3-Clause license.  
You are free to fork this repository under the terms of that license.
If you have suggested changes please start by raising an issue in the issue tracker.
Pull requests are welcome and will be included at the discretion of the author.
Pull requests should be based on the 'develop' branch 
(with the exception of bugfixes where develop has diverged from main).
Bug reports should be made to the issue tracker.

Difficulty in understanding how to use the software is a documentation bug, and should also be raised on the
issue tracker so your question and my response are easily found by others.

Pebbles aims to maintain a respectful and inclusive community and adopts the
[contributor covenant v2.1](code_of_conduct.md)

# Citing Pebbles

Pebbles is currently unpublished. 
The current release can be cited using the Zenodo DOI. 
[![DOI](https://zenodo.org/badge/598950812.svg)](https://zenodo.org/badge/latestdoi/598950812)

# License
Pebbles is released under the BSD 3 Clause License https://opensource.org/license/bsd-3-clause/
