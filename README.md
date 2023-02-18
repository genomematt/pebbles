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

Currently pebbles ignores quality scores and assumes all reads are a single unpaired unique observation.
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

# License
Pebbles is released under the BSD 3 Clause License https://opensource.org/license/bsd-3-clause/
