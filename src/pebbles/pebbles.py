# Pebbles. Sister of BammBamm Flinstone
# A per read bam mutation caller

import re
import pysam
import argparse
import sys
from collections import defaultdict

__version__ = "0.1.3"

def expand_cigar(cigar: str) -> str:
    """Convert a compact cigar string with state counts eg '2S10M1I5M1D1M'
    to a fully expanded string of state operators eg 'SSMMMMMMMMMMIMMMMMDM'
    Supports the operators MIDNSHPX=
    see https://en.wikipedia.org/wiki/Sequence_alignment#Representations

        Arguments:
            o   cigar : a cigar string
        Returns:
            o   string : a fully expanded string of operators
    """
    return "".join([match[1] * int(match[0]) for match in re.findall(r'([0-9]+)([MIDNSHPX=])', cigar)])


def engap(seq: str,
          cigar: str,
          is_reference: bool = False) -> str:
    """Convert a match/delete/insert string and sequence into gapped sequence
    To convert the target sequence swap delete and insert symbols.
        Arguments:
            o   seq : a sequence string
            o   cigar : a cigar string eg 80M5D10M10I
        Returns:
            o   string : gapped seqeunce
    """
    deletion_symbol = 'I' if is_reference else 'D'
    gapped = []
    xcigar = expand_cigar(cigar)
    seq = list(seq)
    for symbol in xcigar:
        if symbol == deletion_symbol:
            gapped.append('-')
        elif is_reference and symbol == 'S':
            gapped.append('-')
        else:
            gapped.append(seq.pop(0))
    return "".join(gapped)


def expand_mdtag(mdtag: str) -> str:
    """Convert a SAM MD tag eg '6^TAG3G2' containing information on
    the reference sequence content not included in the read to a string
    of expanded match tokens '.' and reference states
    eg '......TAG...G..'
    see https://samtools.github.io/hts-specs/SAMv1.pdf

        Arguments:
            o   mdtag : a SAM MD tag string
        Returns:
            o   string : a fully expanded string match and reference states
    """
    mdtag_tokens = re.findall(r'(\d+|\D+)', mdtag)
    result = []
    for token in mdtag_tokens:
        if token[0] in '1234567890':
            result.extend(['.', ] * int(token))
        elif token[0] == '^':
            result.extend(list(token[1:]))
        else:
            result.extend(list(token))
    return ''.join(result)


def call_mutations(refname: str,
                   pos: int,
                   expanded_engapped_md: str,
                   expanded_cigar: str,
                   gapped_read: str) -> list:
    """Call mutations in single SAM/BAM reads to HGVS format
    Assumes single end sequencing or merged paired end sequencing

    Arguments:
        refname: SAM reference sequence name
        pos: 1 based position in reference sequence
        expanded_engapped_md: expanded MD string with deletion states as gaps
        expanded_cigar: expanded string of cigar operators
        gapped_read: expanded read string with deletion states as gaps
    Returns:
        a list of HGVS formatted variant events
    """
    mutations = []
    i = 0
    nonref_bases = 0
    softmasked = 0
    while i < len(gapped_read):
        if expanded_cigar[i] == 'S':
            # softmasked state
            softmasked += 1
            i += 1
            continue
        if expanded_cigar[i] == 'M' and expanded_engapped_md[i] == '.' or \
                expanded_cigar[i] == '=':
            # match state with no variants
            i += 1
            continue
        if expanded_cigar[i] == 'D':
            deleted = ''
            delstart = i + pos + 1 + nonref_bases - softmasked
            while expanded_cigar[i] == 'D':
                deleted += expanded_engapped_md[i]
                i += 1
            mutations.append(f'{refname}:g.{delstart}_{i + pos + nonref_bases - softmasked}del{deleted}')
        if expanded_cigar[i] == 'I':
            inserted = ''
            insstart = i + pos + nonref_bases - softmasked  # base before first event base
            while expanded_cigar[i] == 'I':
                inserted += gapped_read[i]
                i += 1
                nonref_bases += 1
            mutations.append(f'{refname}:g.{insstart}_{insstart + 1}ins{inserted}')
        if (expanded_cigar[i] == 'M' and expanded_engapped_md[i] != '.') or \
                expanded_cigar[i] == 'X':
            mutant = ''
            reference = ''
            substart = i + pos + 1 + nonref_bases - softmasked
            while expanded_cigar[i] in 'MX' and expanded_engapped_md[i] != '.':
                mutant += gapped_read[i]
                reference += expanded_engapped_md[i]
                i += 1
            if len(mutant) == 1:
                mutations.append(f'{refname}:g.{substart}{reference}>{mutant}')
            else:
                mutations.append(f'{refname}:g.{substart}_{i + pos + nonref_bases - softmasked}delins{mutant}')
        i += 1
    return mutations


def call_mutations_from_pysam(pysamfile):
    """A generator function to call variants in all reads in a SAM/BAM
    file to HGVS format, on a per read basis.
    Assumes single end sequencing or merged paired end sequencing
    Arguments:
        pysamfile: a pysam.AlignmentFile object
                   eg SAM pysam.AlignmentFile("data.sam", "r")
                      BAM pysam.AlignmentFile("data.bam", "rb")

    Yields:
        a list of HGVS formatted variant events per read
    """

    for segment in pysamfile:
        if segment.is_qcfail or segment.is_supplementary or segment.is_unmapped:
            continue
        try:
            mdtag = segment.get_tag('MD')
        except KeyError:
            print(f'skipping {segment.qname} as it has no MD tag')
            # skip reads with no MD tag
            continue

        mutations = call_mutations(segment.reference_name,
                                   segment.pos,
                                   engap(expand_mdtag(mdtag),
                                         segment.cigarstring,
                                         is_reference=True),
                                   expand_cigar(segment.cigarstring),
                                   engap(segment.seq, segment.cigarstring))
        yield segment.qname, mutations


def fix_multi_variants(variants):
    if len(variants) > 1:
        return variants[0].split(':g.')[0] + ':g.[' + ";".join([variant.split(':g.')[-1] for variant in variants]) + ']'
    else:
        return variants[0]


def count(pysamfile, max=1):
    counts = defaultdict(int)
    for readname,variants in call_mutations_from_pysam(pysamfile):
        if variants and len(variants) <= max:
            counts[fix_multi_variants(variants)] += 1
    return ''.join([f'{key}\t{counts[key]}\n' for key in counts])


def cli(arguments: str = None):
    parser = argparse.ArgumentParser(description=f"pebbles v{__version__}" + """
                  (                                                                
             /((/ ###((%*                                                       
            /%(((((@(((%&(                                                      
             @(#(((((%                      @,     &#(#&&   *                   
                #%(((    (,              (*       @/@(**#&**                    
                  .*(#(#   /           @( (             %(///%(           @///% 
               / &(((((((((((@.        & @       &    &      (         &////(//%
               #((((#(((((((((((@      @(& %   @,@ /,,       /      &//@////(///
              @((((((@...&%((((((&       (& . %,,,,,,#       #   #//#//&//%/#//*
              (((((((&.....@@#(((#        ,&*,,,,,,,,,  @ ( .. &/////#///%/%/(  
              @(((#((@.....*.%&(((       .,,,##,(@,,,,/#/& # .///#//#//(//((.   
               ((..%(.........#.#/       /,,,,,,,,,,,,(* @ @///&////&//&/#      
                (*........*@...@#           #//.@(*,%&% %//(/////////@.         
                    *,....../*/.&        *,,,,,(,#*(@@&&%/#&*                   
            **@#%%(****%..#,             (,,,(,(,(###/*,,,,@                    
         @((/%******%/%***(                    #%(&##,,/*,,#                    
         &((((@*/@**(...%@.&               #/&,   &&@  (**#*                    
        %.%(((/((,**(...@...         ,%,,,,,(*/*%#(%  @@%/*.&                   
      /*.%......#*%@*&..* ..        *,,(,,,*,/&%,,*,&@//**%                     
         (*.&(**&     ,.,(., ,              (*,&,,,@   ./                       
    A program for calling variants from SAM and BAM files.
    Requires single end or merged reads.

    Created by Matthew Wakefield.
        Copyright (c) 2023  Matthew Wakefield. All rights reserved.

           This program is distributed in the hope that it will be useful,
           but WITHOUT ANY WARRANTY; without even the implied warranty of
           MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

    """, formatter_class=argparse.RawTextHelpFormatter)
    subparsers = parser.add_subparsers(dest='command')
    count = subparsers.add_parser('count',
                                  help='count occurrences of a variant')
    call = subparsers.add_parser('call',
                                  help='count occurrences of a variant')
    count.add_argument('--max',
                       type=int,
                       default=1,
                       help='Maximum number of variants in a read to include in count table'
                       )
    count.add_argument('infile',
                        type=str,
                        help='a SAM or BAM format file of mapped single end reads',
                        )
    call.add_argument('infile',
                        type=str,
                        help='a SAM or BAM format file of mapped single end reads',
                        )
    parser.add_argument('--version',
                        action='store_true',
                        help='print version information and exit')
    if arguments:
        args = parser.parse_args(arguments.split(' '))
    else:
        args = parser.parse_args()
    if args.version:
        print(f"pebbles v{__version__}")
        sys.exit()

    return parser.parse_args()


def main():
    args = cli()
    if args.infile.split('.')[-1].lower() == 'sam':
        infile = pysam.AlignmentFile(args.infile, "r")
    else:
        infile = pysam.AlignmentFile(args.infile, "rb")

    if args.command == 'call':
        print('readname\tvariants')
        for x in call_mutations_from_pysam(infile):
            print("\t".join([str(_) for _ in x]))
    elif args.command == 'count':
        print('variant\tcount')
        print(count(infile, max=args.max), end='')


if __name__ == '__main__':
    main()
