# Pebbles. Sister of BammBamm Flinstone
# A per read bam mutation caller
import argparse
import collections
import re
import shlex
import sys
import traceback
from collections import defaultdict
from itertools import islice
from collections.abc import Iterable, Mapping, Sequence
from typing import Tuple, Optional, Union
from pebbles import VERSION

import pysam


class Logger:
    """Basic Logger compatible with CountESS.core.logger"""

    def __init__(self, stdout=sys.stdout, stderr=sys.stderr, prefix: Optional[str] = None):
        self.stdout = stdout
        self.stderr = stderr
        self.prefix = prefix

    def progress(self, message: str = "Running", percentage: Optional[int] = None):
        if self.prefix:
            message = self.prefix + ": " + message
        if percentage:
            message += f" [{int(percentage):2d}%]"
        self.stdout.write(f"{message}\n")

    def log(self, level: str, message: str, detail: Optional[str] = None):
        if self.prefix:
            message = self.prefix + ": " + message
        if detail:
            message += " " + repr(detail)

        self.stderr.write(message + "\n")

    def info(self, message: str, detail: Optional[str] = None):
        """Log a message at level info"""
        self.log("info", message, detail)

    def warning(self, message: str, detail: Optional[str] = None):
        """Log a message at level warning"""
        self.log("warning", message, detail)

    def error(self, message: str, detail: Optional[str] = None):
        """Log a message at level error"""
        self.log("error", message, detail)

    def exception(self, exception: Exception):
        self.error(str(exception), detail="".join(traceback.format_exception(exception)))

    def clear(self):
        """Clear logs (if possible)"""
        return None


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
    with one character per expanded cigar operator.
    For a read sequence, bases are consumed by M I S = X operators.
    For the reference (an expanded MD tag, is_reference=True), states are
    consumed by M D = X operators, as these are the only ones the MD tag
    describes. All other operators (eg D, N, H, P for reads) are gaps.
        Arguments:
            o   seq : a sequence string
            o   cigar : a cigar string eg 80M5D10M10I
            o   is_reference : seq is an expanded MD tag rather than a read
        Returns:
            o   string : gapped seqeunce
    """
    consuming = 'MD=X' if is_reference else 'MIS=X'
    gapped = []
    xcigar = expand_cigar(cigar)
    seq_list = list(seq)
    for symbol in xcigar:
        if symbol in consuming:
            gapped.append(seq_list.pop(0))
        else:
            gapped.append('-')
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
    nonref_bases = 0  # inserted read bases so far; these columns do not consume reference
    softmasked = 0
    while i < len(gapped_read):
        if expanded_cigar[i] in 'SHP':
            # soft clipped, hard clipped or padding state: no reference consumed
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
            delstart = i + pos + 1 - nonref_bases - softmasked
            while i < len(expanded_cigar) and expanded_cigar[i] == 'D':
                deleted += expanded_engapped_md[i]
                i += 1
            mutations.append(f'{refname}:g.{delstart}_{i + pos - nonref_bases - softmasked}del{deleted}')
            # re-evaluate the state at i, which may be an adjacent event
            continue
        if expanded_cigar[i] == 'I':
            inserted = ''
            insstart = i + pos - nonref_bases - softmasked  # base before first event base
            while i < len(expanded_cigar) and expanded_cigar[i] == 'I':
                inserted += gapped_read[i]
                i += 1
                nonref_bases += 1
            mutations.append(f'{refname}:g.{insstart}_{insstart + 1}ins{inserted}')
            continue
        if (expanded_cigar[i] == 'M' and expanded_engapped_md[i] != '.') or \
                expanded_cigar[i] == 'X':
            mutant = ''
            reference = ''
            substart = i + pos + 1 - nonref_bases - softmasked
            while len(expanded_cigar) > i and expanded_cigar[i] in 'MX' and len(expanded_engapped_md) > i and expanded_engapped_md[i] != '.':
                mutant += gapped_read[i]
                reference += expanded_engapped_md[i]
                i += 1
            if not mutant:
                # X operator without a matching MD mismatch: nothing to call
                i += 1
                continue
            if len(mutant) == 1:
                mutations.append(f'{refname}:g.{substart}{reference}>{mutant}')
            else:
                mutations.append(f'{refname}:g.{substart}_{i + pos - nonref_bases - softmasked}delins{mutant}')
            continue
        # any other operator (eg N, P) does not produce a call
        i += 1
    return mutations


def _warn(message: str, logger: Optional[Logger] = None) -> None:
    """Report a warning via logger if given, otherwise on stderr so that
    it does not corrupt tabular output written to stdout"""
    if logger:
        logger.warning(message)
    else:
        print(message, file=sys.stderr)


def call_mutations_from_pysam(pysamfile: collections.abc.Iterable,
                              min_quality: int =0,
                              logger: Optional[Logger] = None,
                              ) -> Iterable[Tuple[str, list]]:
    """A generator function to call variants in all reads in a SAM/BAM
    file to HGVS format, on a per read basis.
    Assumes single end sequencing or merged paired end sequencing
    qcfail, secondary, supplementary and unmapped segments are ignored,
    as are segments with no stored sequence (SEQ of *)

    Arguments:
        pysamfile: a pysam.AlignmentFile object
                   eg SAM pysam.AlignmentFile("data.sam", "r")
                      BAM pysam.AlignmentFile("data.bam", "rb")
        min_quality:  the minimum mapping quality score for a reported allele. Default 0.
        logger:    a logger object with a .warning() method such as CountESS.core.logger
                   default: None


    Yields:
        a list of HGVS formatted variant events per read
    """

    for segment in pysamfile:
        if segment.is_qcfail or segment.is_secondary or segment.is_supplementary or segment.is_unmapped:
            continue
        if segment.mapping_quality < min_quality:
            continue
        if segment.query_sequence is None:
            _warn(f'skipping {segment.qname} as it has no sequence', logger)
            continue
        try:
            mdtag = segment.get_tag('MD')
        except KeyError:
            _warn(f'skipping {segment.qname} as it has no MD tag', logger)
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


def fix_multi_variants(variants: list) -> str:
    """Convert a list of multiple HGVS genomic variants to a correctly specified allele
    eg ['AY286018:g.16_18delGAC', 'AY286018:g.59A>T'] ->  'AY286018:g.[16_18delGAC;59A>T]'
    Assumes all variants are allelic on the same reference
    """
    if len(variants) > 1:
        return variants[0].split(':g.')[0] + ':g.[' + ";".join([variant.split(':g.')[-1] for variant in variants]) + ']'
    else:
        return variants[0]


def count_dict(pysamfile: Iterable,
               max_variants: int = 1,
               row_limit: Optional[int] = None,
               min_quality: int = 0,
               logger: Optional[Logger] = None,
               ) -> Mapping[str, int]:
    """
    Counts occurrences of alleles in a SAM or BAM file

    Arguments:
        pysamfile    : an iterable of alignment segment objects from pysam
        max_variants : the maximum number of variants in a reported allele. Default 1
        row_limit    : the maximum number of alignments to process. Default None (ie process all)
        min_quality  : the minimum mapping quality score for a reported allele. Default 0
        logger       : a logger object with a .progress() & .warning() method such as CountESS.core.logger
                       default: None


    Returns:
        A dictionary keyed by allele HGVS strings of counts
    """
    counts : dict[str,int] = defaultdict(int)
    number = 0
    for readname, variants in islice(call_mutations_from_pysam(pysamfile, min_quality, logger), row_limit):
        number += 1
        if logger and number % 1000 == 0:
            logger.progress(f"Loading {pysamfile.filename.decode()}")
        if variants and len(variants) <= max_variants:
            counts[fix_multi_variants(variants)] += 1

    if logger:
        logger.progress(f"Loaded {pysamfile.filename.decode()}",100)
        if len(counts) == 0 and row_limit is not None:
            logger.warning(f"No mutations found in first {row_limit} alignments")

    return counts


def count(pysamfile: Iterable,
          max_variants: int = 1,
          min_quality: int = 0,
          ) -> str:
    """
    Counts occurrences of alleles in a SAM or BAM file

    Arguments:
        pysamfile - an iterable of alignment segment objects from pysam
        max_variants - the maximum number of variants in a reported allele. Default 1
        min_quality - the minimum mapping quality score for a reported allele. Default 0

    Returns:
        A tsv formatted text string of allele HGVS description and counts
    """
    counts = count_dict(pysamfile, max_variants=max_variants, min_quality=min_quality)
    return ''.join([f'{key}\t{counts[key]}\n' for key in counts])


def cli(arguments: Optional[Union[str, Sequence[str]]] = None) -> argparse.Namespace:
    """command line interface. Use pebbles -h to see help

    Arguments:
        arguments: a command line string (eg 'count --max 2 in.bam') or a list
                   of arguments. Default None uses sys.argv
    """
    parser = argparse.ArgumentParser(description=f"pebbles v{VERSION}" + """
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
    parser.add_argument('--version',
                        action='version',
                        version=f"pebbles v{VERSION}",
                        help='print version information and exit')
    subparsers = parser.add_subparsers(dest='command', required=True, metavar='{count,call}')
    count = subparsers.add_parser('count',
                                  help='count occurrences of a variant')
    call = subparsers.add_parser('call',
                                 help='count occurrences of a variant')
    count.add_argument('--max',
                       type=int,
                       default=1,
                       help='Maximum number of variants in a read to include in count table'
                       )
    count.add_argument('--min_quality',
                       type=int,
                       default=0,
                       help='Minimum quality score required to count a variant'
                       )
    count.add_argument('infile',
                       type=str,
                       help='a SAM or BAM format file of mapped single end reads (- for stdin)',
                       )
    call.add_argument('--min_quality',
                      type=int,
                      default=0,
                      help='Minimum quality score required to call a variant'
                      )
    call.add_argument('infile',
                      type=str,
                      help='a SAM or BAM format file of mapped single end reads (- for stdin)',
                      )
    if isinstance(arguments, str):
        arguments = shlex.split(arguments)
    return parser.parse_args(arguments)


def main(arguments: Optional[Union[str, Sequence[str]]] = None):
    args = cli(arguments)
    # htslib logs an error when a BAM has no index, which pebbles never needs,
    # and duplicates the open failures reported below, so silence it while opening
    htslib_verbosity = pysam.set_verbosity(0)
    try:
        # htslib detects SAM or BAM from the content, so this also works for stdin ('-')
        infile = pysam.AlignmentFile(args.infile, "r")
    except (OSError, ValueError) as error:
        sys.exit(f"pebbles: error: cannot read alignments from {args.infile}: {error}")
    finally:
        pysam.set_verbosity(htslib_verbosity)

    if args.command == 'call':
        print('readname\tvariants')
        for x in call_mutations_from_pysam(infile, min_quality=args.min_quality):
            print("\t".join([str(_) for _ in x]))
    elif args.command == 'count':
        print('variant\tcount')
        print(count(infile, max_variants=args.max, min_quality=args.min_quality), end='')


if __name__ == '__main__':
    main()
