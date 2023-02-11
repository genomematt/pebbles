# Pebbles. Sister of BammBamm Flinstone
# A per read bam mutation caller

import re
import pysam


def expand_cigar(cigar):
    return "".join([match[1] * int(match[0]) for match in re.findall(r'([0-9]+)([MIDNSHPX=])', cigar)])


def compact_cigar(expanded_cigar):
    state = None
    result = []
    last_state = expanded_cigar[0]
    count = 0
    for state in expanded_cigar:
        if state == last_state:
            count += 1
        else:
            result.append(str(count) + last_state)
            last_state = state
            count = 1
    result.append(str(count) + state)
    return "".join(result)


def engap(seq,
          cigar,
          is_reference: bool = False):
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


def expand_mdtag(mdtag):
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


def compact_expanded_mdtag_tokens(expanded_mdtag_tokens):
    result = []
    count = 0
    in_deletion = False
    for token in expanded_mdtag_tokens:
        if token == '':
            count += 1
            in_deletion = False
        elif count:
            # exiting match
            result.append(str(count))
            count = 0
            if token == '^':
                in_deletion = True
            result.append(token)
        else:
            # in mismatch or deletion states
            if in_deletion and token == '^':
                # adjacent deletions to be merged
                continue
            if in_deletion and (token in 'CAGTN'):
                # have a mismatch adjacent to a deletion
                result.append('0')
                in_deletion = False
            result.append(token)
    if count:
        result.append(str(count))
    return "".join(result).upper()


def call_mutations(refname,
                   pos,
                   expanded_engapped_md,
                   expanded_cigar,
                   gapped_read):
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
        if expanded_cigar[i] == 'M' and expanded_engapped_md[i] == '.':
            # match state with no variants
            i += 1
            continue
        if expanded_cigar[i] == 'D':
            deleted = ''
            delstart = i+pos+1+nonref_bases-softmasked
            while expanded_cigar[i] == 'D':
                deleted += expanded_engapped_md[i]
                i += 1
            mutations.append(f'{refname}:g.{delstart}_{i+pos+nonref_bases-softmasked}del{deleted}')
        if expanded_cigar[i] == 'I':
            inserted = ''
            insstart = i + pos + nonref_bases - softmasked  # base before first event base
            while expanded_cigar[i] == 'I':
                inserted += gapped_read[i]
                i += 1
                nonref_bases += 1
            mutations.append(f'{refname}:g.{insstart}_{insstart+1}ins{inserted}')
        if expanded_cigar[i] == 'M' and expanded_engapped_md[i] != '.':
            mutant = ''
            reference = ''
            substart = i + pos + 1 + nonref_bases - softmasked
            while expanded_cigar[i] == 'M' and expanded_engapped_md[i] != '.':
                mutant += gapped_read[i]
                reference += expanded_engapped_md[i]
                i += 1
            if len(mutant) == 1:
                mutations.append(f'{refname}:g.{substart}{reference}>{mutant}')
            else:
                mutations.append(f'{refname}:g.{substart}_{i+pos+nonref_bases-softmasked}delins{mutant}')
        i += 1
        return mutations


def call_mutations_from_pysam(pysamfile,):
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
        yield (segment.qname, mutations)


if __name__ == '__main__':
    for x in call_mutations_from_pysam(pysam.AlignmentFile("tests/data/map.sam", "r")):
        print(x)
