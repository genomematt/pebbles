# Pebbles. Sister of BammBamm Flinstone
# A per read bam mutation caller

import re
import pysam

# QNAME
# FLAG
# POS
# MAPQ
# CIGAR
# SEQ
# QUAL
# MD


def expand_cigar(cigar):
    return "".join([x[1] * int(x[0]) for x in re.findall(r'([0-9]+)([MIDNSHPX=])', cigar)])

def compact_cigar(expanded_cigar):
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
          is_reference = False):
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

def cigar_trimmer(cigar, trim_from_start=0, trim_from_end=0):
    xcigar = expand_cigar(cigar)
    result = []
    sequence_length = len(xcigar) - xcigar.count('D')
    position_in_sequence = 0
    for state in xcigar:
        if not (state == 'D' or state == 'S'):
            position_in_sequence += 1
        if position_in_sequence > trim_from_start and position_in_sequence <= sequence_length - trim_from_end:
            result.append(state)
    return compact_cigar(result)

def fix_softmasked_expanded_cigar(expanded_cigar, start=0):
    # fix SDS and SIS states
    if isinstance(expanded_cigar, list):
        expanded_cigar = ''.join(expanded_cigar)
    if 'S' in expanded_cigar:
        ##These versions trim insert delete states adjacent to softmask
        ##This produces more valid cigar strings, but may softmask variants
        # initial_SDI_block = re.findall(r'^([SDI]*)',expanded_cigar)[0]
        # final_SDI_block = re.findall(r'([SDI]*$)',expanded_cigar)[0]
        leading_matches = re.findall(r'^([SDI]*S)', expanded_cigar)
        initial_SDI_block = leading_matches[0] if leading_matches else ''
        trailing_matches = re.findall(r'(S[SDI]*$)', expanded_cigar)
        final_SDI_block = trailing_matches[0] if trailing_matches else ''
        initial_softmask_string = initial_SDI_block.replace('D', '').replace('I', 'S')
        start += len(initial_softmask_string)
        final_softmask_string = final_SDI_block.replace('D', '').replace('I', 'S')
        length_in_ref = len(expanded_cigar) - (len(initial_SDI_block) + len(final_SDI_block))
        fixed_expanded_cigar = initial_softmask_string + expanded_cigar[
                                                         len(initial_SDI_block):len(expanded_cigar) - len(
                                                             final_SDI_block)] + final_softmask_string
        ### TODO find any MSM blocks
        ### TODO replace internal S with M
    else:
        fixed_expanded_cigar = expanded_cigar
        length_in_ref = len([x for x in expanded_cigar if x in 'MD'])
    return fixed_expanded_cigar, start, length_in_ref

def gapped_alignment_to_cigar(aligned_reference, aligned_sample, gap='-', snv='M'):
    if len(aligned_reference) != len(aligned_sample):
        raise RuntimeError(
            'Unequal sequences lengths - not correctly aligned \n    {0}\n    {1}'.format(aligned_reference,
                                                                                          aligned_sample))
    xcigar = []
    for i in range(len(aligned_reference)):
        if aligned_reference[i] == aligned_sample[i]:
            if aligned_reference[i] == gap:
                raise RuntimeError('Invalid alignment state \n{0}\n{1}'.format(aligned_reference[i], aligned_sample[i]))
            xcigar.append('M')
        elif aligned_reference[i] in ascii_uppercase and \
                aligned_sample[i] in ascii_uppercase:
            xcigar.append(snv)
        elif aligned_reference[i] in ascii_lowercase:
            if not aligned_sample[i] == gap:
                xcigar.append('S')
            else:
                xcigar.append('D')
        elif aligned_reference[i] == gap:
            xcigar.append('I')
        elif aligned_sample[i] == gap:
            xcigar.append('D')
        else:  # pragma no cover # I cant think of a way to get here but complex so dont want to default to a valid state
            raise RuntimeError('Invalid alignment state \n{0}\n{1}'.format(aligned_reference[i], aligned_sample[i]))
    fixed_xcigar, start, length = fix_softmasked_expanded_cigar(xcigar)
    return compact_cigar(fixed_xcigar), start, length


def expand_mdtag(mdtag):
    mdtag_tokens = re.findall(r'(\d+|\D+)', mdtag)
    result = []
    for x in mdtag_tokens:
        if x[0] in '1234567890':
            result.extend(['.', ] * int(x))
        elif x[0] == '^':
            result.extend(list(x[1:]))
        else:
            result.extend(list(x))
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
                   gapped_read,
                   gapped_quality = None,):
    for i in range(len(gapped_read)):
        print(i,gapped_read[i], expanded_engapped_md[i], expanded_cigar[i])




if __name__ == '__main__':


    samfile = pysam.AlignmentFile("tests/data/map.sam", "r")
    for segment in samfile:
        if segment.is_qcfail or segment.is_supplementary or segment.is_unmapped:
            continue
        try:
            mdtag = segment.get_tag('MD')
        except:
            print(f'skipping {segment.qname} as it has no MD tag')
            #skip reads with no MD tag
            continue
        print(segment.qname)
        #print(mdtag, expand_mdtag(mdtag))
        #print(segment.cigarstring, expand_cigar(segment.cigarstring))
        #print(segment.reference_name, segment.pos)

        call_mutations(segment.reference_name,
                       segment.pos,
                       engap(expand_mdtag(mdtag), segment.cigarstring, is_reference = True),
                       expand_cigar(segment.cigarstring),
                       engap(segment.seq, segment.cigarstring))