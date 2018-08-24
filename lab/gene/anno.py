"""
Module for gene-based annotation
1. a list of genic regions we use
2. variables and functions for the gene-based annotation
"""
import sys
from lab.utils import caller_file_and_line, eprint

"""
/* Annotation value */
The gene-based annotation of one nucleotide is represented as integer which bit length is # of genic regions - 1
('intergenic' is excluded).
Each bit of the integer represents a boolean value for one genic region.

1st significant bit: ORF (MSB)
2nd significant bit: 5'UTR
3rd significant bit: 3'UTR
.
.
.

If the integer value is 0, it means that the nucleotide is intergenic.
"""

__all__ = ['genic_region_list', 'get_anno_priority', 'get_anno_val', 'parse_anno_val']

# constants used in this module
# SS: Splice site, 30nt at the each end of the intron
# intronic: intronic region except the splice site
# promoter: 300nt upstream for each gene
_GENIC_REGIONS = ['ORF', '5UTR', '3UTR', 'ncRNA_exonic', 'SS', 'intronic', 'ncRNA_intronic', 'promoter', 'intergenic']
_BIT_LEN = len(_GENIC_REGIONS) - 1
_REGION_TO_BIT_POS = {genic_region: (i + 1) for i, genic_region in enumerate(_GENIC_REGIONS)}
_BIT_POS_TO_REGION = {(i + 1): genic_region for i, genic_region in enumerate(_GENIC_REGIONS)}


def genic_region_list():
    """
    :return: the list of the all the names of regions used for the gene-based annotation
    """
    return list(_GENIC_REGIONS)


def get_anno_priority(genic_region):
    """
    :return: an integer, priority of the genic region
    """
    return _REGION_TO_BIT_POS[genic_region]


def get_anno_val(genic_region):
    """
    get the annotation value of the input genic region (string)
    :param genic_region: a string that represents the genic region
    :return: an integer that represents the region value of the input genic region
    """
    try:
        bit_pos = _REGION_TO_BIT_POS[genic_region]
        return 2 ** (_BIT_LEN - bit_pos)
    except KeyError:
        eprint("Error in %s: invalid genic region '%s'" % (caller_file_and_line(), genic_region))
        sys.exit()


def parse_anno_val(anno_val):
    """
    parse the annotation value and get information which genic region is used for the annotation
    :param anno_val: an integer
    :return: a dictionary (key: genic region, value: boolean (whether the region is used for the annotation)
    """
    assert 0 <= anno_val < (2 ** _BIT_LEN)  # the maximum bit length of the genic region value is 6

    anno_dict = {genic_region: False for genic_region in _GENIC_REGIONS}

    if anno_val == 0:
        anno_dict['intergenic'] = True
    else:
        for bit_pos in range(1, _BIT_LEN + 1):
            bit = int(anno_val / (2 ** (_BIT_LEN - bit_pos))) % 2

            if bit == 1:
                genic_region = _BIT_POS_TO_REGION[bit_pos]
                anno_dict[genic_region] = True

    return anno_dict
