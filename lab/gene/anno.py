"""
Module for gene-based annotation
1. a list of genic regions we use
2. variables and functions for the gene-based annotation
"""
import sys
from lab.utils import caller_file_and_line, eprint

"""
/* genic region value */
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

__all__ = ['genic_region_list', 'get_genic_region_val', 'parse_genic_region_val']

# constants used in this module
_GENIC_REGIONS = ['ORF', '5UTR', '3UTR', 'ncRNA_exonic',
                  'SS', 'SS-30nt', 'SS-50nt', 'intronic',
                  'ncRNA_intronic', 'intergenic']
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


def get_genic_region_val(genic_region):
    """
    get the genic region value of the input genic region (string)
    :param genic_region: a string that represents the genic region
    :return: an integer that represents the region value of the input genic region
    """
    try:
        bit_pos = _REGION_TO_BIT_POS[genic_region]
        return 2 ** (_BIT_LEN - bit_pos)
    except KeyError:
        eprint("Error in %s: invalid genic region '%s'" % (caller_file_and_line(), genic_region))
        sys.exit()


def parse_genic_region_val(region_val):
    """
    parse the genic region value and get information which genic region is used for the annotation
    :param region_val: an integer that represents a genic region value
    :return: a dictionary (key: genic region, value: boolean (whether the region is used for the annotation)
    """
    assert 0 <= region_val < (2 ** _BIT_LEN)  # the maximum bit length of the genic region value is 6

    genic_region_dict = {genic_region: False for genic_region in _GENIC_REGIONS}

    if region_val == 0:
        genic_region_dict['intergenic'] = True
    else:
        for bit_pos in range(1, _BIT_LEN + 1):
            bit = int(region_val / (2 ** (_BIT_LEN - bit_pos))) % 2

            if bit == 1:
                genic_region = _BIT_POS_TO_REGION[bit_pos]
                genic_region_dict[genic_region] = True

    # deal with the case of splicing sites
    if genic_region_dict['SS']:
        genic_region_dict['SS-30nt'] = True

    if genic_region_dict['SS-30nt']:
        genic_region_dict['SS-50nt'] = True

    if genic_region_dict['SS-50nt']:
        genic_region_dict['intronic'] = True

    return genic_region_dict
