"""
Gene utilities

1. variables and functions for the gene-based annotation
"""

GENIC_REGION_LIST = ['ORF', '5UTR', '3UTR', 'ncRNA_exonic', 'intronic', 'ncRNA_intronic', 'intergenic']

"""
/* genic region value */
The gene-based annotation of one nucleotide is represented as integer which bit length is 6.
Each bit of the integer represents a boolean value for one genic region.

1st bit: ORF
2nd bit: 5'UTR
3rd bit: 3'UTR
4th bit: ncRNA exonic
5th bit: intron
6th bit: ncRNA intronic

For example, if a nucleotide is annotated as ORF and intron, an integer value that represents gene-based annotation
of the nucleotide will be 34 (0b100010). If the integer value is 0, it means that the nucleotide is intergenic.
"""

_region_to_bit_pos = {'ORF': 1,  # MSB
                      '5UTR': 2,
                      '3UTR': 3,
                      'ncRNA_exonic': 4,
                      'intronic': 5,
                      'ncRNA_intronic': 6
                      }

_bit_pos_to_region = {(i + 1): genic_region for i, genic_region in enumerate(GENIC_REGION_LIST)}


def parse_genic_region_val(region_val):
    """
    parses the genic region value and get information which genic region is annotated in the region
    :param region_val: an integer that represents a genic region value
    :return: a dictionary (key: genic region, value: boolean (whether the region is annotated or not)
    """
    assert region_val < 64  # the maximum bit length of the genic region value is 6

    genic_region_dict = {genic_region: False for genic_region in GENIC_REGION_LIST}

    if region_val == 0:
        genic_region_dict['intergenic'] = True
    else:
        for bit_pos in range(1, 7):
            bit = int(region_val / (2 ** (6 - bit_pos))) % 2

            if bit == 1:
                genic_region = _bit_pos_to_region[bit_pos]
                genic_region_dict[genic_region] = True

    return genic_region_dict

# END: the function 'parse_genic_region_val'
