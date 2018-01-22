from lib_utils import *

import sys, os

class NarrowPeak:
    """ The object of this class represents one entry of a file that has narrow peaks BED format """
    def __init__(self):
        self.chrID = None
        self.start = 0
        self.end = 0
        self.name = '.'
        self.score = 0
        self.strand = '.'
        self.sig_val = 0.0
        self.p_val = -1.0
        self.q_val = -1.0
        self.point_source = -1

        """ info for gene-based annotation of the peak """
        self.genic_region_to_size = {}

    def __str__(self):
        return '%s\t%d\t%d\t%s\t%d\t%s\t%f\t%f\t%f\t%d' % \
               (self.chrID, self.start, self.end, self.name, self.score,
                self.strand, self.sig_val, self.p_val, self.q_val, self.point_source)

    @staticmethod
    def parse_peak_file():
        # File Format:
        # Column:        0       | 1         | 2        | 3               | 4      | 5       |
        # ID:            ChrID   | start     | end      | name            | score  | strand  |
        # Example:       chr14   | 56879239  | 56879435 | ILF3_K562_rep02 | 1000   | -       |
        # Column:        6            | 7         | 8        | 9
        # ID:            signal_value | p-value   | q-value  | point_source
        # Example:       1.29065      | 0.198802  | -1       | -1

        pass

    def parse_genic_region_info(self, region_code_list):
        """
        This code makes up the 'genic_region_to_size' attribute.

        :param region_code_list: a list that has a same length with the peak and
                                 consists of codes representing each genic region
                                 (100, 101, ..., 107)

        ** codes for each genic region **
        100: ORF
        101: 5UTR
        102: 3UTR
        103: UTR (both 5UTR and 3UTR)
        104: ncRNA exonic
        105: intronic
        106: ncRNA intronic
        107: intergenic
        """
        assert len(region_code_list) == (self.end - self.start)

        code_to_region = {100: 'ORF', 101: '5UTR', 102: '3UTR', 103: 'UTR', 104: 'ncRNA_exonic',
                          105: 'intronic', 106: 'ncRNA_intronic', 107: 'intergenic'}

        for code in region_code_list:
            try:
                genic_region = code_to_region[code]
            except KeyError:
                print('Error in %s: invalid code for genic region' % caller_file_and_line())
                sys.exit()

            if genic_region not in self.genic_region_to_size:
                self.genic_region_to_size[genic_region] = 0

            self.genic_region_to_size[genic_region] += 1

