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
