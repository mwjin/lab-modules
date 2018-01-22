import sys, os

class NarrowPeak:
    """ The object of this class represents one entry of a file that has narrow peaks BED format """
    def __init__(self):
        pass

    def __str__(self):
        pass

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


