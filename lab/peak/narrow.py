

class NarrowPeak:
    """ The object of this class represents one entry of a file that has Narrow Peaks format """
    def __init__(self, *args, **kwargs):
        """
        :param args: Its length must be 0 or 3. If the length is 3, args must consist of the following values
            1. the chromosome ID (string)
            2. start position (integer)
            3. end position (integer)
            If the length is 0, default values will be entered.

        :param kwargs: The key values must be matched with one of the attributes in this class.
        """
        argc = len(args)
        assert argc == 0 or argc == 3

        if argc == 0:
            self.chrom = None
            self.start = 0  # 0-based
            self.end = 0
        else:
            self.chrom = args[0]
            self.start = int(args[1])
            self.end = int(args[2])

        self.name = '.'
        self.score = '0'
        self.strand = '.'
        self.sig_val = '0.0'
        self.p_val = '-1.0'
        self.q_val = '-1.0'
        self.point_src = '-1'

        for attr in kwargs:
            if hasattr(self, attr):
                setattr(self, attr, kwargs[attr])
            else:
                raise AttributeError('Error: %s has no attribute named as %s.' % (self.__class__, attr))

    # END: __init__

    def __str__(self):
        return '%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % \
               (self.chrom, self.start, self.end, self.name, self.score,
                self.strand, self.sig_val, self.p_val, self.q_val, self.point_src)

    # END: __str__

    def get_position(self):
        """
        :return: a tuple (start, end)
        """
        return self.start, self.end

    # END: the function 'get_position'

    def parse_peak_entry(self, peak_entry):
        """
        :param peak_entry: an entry from a narrow peak bed file
        """
        # File Format
        # Column:        0       | 1         | 2        | 3               | 4      | 5       |
        # ID:            chrom   | start     | end      | name            | score  | strand  |
        # Example:       chr14   | 56879239  | 56879435 | ILF3_K562_rep02 | 1000   | -       |
        # Column:        6            | 7         | 8        | 9
        # ID:            signal_value | p-value   | q-value  | point_src
        # Example:       1.29065      | 0.198802  | -1       | -1

        fields = peak_entry.split('\t')

        self.chrom = fields[0]
        self.start = int(fields[1])
        self.end = int(fields[2])
        self.name = fields[3]
        self.score = fields[4]
        self.strand = fields[5]
        self.sig_val = fields[6]
        self.p_val = fields[7]
        self.q_val = fields[8]
        self.point_src = fields[9]

    # END: the function 'parse_peak_entry'

    @staticmethod
    def parse_peak_file(peak_filename):
        """
        :param peak_filename: a file which has a narrow peak bed format
        :return: a list of NarrowPeak objects
        """
        peak_file = open(peak_filename, 'r')

        peak_list = []

        for line in peak_file.readlines():
            peak = NarrowPeak()
            peak.parse_peak_entry(line.strip())
            peak_list.append(peak)

        peak_file.close()

        return peak_list

    # END: the function 'parse_peak_file'

# END: the definition of the class 'NarrowPeak'
