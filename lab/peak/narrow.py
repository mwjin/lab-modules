"""
'NarrowPeak' class that represents one entry of the 'NarrowPeak' BED file
"""
__all__ = ['NarrowPeak']


class NarrowPeak:
    """ The object of this class represents one entry of a file that has Narrow Peaks format """
    def __init__(self, *args, **kwargs):
        """
        :param args: Its length must be 0, 3, or 4. If the length is 3 or 4, args must consist of the following values
            1. the chromosome ID (string)
            2. start position (integer)
            3. end position (integer)
            4. strand (+ or -) (string) (if the length is 3, there is no input for strand)
            If the length is 0, default values will be entered.

        :param kwargs: The key values must be matched with one of the attributes in this class.
        """
        argc = len(args)
        assert argc == 0 or argc == 3 or argc == 4

        self.chrom = None
        self.start = 0  # 0-based
        self.end = 0
        self.strand = '.'

        if argc != 0:  # 3 or 4
            self.chrom = args[0]
            self.start = int(args[1])
            self.end = int(args[2])

            if argc == 4:
                self.strand = args[3]

        self.name = '.'
        self.score = '0'
        self.sig_val = '0.0'
        self.p_val = '-1.0'
        self.q_val = '-1.0'
        self.point_src = '-1'

        for attr in kwargs:
            if hasattr(self, attr):
                setattr(self, attr, kwargs[attr])
            else:
                raise AttributeError('Error: %s has no attribute named as %s.' % (self.__class__.__name__, attr))

    # END: __init__

    def __str__(self):
        return '%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % \
               (self.chrom, self.start, self.end, self.name, self.score,
                self.strand, self.sig_val, self.p_val, self.q_val, self.point_src)

    # END: __str__

    def __eq__(self, other):
        self_info = (self.chrom, self.start, self.end, self.name, self.score, self.strand)
        other_info = (other.chrom, other.start, other.end, other.name, other.score, other.strand)

        return self_info == other_info

    # END: __eq__

    def get_position(self):
        """
        :return: a tuple (start, end)
        """
        return self.start, self.end

    # END: the function 'get_position'

    def merge(self, peak):
        """
        merge the position of this peak with the input peak and return the peak which has the merged position
        :param peak: an object of the same class which chromosomal position and strand are same with this peak
        :return: an object of the 'NarrowPeak' or None if there is a gap between the two peaks
        """
        assert peak.__class__.__name__ == self.__class__.__name__
        assert peak.chrom == self.chrom
        assert peak.strand == self.strand

        pos_list = []

        peak_start, peak_end = peak.get_position()
        pos_list.append((peak_start, peak_end))
        pos_list.append((self.start, self.end))
        pos_list.sort(key=lambda pos: (pos[0], pos[1]))

        if pos_list[1][0] - pos_list[0][1] > 0:
            return None
        else:
            return NarrowPeak(self.chrom, pos_list[0][0], pos_list[1][1], self.strand)

    # END: the function 'merge'

    def intersect(self, peak):
        """
        intersect the position of this peak with the input peak and return the peak which has the intersected position
        :param peak: an object of the same class which chromosomal position and strand are same with this peak
        :return: an object of the 'NarrowPeak' or None if there is no overlap
        """
        assert peak.__class__.__name__ == self.__class__.__name__
        assert peak.chrom == self.chrom
        assert peak.strand == self.strand

        pos_list = []

        peak_start, peak_end = peak.get_position()
        pos_list.append((peak_start, peak_end))
        pos_list.append((self.start, self.end))
        pos_list.sort(key=lambda pos: (pos[0], pos[1]))

        if pos_list[1][0] - pos_list[0][1] >= 0:
            return None
        else:
            return NarrowPeak(self.chrom, pos_list[1][0], pos_list[0][1], self.strand)

    # END: the function 'intersect'

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
        field_num = len(fields)
        assert field_num == 6 or field_num == 10

        # default
        self.chrom = fields[0]
        self.start = int(fields[1])
        self.end = int(fields[2])
        self.name = fields[3]
        self.score = fields[4]
        self.strand = fields[5]

        # optional
        if field_num == 10:
            self.sig_val = fields[6]
            self.p_val = fields[7]
            self.q_val = fields[8]
            self.point_src = fields[9]

    # END: the function 'parse_peak_entry'

    @classmethod
    def parse_peak_file(cls, peak_filename):
        """
        :param peak_filename: a file which has a narrow peak bed format
        :return: a list of NarrowPeak objects
        """
        peak_file = open(peak_filename, 'r')

        peak_list = []

        for line in peak_file.readlines():
            peak = cls()
            peak.parse_peak_entry(line.strip())
            peak_list.append(peak)

        peak_file.close()

        return peak_list

    # END: the function 'parse_peak_file'

# END: the definition of the class 'NarrowPeak'
