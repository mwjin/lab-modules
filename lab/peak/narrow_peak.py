from gene.anno import *


class NarrowPeak:
    """ The object of this class represents one entry of a file that has narrow peaks BED format """
    def __init__(self):
        self.chrom = None
        self.start = 0
        self.end = 0

        # following attributes will be string because of the bed files after merged.
        self.name = '.'
        self.score = '0'
        self.strand = '.'
        self.sig_val = '0.0'
        self.p_val = '-1.0'
        self.q_val = '-1.0'
        self.point_src = '-1'

        """ attributes for gene-based annotation of the peak """
        self.genic_region_to_size = {genic_region: 0 for genic_region in GENIC_REGION_LIST}
        self.genic_region_to_var_cnt = {genic_region: 0 for genic_region in GENIC_REGION_LIST}

        """ attributes for seeing a distribution of variants """
        self.var_pos_to_cnt = {}  # positions of variants (0-based) to their counts
        self.var_pos_to_region_val = {}  # value: genic region value

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

    def get_genic_region_to_size(self):
        """
        :return: a dictionary that documents the size of each genic region on this peak.
                 (key: a genic region, value: a size of the genic region (integer))
        """
        return self.genic_region_to_size

    # END: the function 'get_genic_region_to_size'

    def get_genic_region_to_var_cnt(self):
        """
        :return: a dictionary that documents the number of variants of each genic region on this peak.
                 (key: a genic region, value: the number of variants on the genic region (integer))
        """
        return self.genic_region_to_var_cnt

    # END: the function 'get_genic_region_to_var_cnt'

    @staticmethod
    def parse_peak_file(peak_filename):
        """
        :param peak_filename: a file which has a narrow peak bed format
        :return: a list of NarrowPeak objects
        """
        # File Format:
        # Column:        0       | 1         | 2        | 3               | 4      | 5       |
        # ID:            chrom   | start     | end      | name            | score  | strand  |
        # Example:       chr14   | 56879239  | 56879435 | ILF3_K562_rep02 | 1000   | -       |
        # Column:        6            | 7         | 8        | 9
        # ID:            signal_value | p-value   | q-value  | point_src
        # Example:       1.29065      | 0.198802  | -1       | -1

        peak_file = open(peak_filename, 'r')

        peak_list = []

        for line in peak_file.readlines():
            peak = NarrowPeak()
            peak.parse_peak_entry(line.strip())
            peak_list.append(peak)

        peak_file.close()

        return peak_list

    # END: the function 'parse_peak_file'

    def parse_peak_entry(self, peak_entry):
        """
        :param peak_entry: an entry from a narrow peak bed file
        """

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

    def set_genic_region_size(self, genic_region_val_list, repr_on=False):
        """
        This code makes up the 'genic_region_to_size' attribute.
        :param genic_region_val_list: a list of genic region values (see gene.utils)
        :param repr_on: if it is true, consider only the representative genic region
                        when making up the self.genic_region_to_size.

        * representative genic region: a genic region which has the highest priority among genic region candidates
        """
        # make a statistics for genic regions
        if repr_on:
            for region_val in genic_region_val_list:
                genic_region_to_bool = parse_genic_region_val(region_val)

                for genic_region in GENIC_REGION_LIST:
                    if genic_region_to_bool[genic_region]:
                        self.genic_region_to_size[genic_region] += 1

                        # 5UTR and 3UTR have same priority.
                        if genic_region.startswith('5') and genic_region_to_bool['3UTR'] is True:
                            self.genic_region_to_size['3UTR'] += 1

                        break
        else:
            for region_val in genic_region_val_list:
                genic_region_to_bool = parse_genic_region_val(region_val)

                for genic_region in GENIC_REGION_LIST:
                    if genic_region_to_bool[genic_region]:
                        self.genic_region_to_size[genic_region] += 1

    # END: the function 'set_genic_region_size'

    def put_variant(self, variant, repr_on=False):
        """
        :param variant: an object of the class 'VCFData'
        :param repr_on: if it is true, consider only the representative genic region
                        when making up the self.genic_region_to_var_cnt.

        * representative genic region: a genic region which has the highest priority among genic region candidates
        """

        assert variant.__class__.__name__ == 'VCFData'
        assert self.start <= (variant.pos - 1) < self.end

        var_pos = variant.pos - 1  # 1-based -> 0-based

        if var_pos not in self.var_pos_to_cnt:
            self.var_pos_to_cnt[var_pos] = 0

        self.var_pos_to_cnt[var_pos] += 1

        if var_pos not in self.var_pos_to_region_val:
            var_region_val = variant.get_strand_region_val(self.strand)
            self.var_pos_to_region_val[var_pos] = var_region_val

            genic_region_to_bool = parse_genic_region_val(var_region_val)

            if repr_on:
                for genic_region in GENIC_REGION_LIST:
                    if genic_region_to_bool[genic_region]:
                        self.genic_region_to_var_cnt[genic_region] += 1

                        if genic_region.startswith('5') and genic_region_to_bool['3UTR'] is True:
                            self.genic_region_to_var_cnt['3UTR'] += 1

                        break
            else:
                for genic_region in GENIC_REGION_LIST:
                    if genic_region_to_bool[genic_region]:
                        self.genic_region_to_var_cnt[genic_region] += 1

        else:
            assert self.var_pos_to_region_val[var_pos] == variant.get_var_genic_region(self.strand)

    # END: the function 'put_variant'
# END: the definition of the class 'NarrowPeak'
