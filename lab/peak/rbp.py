from lab.peak.narrow import NarrowPeak
from lab.gene.anno import genic_region_list, parse_genic_region_val

__all__ = ['RBPPeak']


class RBPPeak(NarrowPeak):
    """ The object in this class represents one peak calling of the binding sites of RNA-binding proteins """
    _genic_regions = genic_region_list()  # for the attributes in this class

    def __init__(self):
        super().__init__()

        # attributes for the gene-based annotation of the peak
        self.genic_region_to_size = {genic_region: 0 for genic_region in self._genic_regions}
        self.genic_region_to_var_cnt = {genic_region: 0 for genic_region in self._genic_regions}

        # attributes for the variants on this peak
        self.var_pos_to_cnt = {}  # positions of variants (0-based) to their counts
        self.var_pos_to_region_val = {}  # positions of variants to their genic region values

    # END: __init__

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

    def set_genic_region_size(self, genic_region_val_list, repr_only=False):
        """
        This code makes up the 'genic_region_to_size' attribute.
        :param genic_region_val_list: a list of genic region values (see gene.utils)
        :param repr_only: if it is true, consider only the representative genic region
                        when making up the self.genic_region_to_size.

        * representative genic region: a genic region which has the highest priority among genic region candidates
        """
        # make a statistics for genic regions
        if repr_only:
            for region_val in genic_region_val_list:
                genic_region_to_bool = parse_genic_region_val(region_val)

                for genic_region in self._genic_regions:
                    if genic_region_to_bool[genic_region]:
                        self.genic_region_to_size[genic_region] += 1

                        # 5UTR and 3UTR have same priority.
                        if genic_region.startswith('5') and genic_region_to_bool['3UTR'] is True:
                            self.genic_region_to_size['3UTR'] += 1

                        break
        else:
            for region_val in genic_region_val_list:
                genic_region_to_bool = parse_genic_region_val(region_val)

                for genic_region in self._genic_regions:
                    if genic_region_to_bool[genic_region]:
                        self.genic_region_to_size[genic_region] += 1

    # END: the function 'set_genic_region_size'

    def put_variant(self, variant, repr_only=False):
        """
        :param variant: an object of the class 'VCFData'
        :param repr_only: if it is true, consider only the representative genic region
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

            if repr_only:
                for genic_region in self._genic_regions:
                    if genic_region_to_bool[genic_region]:
                        self.genic_region_to_var_cnt[genic_region] += 1

                        if genic_region.startswith('5') and genic_region_to_bool['3UTR'] is True:
                            self.genic_region_to_var_cnt['3UTR'] += 1

                        break
            else:
                for genic_region in self._genic_regions:
                    if genic_region_to_bool[genic_region]:
                        self.genic_region_to_var_cnt[genic_region] += 1

        else:
            assert self.var_pos_to_region_val[var_pos] == variant.get_var_genic_region(self.strand)

    # END: the function 'put_variant'

# END: the definition of the class 'RBPPeak'
