from lab.peak.narrow import NarrowPeak
from lab.gene.anno import genic_region_list, parse_anno_val

__all__ = ['RBPPeak']


class RBPPeak(NarrowPeak):
    """ The object in this class represents one peak calling of the binding sites of RNA-binding proteins """
    _genic_regions = genic_region_list()  # for the attributes in this class

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # attributes for the gene-based annotation of the peak
        self.anno_vals = []
        self.genic_region_to_size = {genic_region: 0 for genic_region in self._genic_regions}
        self.genic_region_to_var_cnt = {genic_region: 0 for genic_region in self._genic_regions}
        self.repr_genic_region_to_size = {genic_region: 0 for genic_region in self._genic_regions}

        # attributes for the variants on this peak
        self.var_pos_to_cnt = {}  # positions of variants (0-based) to their counts
        self.var_pos_to_anno_val = {}  # positions of variants to their genic region values
        self.var_pos_to_genes = {}  # positions of variants to their associated genes (type: dictionary)

    def combine(self, other):
        """
        Combine the distribution of variants on other peak with this peak.
        :param other: an object of 'RBPPeak'
        """
        assert other.__class__.__name__ == self.__class__.__name__
        assert other == self

        for var_pos in other.var_pos_to_cnt:
            var_cnt = other.get_var_cnt_in_pos(var_pos)
            anno_val = other.get_anno_val_in_pos(var_pos)
            genes = other.get_genes_in_pos(var_pos)

            if var_pos not in self.var_pos_to_cnt:
                self.var_pos_to_anno_val[var_pos] = anno_val
                self.var_pos_to_genes[var_pos] = genes
                self.var_pos_to_cnt[var_pos] = var_cnt
            else:
                assert self.var_pos_to_anno_val[var_pos] == anno_val
                assert self.var_pos_to_genes[var_pos] == genes
                self.var_pos_to_cnt[var_pos] += var_cnt

    def cut(self, start, end):
        """
        cut the peak and return new object
        :param start: a start position of the new object
        :param end: an end position of the new object
        :return: a 'RBPPeak' object
        """
        assert self.start <= start < end <= self.end
        new_peak = RBPPeak(self.chrom, self.start, self.end, self.strand)

        var_pos_list = self.get_var_pos_list()

        for var_pos in var_pos_list:
            if start <= var_pos < end:
                var_cnt = self.var_pos_to_cnt[var_pos]
                var_anno_val = self.var_pos_to_anno_val[var_pos]
                assoc_genes = self.var_pos_to_genes[var_pos]

                new_peak.var_pos_to_cnt[var_pos] = var_cnt
                new_peak.var_pos_to_anno_val[var_pos] = var_anno_val
                new_peak.var_pos_to_genes[var_pos] = assoc_genes

                genic_region_dict = parse_anno_val(var_anno_val)

                for genic_region in genic_region_dict:
                    if genic_region_dict[genic_region]:
                        new_peak.genic_region_to_var_cnt[genic_region] += var_cnt
            elif var_pos >= end:
                break

        return new_peak

    def get_genic_region_to_size(self, only_repr=False):
        """
        :param: if True, then return repr_genic_region_to_size
        :return: a dictionary that documents the size of each genic region on this peak.
                 (key: a genic region, value: a size of the genic region (integer))
        """
        if only_repr:
            return self.repr_genic_region_to_size
        else:
            return self.genic_region_to_size

    def get_genic_region_to_var_cnt(self):
        """
        :return: a dictionary that documents the number of variants of each genic region on this peak.
                 (key: a genic region, value: the number of variants on the genic region (integer))
        """
        return self.genic_region_to_var_cnt

    def get_var_pos_list(self):
        """
        :return: a list of the positions (0-based) of all variants on this peak
        """
        return sorted(list(self.var_pos_to_cnt.keys()))

    def get_var_cnt(self):
        """
        :return: the number of variants in this peak
        """
        var_cnt = 0

        for var_pos in self.var_pos_to_cnt:
            var_cnt += self.var_pos_to_cnt[var_pos]

        return var_cnt

    def get_var_cnt_in_pos(self, var_pos):
        """
        :param var_pos: a variant position (0-based) inside the peak
        :return: the number of variants on the position
        """
        var_cnt = self.var_pos_to_cnt.get(var_pos)
        assert var_cnt is not None

        return var_cnt

    def get_anno_val_in_pos(self, var_pos):
        """
        :param var_pos: a variant position (0-based) inside the peak
        :return: a genic region value of the position
        """
        anno_val = self.var_pos_to_anno_val.get(var_pos)
        assert anno_val is not None

        return anno_val

    def get_genes_in_pos(self, var_pos):
        """
        :param var_pos: a variant position (0-based) inside the peak
        :return: a dictionary (mapping route: gene symbol -> ID -> genic region)
        """
        return self.var_pos_to_genes[var_pos]

    def get_anno_vals(self):
        """
        Return the annotation values of this peak (each element corresponds to each nucleotide on this peak).
        If it is not set, return an empty list.
        """
        return self.anno_vals

    def set_genic_region_size(self, anno_val_list):
        """
        This code makes up the 'genic_region_to_size' and 'repr_genic_region_to_size' attribute.
        * representative genic region: a genic region which has the highest priority among genic region candidates
        :param anno_val_list: a list of genic region values (see gene.utils)
        """
        self.anno_vals = anno_val_list

        # make a statistics for genic regions
        for anno_val in anno_val_list:
            anno_dict = parse_anno_val(anno_val)

            for genic_region in self._genic_regions:
                if anno_dict[genic_region]:
                    self.genic_region_to_size[genic_region] += 1

        # make a statistics for representative genic regions
        for anno_val in anno_val_list:
            anno_dict = parse_anno_val(anno_val)

            for genic_region in self._genic_regions:
                if anno_dict[genic_region]:
                    self.repr_genic_region_to_size[genic_region] += 1

                    # 5UTR and 3UTR have same priority.
                    if genic_region.startswith('5') and anno_dict['3UTR'] is True:
                        self.repr_genic_region_to_size['3UTR'] += 1
                    break

    def put_variant(self, variant, only_repr=False):
        """
        :param variant: an object of the class 'VCFData'
        :param only_repr: if it is true, consider only the representative genic region
                        when making up the self.genic_region_to_var_cnt.

        * representative genic region: a genic region which has the highest priority among genic region candidates
        """

        assert variant.__class__.__name__ == 'VCFData'
        assert self.start <= variant.pos < self.end

        var_pos = variant.pos

        # variant counting
        if var_pos not in self.var_pos_to_cnt:
            self.var_pos_to_cnt[var_pos] = 1
        else:
            self.var_pos_to_cnt[var_pos] += 1

        # gene-based annotation of the variants
        if var_pos not in self.var_pos_to_anno_val:
            var_anno_val = variant.get_anno_val(self.strand)
            self.var_pos_to_anno_val[var_pos] = var_anno_val
        else:
            assert self.var_pos_to_anno_val[var_pos] == variant.get_anno_val(self.strand)

        anno_dict = parse_anno_val(self.var_pos_to_anno_val[var_pos])

        if only_repr:
            for genic_region in self._genic_regions:
                if anno_dict[genic_region]:
                    self.genic_region_to_var_cnt[genic_region] += 1

                    if genic_region.startswith('5') and anno_dict['3UTR'] is True:
                        self.genic_region_to_var_cnt['3UTR'] += 1

                    break
        else:
            for genic_region in self._genic_regions:
                if anno_dict[genic_region]:
                    self.genic_region_to_var_cnt[genic_region] += 1

        # save the information of the variant-associated genes
        if var_pos not in self.var_pos_to_genes:
            assoc_genes = variant.get_assoc_genes(self.strand)
            self.var_pos_to_genes[var_pos] = assoc_genes
        else:
            assert self.var_pos_to_genes[var_pos] == variant.get_assoc_genes(self.strand)
