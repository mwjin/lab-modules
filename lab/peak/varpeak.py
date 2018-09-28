from lab.peak.narrow import NarrowPeak
from lab.gene.anno import genic_region_list, parse_anno_val

__all__ = ['VarPeak']


class VarPeak(NarrowPeak):
    """ The object of this class is a NarrowPeak's child class that can contain mutations """
    _genic_regions = genic_region_list()  # for the attributes in this class

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # attributes for the variants on this peak
        self.var_pos_to_cnt = {}  # positions of variants (0-based) to their counts
        self.var_pos_to_anno_val = {}  # positions of variants to their annotation values
        self.var_pos_to_genes = {}  # positions of variants to their associated genes (type: dictionary)

        # attributes for the gene-based annotation of the peak
        self.anno_vals = []
        self.genic_region_to_size = {genic_region: 0 for genic_region in self._genic_regions}
        self.genic_region_to_var_cnt = {genic_region: 0 for genic_region in self._genic_regions}
        self.repr_genic_region_to_size = {genic_region: 0 for genic_region in self._genic_regions}
        self.repr_genic_region_to_var_cnt = {genic_region: 0 for genic_region in self._genic_regions}

    def combine(self, other):
        """
        Combine the distribution of variants on other peak with this peak.
        :param other: an object of 'VarPeak'
        """
        assert other.__class__.__name__ == self.__class__.__name__
        assert other == self

        for var_pos in other.var_pos_to_cnt:
            # combine the attributes for the variants on this peak
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

            # combine the attributes for the gene-based annotation of this peak
            anno_dict = parse_anno_val(anno_val)

            for genic_region in self._genic_regions:
                if anno_dict[genic_region]:
                    self.genic_region_to_var_cnt[genic_region] += 1

            for genic_region in self._genic_regions:
                if anno_dict[genic_region]:
                    self.repr_genic_region_to_var_cnt[genic_region] += 1

                    if genic_region.startswith('5') and anno_dict['3UTR'] is True:  # both UTRs have same priority
                        self.repr_genic_region_to_var_cnt['3UTR'] += 1
                    break

    def cut(self, start, end):
        """
        cut the peak and return new object
        :param start: a start position of the new object
        :param end: an end position of the new object
        :return: a 'VarPeak' object
        """
        assert self.start <= start < end <= self.end
        new_peak = VarPeak(self.chrom, start, end, self.strand)

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

    def merge(self, peak):
        pass

    def intersect(self, peak):
        pass

    def get_genic_region_to_size(self, only_repr=False):
        """
        :param only_repr: if True, then return repr_genic_region_to_size
        :return: a dictionary that documents the size of each genic region on this peak.
                 (key: a genic region, value: a size of the genic region (integer))
        """
        if only_repr:
            return self.repr_genic_region_to_size
        else:
            return self.genic_region_to_size

    def get_genic_region_to_var_cnt(self, only_repr=False):
        """
        :param only_repr: if True, then return repr_genic_region_to_var_cnt
        :return: a dictionary that documents the number of variants of each genic region on this peak.
                 (key: a genic region, value: the number of variants on the genic region (integer))
        """
        if only_repr:
            return self.repr_genic_region_to_var_cnt
        else:
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

    def gene_based_annotation(self, anno_val_list):
        """
        This function enters the value of the attribute 'anno_vals' and
        makes up the 'genic_region_to_size' and 'repr_genic_region_to_size' attribute.
        * representative genic region: a genic region which has the highest priority among genic region candidates
        :param anno_val_list: a list of anntation values corresponding to this peak (see gene.utils)
        """
        assert len(anno_val_list) == self.get_size()
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

    def put_variant(self, variant):
        """
        :param variant: an object of the class 'VCFData'
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

        # save the information of the variant-associated genes
        if var_pos not in self.var_pos_to_genes:
            assoc_genes = variant.get_assoc_genes(self.strand)
            self.var_pos_to_genes[var_pos] = assoc_genes
        else:
            assert self.var_pos_to_genes[var_pos] == variant.get_assoc_genes(self.strand)

        # make up the attributes for the gene-based annotation
        anno_dict = parse_anno_val(self.var_pos_to_anno_val[var_pos])

        for genic_region in self._genic_regions:
            if anno_dict[genic_region]:
                self.genic_region_to_var_cnt[genic_region] += 1

        for genic_region in self._genic_regions:
            if anno_dict[genic_region]:
                self.repr_genic_region_to_var_cnt[genic_region] += 1

                if genic_region.startswith('5') and anno_dict['3UTR'] is True:  # both UTRs have same priority
                    self.repr_genic_region_to_var_cnt['3UTR'] += 1
                break

