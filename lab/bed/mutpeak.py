from lab.bed.narrow import NarrowPeak
from lab.gene.anno import genic_region_list, parse_anno_val, get_repr_anno

__all__ = ['MutPeak']


class MutPeak(NarrowPeak):
    """ The object of this class is a NarrowPeak's child class that can contain mutations """
    _genic_regions = genic_region_list()  # for the attributes in this class

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # attributes for the variants on this bed
        self.var_pos_to_cnt = {}  # positions of variants (0-based) to their counts
        self.var_pos_to_var_set = {}  # positions of variants to a set of their repr strings
        self.var_pos_to_anno_val = {}  # positions of variants to their annotation values
        self.var_pos_to_genes = {}  # positions of variants to their associated genes (type: dictionary)

        # attributes for the gene-based annotation of the bed
        self.anno_vals = []
        self.genic_region_to_size = {genic_region: 0 for genic_region in self._genic_regions}
        self.genic_region_to_var_cnt = {genic_region: 0 for genic_region in self._genic_regions}
        self.repr_genic_region_to_size = {genic_region: 0 for genic_region in self._genic_regions}
        self.repr_genic_region_to_var_cnt = {genic_region: 0 for genic_region in self._genic_regions}

        # attributes for the mutation types of the variants on this peak
        self.var_is_non_synonymous_dict = {}  # key: a repr string of a SNV object, value: is_non_synonymous (a boolean)
        self.var_to_mut_type_dict = {}  # key: a repr string of a SNV object, value: '_mut_type_dict' attr of the 'SNV'

    def combine(self, other):
        """
        Combine the distribution of variants on other bed with this bed.
        :param other: an object of 'VarPeak'
        """
        assert other.__class__.__name__ == self.__class__.__name__
        assert other == self

        for var_pos in other.var_pos_to_cnt:
            # combine the attributes for the variants on this bed
            var_cnt = other.get_var_cnt_at_pos(var_pos)
            anno_val = other.get_anno_val_at_pos(var_pos)
            genes = other.get_genes_at_pos(var_pos)

            if var_pos not in self.var_pos_to_cnt:
                self.var_pos_to_anno_val[var_pos] = anno_val
                self.var_pos_to_genes[var_pos] = genes
                self.var_pos_to_cnt[var_pos] = var_cnt
            else:
                assert self.var_pos_to_anno_val[var_pos] == anno_val
                assert self.var_pos_to_genes[var_pos] == genes
                self.var_pos_to_cnt[var_pos] += var_cnt

            # combine the attributes for the gene-based annotation of this bed
            anno_dict = parse_anno_val(anno_val)
            repr_is_multi, repr_genic_region = get_repr_anno(anno_val)

            for genic_region in self._genic_regions:
                if anno_dict[genic_region]:
                    self.genic_region_to_var_cnt[genic_region] += var_cnt

            if repr_is_multi:
                self.repr_genic_region_to_var_cnt['5UTR'] += var_cnt
                self.repr_genic_region_to_var_cnt['3UTR'] += var_cnt
            else:
                self.repr_genic_region_to_var_cnt[repr_genic_region] += var_cnt

    def cut(self, start, end):
        """
        cut the bed and return new object
        :param start: a start position of the new object
        :param end: an end position of the new object
        :return: a 'VarPeak' object
        """
        assert self.start <= start < end <= self.end

        rel_start = start - self.start
        rel_end = end - self.start

        cut_peak = MutPeak(self.chrom, start, end, self.strand)
        cut_anno_vals = self.anno_vals[rel_start:rel_end]
        cut_peak.gene_based_anno(cut_anno_vals)

        var_pos_list = self.get_var_pos_list()

        for var_pos in var_pos_list:
            if start <= var_pos < end:
                # make up the attributes for the variant distribution of this bed
                var_cnt = self.var_pos_to_cnt[var_pos]
                var_anno_val = self.var_pos_to_anno_val[var_pos]
                assoc_genes = self.var_pos_to_genes[var_pos]

                cut_peak.var_pos_to_cnt[var_pos] = var_cnt
                cut_peak.var_pos_to_anno_val[var_pos] = var_anno_val
                cut_peak.var_pos_to_genes[var_pos] = assoc_genes

                # make up the attributes for the gene-based annotation of this bed
                anno_dict = parse_anno_val(var_anno_val)
                repr_is_multi, repr_genic_region = get_repr_anno(var_anno_val)

                for genic_region in anno_dict:
                    if anno_dict[genic_region]:
                        cut_peak.genic_region_to_var_cnt[genic_region] += var_cnt

                if repr_is_multi:
                    cut_peak.repr_genic_region_to_var_cnt['5UTR'] += var_cnt
                    cut_peak.repr_genic_region_to_var_cnt['3UTR'] += var_cnt
                else:
                    cut_peak.repr_genic_region_to_var_cnt[repr_genic_region] += var_cnt

            elif var_pos >= end:
                break

        return cut_peak

    # TODO: implementation of the overrided methods 'merge' and 'cut'

    def put_variant(self, variant):
        """
        Make up the distribution of variants on this bed
        :param variant: an object of the class 'SNV'
        * representative genic region: a genic region which has the highest priority among genic region candidates
        """
        assert variant.__class__.__name__ == 'SNV'
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
        repr_is_multi, repr_genic_region = get_repr_anno(self.var_pos_to_anno_val[var_pos])

        for genic_region in self._genic_regions:
            if anno_dict[genic_region]:
                self.genic_region_to_var_cnt[genic_region] += 1

        if repr_is_multi:
            self.repr_genic_region_to_var_cnt['5UTR'] += 1
            self.repr_genic_region_to_var_cnt['3UTR'] += 1
        else:
            self.repr_genic_region_to_var_cnt[repr_genic_region] += 1

    def gene_based_anno(self, anno_val_list):
        """
        This function enters the value of the attribute 'anno_vals' and
        makes up the 'genic_region_to_size' and 'repr_genic_region_to_size' attribute.
        * representative genic region: a genic region which has the highest priority among genic region candidates
        :param anno_val_list: a list of anntation values corresponding to this bed (see gene.utils)
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
            repr_is_multi, repr_genic_region = get_repr_anno(anno_val)

            if repr_is_multi:
                self.repr_genic_region_to_size['5UTR'] += 1
                self.repr_genic_region_to_size['3UTR'] += 1
            else:
                self.repr_genic_region_to_size[repr_genic_region] += 1

    def get_var_pos_list(self):
        """
        :return: a list of the positions (0-based) of all variants on this bed
        """
        return sorted(list(self.var_pos_to_cnt.keys()))

    def get_var_cnt(self):
        """
        :return: the number of variants in this bed
        """
        var_cnt = 0

        for var_pos in self.var_pos_to_cnt:
            var_cnt += self.var_pos_to_cnt[var_pos]

        return var_cnt

    def get_var_cnt_at_pos(self, var_pos):
        """
        :param var_pos: a variant position (0-based) on this peak
        :return: the number of variants at this position (If the position is unavailable, return None)
        """
        return self.var_pos_to_cnt.get(var_pos)

    def get_variants(self):
        """
        :return: the set of the repr strings of the variants on this peak
        """
        var_repr_set = set()

        for var_pos in self.var_pos_to_var_set:
            var_repr_set.update(self.var_pos_to_var_set[var_pos])

        return var_repr_set

    def get_variants_at_pos(self, var_pos):
        """
        :param var_pos: a variant position (0-based) on this peak
        :return: a set of repr strings of variants at this position (If the position is unavailable, return None)
        """
        return self.var_pos_to_var_set.get(var_pos)

    def get_anno_val_at_pos(self, var_pos):
        """
        :param var_pos: a variant position (0-based) inside the bed
        :return: a genic region value of the position
        """
        anno_val = self.var_pos_to_anno_val.get(var_pos)
        assert anno_val is not None

        return anno_val

    def get_genes_at_pos(self, var_pos):
        """
        :param var_pos: a variant position (0-based) inside the bed
        :return: a dictionary (mapping route: gene symbol -> ID -> genic region)
        """
        return self.var_pos_to_genes[var_pos]

    def get_anno_vals(self):
        """
        Return the annotation values of this bed (each element corresponds to each nucleotide on this bed).
        If it is not set, return an empty list.
        """
        return self.anno_vals

    def get_genic_region_to_size(self, only_repr=False):
        """
        :param only_repr: if True, then return repr_genic_region_to_size
        :return: a dictionary that documents the size of each genic region on this bed.
                 (key: a genic region, value: a size of the genic region (integer))
        """
        if only_repr:
            return self.repr_genic_region_to_size
        else:
            return self.genic_region_to_size

    def get_genic_region_to_var_cnt(self, only_repr=False):
        """
        :param only_repr: if True, then return repr_genic_region_to_var_cnt
        :return: a dictionary that documents the number of variants of each genic region on this bed.
                 (key: a genic region, value: the number of variants on the genic region (integer))
        """
        if only_repr:
            return self.repr_genic_region_to_var_cnt
        else:
            return self.genic_region_to_var_cnt
