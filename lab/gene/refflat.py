from lab.utils import eprint
from lab.genome.genome import get_seq, reverse_complement
from lab.gene.anno import genic_region_list

__all__ = ['RefFlat']


class RefFlat:
    """ represents each entry of the RefFlat file """
    def __init__(self):
        # default information of RefFlat
        self.symbol = None
        self.id = None  # NR: ncRNA, NM: mRNA
        self.chrom = None
        self.strand = None
        self.tx_start = 0
        self.tx_end = 0
        self.cds_start = 0
        self.cds_end = 0
        self.exon_cnt = 0
        self.exon_starts = []
        self.exon_ends = []

        # information gotten after parsing the default information
        self.exons_size = 0
        self.seq_5utr = ''
        self.seq_3utr = ''
        self.seq_orf = ''

    def __str__(self):
        """
        represents the object as the RefFlat entry
        """
        exon_starts = ""
        exon_ends = ""

        for i in range(self.exon_cnt):
            exon_starts += "%s," % self.exon_starts[i]
            exon_ends += "%s," % self.exon_ends[i]

        return "%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%s\t%s" % (self.symbol, self.id, self.chrom, self.strand,
                                                               self.tx_start, self.tx_end, self.cds_start, self.cds_end,
                                                               self.exon_cnt, exon_starts, exon_ends)

    def parse_refflat_line(self, refflat_line):
        refflat_fields = refflat_line.strip().split('\t')

        self.symbol = refflat_fields[0]
        self.id = refflat_fields[1]
        self.chrom = refflat_fields[2]
        self.strand = refflat_fields[3]
        self.tx_start = int(refflat_fields[4])
        self.tx_end = int(refflat_fields[5])

        if self.tx_start > self.tx_end:
            eprint("[ERROR] in %s: Tx end point is ahead of tx start." % self.symbol)
            return

        self.cds_start = int(refflat_fields[6])
        self.cds_end = int(refflat_fields[7])

        if self.cds_start > self.cds_end:
            eprint("[ERROR] in %s: CDS end point is ahead of CDS start." % self.symbol)
            return

        self.exon_cnt = int(refflat_fields[8])
        self.exon_starts = [int(exon_start) for exon_start in refflat_fields[9].strip(',').split(',')]
        self.exon_starts.sort()

        if self.tx_start != self.exon_starts[0]:
            eprint("[ERROR] in %s: Tx start point is different with exon start point." % self.symbol)
            return

        if self.exon_cnt != len(self.exon_starts):
            eprint("[ERROR] in %s: Exon count is different with the number of exon starts." % self.symbol)
            return

        self.exon_ends = [int(exon_end) for exon_end in refflat_fields[10].strip(',').split(',')]
        self.exon_ends.sort()

        if self.tx_end != self.exon_ends[-1]:
            eprint("[ERROR] in %s: Tx end point is different with exon end point." % self.symbol)
            return

        if self.exon_cnt != len(self.exon_ends):
            eprint("[ERROR] in %s: Exon count is different with the number of exon ends." % self.symbol)
            return

        valid_exon = self._check_exon_overlap()

        if not valid_exon:
            return

    def parse_exon_seq(self):
        """
        By parsing the default information from RefFlat,
        get the total size of the exons and the sequences of ORF and UTR.
        :return: a boolean value (if False, the default information is invalid)
        """
        seq_exons = ""

        for i in range(self.exon_cnt):
            one_exon_size = self.exon_ends[i] - self.exon_starts[i]

            if one_exon_size < 0:
                eprint("[ERROR] in {0}: exon{1} end point is ahead of exon{1} start.".format(self.symbol, i + 1))
                return False

            seq_exon = get_seq(self.chrom, self.exon_starts[i], self.exon_ends[i])

            seq_exons += seq_exon
            self.exons_size += one_exon_size

        cds_start_offset = 0  # inclusive

        for i in range(self.exon_cnt):
            if self.cds_start < self.exon_ends[i]:
                cds_start_offset += (self.cds_start - self.exon_starts[i])
                break
            else:
                cds_start_offset += (self.exon_ends[i] - self.exon_starts[i])

        cds_end_offset = self.exons_size  # exclusive

        for i in range(self.exon_cnt - 1, -1, -1):  # reverse for loop
            if self.cds_end >= self.exon_starts[i]:
                cds_end_offset -= (self.exon_ends[i] - self.cds_end)
                break
            else:
                cds_end_offset -= (self.exon_ends[i] - self.exon_starts[i])

        if self.strand == '+':
            self.seq_5utr = seq_exons[0:cds_start_offset]
            self.seq_orf = seq_exons[cds_start_offset:cds_end_offset]
            self.seq_3utr = seq_exons[cds_end_offset:self.exons_size]

        elif self.strand == '-':
            self.seq_5utr = reverse_complement(seq_exons[cds_end_offset:self.exons_size])
            self.seq_orf = reverse_complement(seq_exons[cds_start_offset:cds_end_offset])
            self.seq_3utr = reverse_complement(seq_exons[0:cds_start_offset])

        else:
            eprint("[ERROR] Invalid strand %s" % self.strand)
            return False

        return True

    def has_wrong_orf(self):
        if not self._has_start_codon():
            return True

        if not self._has_stop_codon():
            return True

        if len(self.seq_orf) % 3 != 0:
            return True

        if self._has_internal_stop_codon():
            return True

        return False

    def is_nmd_candidate(self):
        """
        NMD: Nonsense Mediated Decay
        NMD candidate: The boundary between orf and 3'utr is ahead of the last exon junction + 50nt upstream.
        """
        if self.strand == '+':
            last_exon_size = self.exon_ends[-1] - self.exon_starts[-1]
        else:  # '-'
            last_exon_size = self.exon_ends[0] - self.exon_starts[0]

        size_3utr = len(self.seq_3utr)

        if size_3utr > (last_exon_size + 50):
            return True
        else:
            return False

    def find_genic_region(self, pos):
        """
        :param pos: 0-based position
        :return: a genic region. if intergenic, return None
        """
        if self.tx_start - 300 <= pos < self.tx_start:
            return 'promoter'
        elif self.tx_start <= pos < self.tx_end:
            idx = -1  # even number (include 0): exonic, odd number: intronic

            for i in range(self.exon_cnt):
                if pos < self.exon_starts[i]:
                    break
                elif self.exon_starts[i] <= pos < self.exon_ends[i]:
                    idx = 2 * i
                else:
                    idx = 2 * i + 1

            assert idx != -1 and idx < 2 * self.exon_cnt - 1  # not the outside of the gene

            is_mrna = (self.id[:2] == 'NM')

            if idx % 2 == 1:
                if is_mrna:
                    exon_end_idx = int(idx / 2)
                    intron_start = self.exon_ends[exon_end_idx]
                    intron_end = self.exon_starts[exon_end_idx + 1]

                    space_from_jct = min((pos - intron_start + 1), (intron_end - pos))

                    if space_from_jct <= 30:
                        return 'SS'
                    else:
                        return 'intronic'
                else:
                    return 'ncRNA_intronic'
            else:
                if is_mrna:
                    if pos < self.cds_start:
                        return '5UTR' if self.strand == '+' else '3UTR'
                    elif pos >= self.cds_end:
                        return '3UTR' if self.strand == '+' else '5UTR'
                    else:
                        return 'ORF'
                else:
                    return 'ncRNA_exonic'
        else:
            eprint("[ERROR] The position %d is not in the gene %s(%s)." % (pos, self.symbol, self.id))
            return None

    def is_non_synonymous(self, pos, ref_nuc, alt_nuc):
        """
        Determine whether the variant on this gene is non-synonymous or not.
        (It is assumed that the variant is on a same chromosome with this gene)
        :param pos: 0-based position on the genome
        :param ref_nuc: a reference nucleotide of this variant
        :param alt_nuc: a alternative nucleotide of this variants
        :return: a boolean value. If True, this variant is non-synonymous.
        """
        genic_region = self.find_genic_region(pos)

        if genic_region is not None and genic_region == 'ORF':
            exon_idx = -1

            for i in range(self.exon_cnt):
                if self.exon_starts[i] <= pos < self.exon_ends[i]:
                    exon_idx = i
                    break

            assert exon_idx != -1

            # find the relative position of the variant on mRNA
            rel_pos = pos

            for j in range(exon_idx):
                intron_size = self.exon_starts[j + 1] - self.exon_ends[j]
                rel_pos -= intron_size

            # find the relative position of the variant on CDS of mRNA
            if self.strand == '+':
                cds_rel_pos = rel_pos - len(self.seq_5utr)
            else:  # self.strand == '-'
                cds_rel_pos = rel_pos - len(self.seq_3utr)
                cds_rel_pos = len(self.seq_orf) - cds_rel_pos -1  # reverse complement

            # TODO: Implementation of codes for checking the modification of AA
        else:
            return False

    def get_genic_region_dist(self, start, end):
        """
        Return the dictionary that contains the size of each genic region in the region (start, end)
        :param start: an integer (0-based)
        :param end: an integer
        :return: a dictionary (key: genic region, value: size)
        """
        if self.tx_start <= start < end <= self.tx_end:
            # even number (include 0): exonic, odd number: intronic
            start_idx = -1
            end_idx = -1

            for i in range(self.exon_cnt):
                if start < self.exon_starts[i]:
                    break
                elif self.exon_starts[i] <= start < self.exon_ends[i]:
                    start_idx = 2 * i
                else:
                    start_idx = 2 * i + 1

            for i in range(self.exon_cnt):
                if end <= self.exon_starts[i]:
                    break
                elif self.exon_starts[i] < end <= self.exon_ends[i]:
                    end_idx = 2 * i
                else:
                    end_idx = 2 * i + 1

            assert start_idx != -1 and start_idx < 2 * self.exon_cnt - 1
            assert end_idx != -1 and end_idx < 2 * self.exon_cnt - 1

            is_mrna = (self.id.startswith('NM'))
            is_top_strand = (self.strand == '+')

            # initialization
            genic_regions = genic_region_list()
            region_to_size = {genic_region: 0 for genic_region in genic_regions}

            for i in range(start_idx, end_idx + 1):
                if i % 2 == 0:  # exon
                    exon_idx = int(i / 2)
                    start_pos = self.exon_starts[exon_idx]
                    end_pos = self.exon_ends[exon_idx]

                    if start_pos < start:
                        start_pos = start

                    if end_pos > end:
                        end_pos = end

                    if is_mrna:
                        left_utr = 0
                        orf = 0
                        right_utr = 0

                        if end_pos <= self.cds_start:
                            left_utr = end_pos - start_pos
                        elif start_pos >= self.cds_end:
                            right_utr = end_pos - start_pos
                        else:
                            if start_pos < self.cds_start:
                                left_utr = self.cds_start - start_pos
                                start_pos = self.cds_start

                            if end_pos > self.cds_end:
                                right_utr = end_pos - self.cds_end
                                end_pos = self.cds_end

                            orf = end_pos - start_pos

                        if is_top_strand:
                            region_to_size['5UTR'] += left_utr
                            region_to_size['3UTR'] += right_utr
                        else:
                            region_to_size['5UTR'] += right_utr
                            region_to_size['3UTR'] += left_utr

                        region_to_size['ORF'] += orf
                    else:
                        region_to_size['ncRNA_exonic'] += (end_pos - start_pos)

                else:  # intron
                    intron_idx = int(i / 2)
                    intron_start = self.exon_ends[intron_idx]
                    intron_end = self.exon_starts[intron_idx + 1]

                    start_pos = intron_start
                    end_pos = intron_end

                    if start_pos < start:
                        start_pos = start

                    if end_pos > end:
                        end_pos = end

                    intron_size = (end_pos - start_pos)

                    if is_mrna:
                        region_to_size['intronic'] += intron_size

                        left_ss_end = intron_start + 30
                        right_ss_start = intron_end - 30

                        left_ss_size = left_ss_end - start_pos
                        right_ss_size = end_pos - right_ss_start

                        if left_ss_size < 0:
                            left_ss_size = 0

                        if right_ss_size < 0:
                            right_ss_size = 0

                        ss_size = left_ss_size + right_ss_size

                        if ss_size > intron_size:
                            ss_size = intron_size

                        region_to_size['SS'] += ss_size
                        region_to_size['intronic'] -= ss_size

                    else:
                        region_to_size['ncRNA_intronic'] += intron_size

            return region_to_size

        else:
            return None

    def _check_exon_overlap(self):
        exon_positions = []

        for i in range(self.exon_cnt):
            exon_positions.append(self.exon_starts[i])
            exon_positions.append(self.exon_ends[i])

        num_iter = 2 * self.exon_cnt - 1

        for i in range(num_iter):
            if exon_positions[i] > exon_positions[i + 1]:
                eprint("[ERROR] %s has a exon overlap." % self.symbol)
                return False

        return True

    def _has_start_codon(self):
        if self.seq_orf[:3] == 'ATG':
            return True
        else:
            return False

    def _has_stop_codon(self):
        stop_codons = ['TAA', 'TAG', 'TGA']

        if self.seq_orf[-3:] in stop_codons:
            return True
        else:
            return False

    def _has_internal_stop_codon(self):
        stop_codons = ['TAA', 'TAG', 'TGA']
        orf_size = len(self.seq_orf)

        for i in range(0, orf_size - 3, 3):
            codon = self.seq_orf[i:i + 3]

            if codon in stop_codons:
                return True

        return False
