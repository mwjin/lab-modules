from genome.reader import read_partial_seq, reverse_complement

class RefFlat:
    """ save the information from refFlat file """
    def __init__(self):
        self.symbol = None
        self.id = None  # NR: ncRNA, NM: mRNA
        self.chrID = None
        self.strand = None
        self.tx_start = 0
        self.tx_end = 0
        self.cds_start = 0
        self.cds_end = 0
        self.exon_cnt = 0
        self.exon_starts = []
        self.exon_ends = []
        self.exons_size = 0
        self.seq_5UTR = 0
        self.seq_3UTR = 0
        self.seq_ORF = 0

    def __str__(self):
        exon_starts = ""
        exon_ends = ""

        for i in range(self.exon_cnt):
            exon_starts += "%s," % self.exon_starts[i]
            exon_ends += "%s," % self.exon_ends[i]

        return "%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%s\t%s" % (self.symbol, self.id, self.chrID, self.strand,
                                                               self.tx_start, self.tx_end, self.cds_start, self.cds_end,
                                                               self.exon_cnt, exon_starts, exon_ends)

    def parse_refFlat_line(self, refFlat_line):
        refFlat_fields = refFlat_line.strip().split('\t')

        self.symbol = refFlat_fields[0]
        self.id = refFlat_fields[1]
        self.chrID = refFlat_fields[2]
        self.strand = refFlat_fields[3]
        self.tx_start = int(refFlat_fields[4])
        self.tx_end = int(refFlat_fields[5])

        if self.tx_start > self.tx_end:
            print("Error in %s: tx end point is ahead of tx start." % (self.symbol))
            return

        self.cds_start = int(refFlat_fields[6])
        self.cds_end = int(refFlat_fields[7])

        if self.cds_start > self.cds_end:
            print("Error in %s: cds end point is ahead of cds start." % (self.symbol))
            return

        self.exon_cnt = int(refFlat_fields[8])
        self.exon_starts = [ int(exon_start) for exon_start in refFlat_fields[9].strip(',').split(',') ]
        self.exon_starts.sort()

        if self.tx_start != self.exon_starts[0]:
            print("Error in %s: Tx start point is different with exon start point." % (self.symbol))
            return

        if self.exon_cnt != len(self.exon_starts):
            print("Error in %s: Exon count is different with the number of exon starts." % (self.symbol))
            return

        self.exon_ends = [ int(exon_end) for exon_end in refFlat_fields[10].strip(',').split(',') ]
        self.exon_ends.sort()

        if self.tx_end != self.exon_ends[-1]:
            print("Error in %s: Tx end point is different with exon end point." % (self.symbol))
            return

        if self.exon_cnt != len(self.exon_ends):
            print("Error in %s: Exon count is different with the number of exon ends." % (self.symbol))
            return

        valid_exon = self._check_exon_overlap()

        if not valid_exon:
            return

    def _check_exon_overlap(self):
        exon_positions = []

        for i in range(self.exon_cnt):
            exon_positions.append(self.exon_starts[i])
            exon_positions.append(self.exon_ends[i])

        num_iter = 2 * self.exon_cnt - 1

        for i in range(num_iter):
            if exon_positions[i] > exon_positions[i + 1]:
                print("Error: %s has a exon overlap." % self.symbol)
                return False

        return True

    def parse_exon_seq(self):
        seq_exons = ""

        for i in range(self.exon_cnt):
            one_exon_size = self.exon_ends[i] - self.exon_starts[i]

            if one_exon_size < 0:
                print("Error in {0}: exon{1} end point is ahead of exon{1} start.".format(self.symbol, i + 1))
                return False

            seq_exon = read_partial_seq(self.chrID, self.exon_starts[i], self.exon_ends[i])

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

        for i in range(self.exon_cnt - 1, -1, -1):
            if self.cds_end >= self.exon_starts[i]:
                cds_end_offset -= (self.exon_ends[i] - self.cds_end)
                break
            else:
                cds_end_offset -= (self.exon_ends[i] - self.exon_starts[i])

        if self.strand == '+':
            self.seq_5UTR = seq_exons[0 : cds_start_offset]
            self.seq_ORF = seq_exons[cds_start_offset : cds_end_offset]
            self.seq_3UTR = seq_exons[cds_end_offset : self.exons_size]

        elif self.strand == '-':
            self.seq_5UTR = reverse_complement(seq_exons[cds_end_offset : self.exons_size])
            self.seq_ORF = reverse_complement(seq_exons[cds_start_offset : cds_end_offset])
            self.seq_3UTR = reverse_complement(seq_exons[0 : cds_start_offset])

        else:
            print("Error: invalid strand %s" % self.strand)
            return False

        return True

    def has_wrong_ORF(self):
        if not self._has_start_codon():
            return True

        if not self._has_stop_codon():
            return True

        if len(self.seq_ORF) % 3 != 0:
            return True

        if self._has_internal_stop_codon():
            return True

        return False

    def _has_start_codon(self):
        if self.seq_ORF[:3] == 'ATG':
            return True
        else:
            return False

    def _has_stop_codon(self):
        stop_codons = ['TAA', 'TAG', 'TGA']

        if self.seq_ORF[-3:] in stop_codons:
            return True
        else:
            return False

    def _has_internal_stop_codon(self):
        stop_codons = ['TAA', 'TAG', 'TGA']
        ORF_size = len(self.seq_ORF)

        for i in range(0, ORF_size - 3, 3):
            codon = self.seq_ORF[i : i + 3]

            if codon in stop_codons:
                return True

        return False

    def is_NMD_candidate(self):
        """ The boundary between ORF and 3'UTR is ahead of the last exon junction + 50nt upstream. """
        if self.strand == '+':
            last_exon_size = self.exon_ends[-1] - self.exon_starts[-1]
        else:  # '-'
            last_exon_size = self.exon_ends[0] - self.exon_starts[0]

        size_3UTR = len(self.seq_3UTR)

        if size_3UTR > (last_exon_size + 50):
            return True
        else:
            return False

    def find_genic_region(self, pos):
        """
        :param pos: should be 0-based.
        :return: genic region (intronic, 5UTR, exonic(ORF), 3UTR). if intergenic, return None
        """
        exon_positions = []

        for i in range(self.exon_cnt):
            exon_positions.append(self.exon_starts[i])
            exon_positions.append(self.exon_ends[i])

        index = 0

        for i in range(self.exon_cnt * 2):
            if exon_positions[i] > pos:
                index = i
                break

        is_mRNA = (self.id[:2] == 'NM')

        if index == 0:
            print("Error: position %d is not in the rep-genes %s." % (pos, self.symbol))
            return None

        elif index % 2 == 0:
            return 'intronic' if is_mRNA else 'ncRNA_intronic'

        else:
            if is_mRNA:
                if pos < self.cds_start:
                    return '5UTR' if self.strand == '+' else '3UTR'
                elif pos >= self.cds_end:
                    return '3UTR' if self.strand == '+' else '5UTR'
                else:
                    return 'ORF'
            else:
                return 'ncRNA_exonic'
