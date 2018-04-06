from lab.utils import eprint
from lab.genome.seq import get_seq, reverse_complement


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
            eprint("Error in %s: tx end point is ahead of tx start." % self.symbol)
            return

        self.cds_start = int(refflat_fields[6])
        self.cds_end = int(refflat_fields[7])

        if self.cds_start > self.cds_end:
            eprint("Error in %s: cds end point is ahead of cds start." % self.symbol)
            return

        self.exon_cnt = int(refflat_fields[8])
        self.exon_starts = [int(exon_start) for exon_start in refflat_fields[9].strip(',').split(',')]
        self.exon_starts.sort()

        if self.tx_start != self.exon_starts[0]:
            eprint("Error in %s: Tx start point is different with exon start point." % self.symbol)
            return

        if self.exon_cnt != len(self.exon_starts):
            eprint("Error in %s: Exon count is different with the number of exon starts." % self.symbol)
            return

        self.exon_ends = [int(exon_end) for exon_end in refflat_fields[10].strip(',').split(',')]
        self.exon_ends.sort()

        if self.tx_end != self.exon_ends[-1]:
            eprint("Error in %s: Tx end point is different with exon end point." % self.symbol)
            return

        if self.exon_cnt != len(self.exon_ends):
            eprint("Error in %s: Exon count is different with the number of exon ends." % self.symbol)
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
                eprint("Error: %s has a exon overlap." % self.symbol)
                return False

        return True

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
                eprint("Error in {0}: exon{1} end point is ahead of exon{1} start.".format(self.symbol, i + 1))
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
            eprint("Error: invalid strand %s" % self.strand)
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
        :return: genic region (intronic, 5utr, exonic(orf), 3utr). if intergenic, return None
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

        is_mrna = (self.id[:2] == 'NM')

        if index == 0:
            eprint("Error: position %d is not in the rep-genes %s." % (pos, self.symbol))
            return None

        elif index % 2 == 0:
            return 'intronic' if is_mrna else 'ncRNA_intronic'

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
