import os
import sys
import re

re_nonchr = re.compile('[^a-zA-Z]')


class Fasta:
    def __init__(self, genome_filename):
        # V-S Check: File Existence
        if not os.path.isfile(genome_filename):
            sys.exit('Error: The genome file does not exist.')
        
        self.genome_file = open(genome_filename, 'r')
        self.chrom_list = []
        self.chrlen_list = []
        self.offset_list = []
        self.line_len_list = []
        self.line_len_with_blank_list = []

        # V-S Check: File Existence
        if not os.path.isfile('%s.fai' % genome_filename):
            sys.exit('.fai file does not exist')

        genome_idx_file = open('%s.fai' % genome_filename, 'r')

        for line in genome_idx_file:
            fields = line.strip('\n').split()  # Goes backwards, -1 skips the new line character

            self.chrom_list.append(fields[0])
            self.chrlen_list.append(int(fields[1]))
            self.offset_list.append(int(fields[2]))
            self.line_len_list.append(int(fields[3]))
            self.line_len_with_blank_list.append(int(fields[4]))
        # END: for loop 'line'

        genome_idx_file.close()
    # END: __init__

    def __del__(self):
        self.genome_file.close()

    def fetch_seq(self, chrom, start=None, end=None, strand='+'):
        assert chrom in self.chrom_list, chrom
        chr_idx = self.chrom_list.index(chrom)

        if start is None:
            start = 0

        if end is None:
            end = self.chrlen_list[chr_idx]

        try:
            assert (0 <= start) and (start < end) and (end <= self.chrlen_list[chr_idx])
        except AssertionError:
            print('Fasta fetch assertion error', chrom, start, end)
            sys.exit()

        blank_cnt = self.line_len_with_blank_list[chr_idx] - self.line_len_list[chr_idx]

        start = int(start + (start / self.line_len_list[chr_idx]) * blank_cnt)  # Start Fetch Position
        end = int(end + (end / self.line_len_list[chr_idx]) * blank_cnt)  # End Fetch Position

        self.genome_file.seek(self.offset_list[chr_idx] + start)            # Get Sequence

        seq = re.sub(re_nonchr, '', self.genome_file.read(end - start))

        if strand == '+':
            return seq
        elif strand == '-':
            return self._reverse_complement(seq)
        else:
            sys.exit('Error: invalid strand')
    # END: fetch_seq

    @staticmethod
    def _reverse_complement(seq):
        base_to_comp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
                        'a': 't', 'c': 'g', 'g': 'c', 't': 'a',
                        'N': 'N', '.': '.'}

        comp_seq = ''

        for base in seq:
            comp_seq += base_to_comp[base]

        return comp_seq[::-1]  # reverse

    # END: _reverse_complement
# END: 'Fasta' class
