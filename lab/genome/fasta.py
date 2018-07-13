import os
import re

__all__ = ['Fasta']


class Fasta:
    def __init__(self, genome_file_path):
        """
        :param genome_file_path: a path of the genome (.fa)
        """
        if genome_file_path is None:
            raise AssertionError('ERROR: the argument for the path of the genome file is necessary.')

        if not os.path.isfile(genome_file_path):
            raise IOError('ERROR: the genome file \'%s\' does not exist.' % genome_file_path)
        
        self.genome_file = open(genome_file_path, 'r')
        self.chroms = []
        self.chr_lens = []
        self.offsets = []
        self.line_lens = []
        self.line_lens_with_blank = []

        genome_idx_file_path = '%s.fai' % genome_file_path

        if not os.path.isfile(genome_idx_file_path):
            raise IOError('ERROR: the genome index file \'%s.fa\' does not exist' % genome_file_path)

        # Parse the index file
        genome_idx_file = open('%s.fai' % genome_file_path, 'r')

        for line in genome_idx_file:
            fields = line.strip('\n').split()  # Goes backwards, -1 skips the new line character

            self.chroms.append(fields[0])
            self.chr_lens.append(int(fields[1]))
            self.offsets.append(int(fields[2]))
            self.line_lens.append(int(fields[3]))
            self.line_lens_with_blank.append(int(fields[4]))

        genome_idx_file.close()

    def __del__(self):
        self.genome_file.close()

    def fetch_seq(self, chrom, start=None, end=None, strand='+'):
        """
        :param chrom: a chromosome ID
        :param start: a start position of the region
        :param end: an end position of the region
        :param strand: '+' if the region is on the top strand, else '-'
        :return: the sequence of the input region
        """
        if chrom not in self.chroms:
            raise ValueError('ERROR: invalid chromosome ID \'%s\'' % chrom)

        chr_idx = self.chroms.index(chrom)

        if start is None:
            start = 0

        if end is None:
            end = self.chr_lens[chr_idx]

        if not ((0 <= start) and (start < end) and (end <= self.chr_lens[chr_idx])):
            raise AssertionError('ERROR: invalid region %s:%s-%s' % (chrom, start, end))

        blank_cnt = self.line_lens_with_blank[chr_idx] - self.line_lens[chr_idx]

        start = int(start + (start / self.line_lens[chr_idx]) * blank_cnt)  # Start Fetch Position
        end = int(end + (end / self.line_lens[chr_idx]) * blank_cnt)  # End Fetch Position

        self.genome_file.seek(self.offsets[chr_idx] + start)            # Get Sequence

        re_nonchr = re.compile('[^a-zA-Z]')
        seq = re.sub(re_nonchr, '', self.genome_file.read(end - start))

        if strand == '+':
            return seq
        elif strand == '-':
            return self.reverse_complement(seq)
        else:
            raise ValueError('ERROR: invalid strand %s' % strand)

    @staticmethod
    def reverse_complement(seq):
        base_to_comp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
                        'a': 't', 'c': 'g', 'g': 'c', 't': 'a',
                        'N': 'N', 'n': 'n', '.': '.'}

        comp_seq = ''

        for base in seq:
            comp_seq += base_to_comp[base]

        return comp_seq[::-1]  # reverse
