import os
import re

__all__ = ['set_genome']


class _Genome:
    def __init__(self):
        self.genome_file = None
        self.chroms = []
        self.chrom_to_size = {}
        self.chrom_to_offset = {}
        self.chrom_to_line_len = {}
        self.chrom_to_line_len_with_blk = {}

    def __del__(self):
        if self.genome_file is not None:
            self.genome_file.close()

    def parse_file(self, genome_file_path):
        """
        :param genome_file_path: a path of the genome (.fa)
                                 we suppose that this path is entered through the function 'set_genome'
        """
        self.genome_file = open(genome_file_path, 'r')

        # Parse the index file
        genome_idx_file_path = '%s.fai' % genome_file_path
        genome_idx_file = open(genome_idx_file_path, 'r')

        for line in genome_idx_file:
            fields = line.strip('\n').split()  # Goes backwards, -1 skips the new line character

            self.chroms.append(fields[0])
            self.chrom_to_size[fields[0]] = int(fields[1])
            self.chrom_to_offset[fields[0]] = int(fields[2])
            self.chrom_to_line_len[fields[0]] = int(fields[3])
            self.chrom_to_line_len_with_blk[fields[0]] = int(fields[4])

        genome_idx_file.close()

    def fetch_seq(self, chrom, start=None, end=None):
        """
        :param chrom: a chromosome ID
        :param start: a start position of the region
        :param end: an end position of the region
        :return: the sequence of the input region
        """
        if self.genome_file is None:
            raise AssertionError('ERROR: the genome data does not exist. call \'parse_file\' first.')

        if chrom not in self.chroms:
            raise ValueError('ERROR: invalid chromosome ID \'%s\'' % chrom)

        chr_size = self.chrom_to_size[chrom]
        offset = self.chrom_to_offset[chrom]
        line_len = self.chrom_to_line_len[chrom]
        line_len_with_blk = self.chrom_to_line_len_with_blk[chrom]

        if start is None:
            start = 0

        if end is None:
            end = chr_size

        if not (0 <= start < end <= chr_size):
            raise AssertionError('ERROR: invalid region %s:%s-%s' % (chrom, start, end))

        blank_cnt = line_len_with_blk - line_len

        start = int(start + (start / line_len) * blank_cnt)  # Start Fetch Position
        end = int(end + (end / line_len) * blank_cnt)  # End Fetch Position

        self.genome_file.seek(offset + start)            # Get Sequence

        re_nonchr = re.compile('[^a-zA-Z]')
        seq = re.sub(re_nonchr, '', self.genome_file.read(end - start))

        return seq

    def get_chr_size(self, chrom):
        """
        Return the size of the chromosome
        """
        return self.chrom_to_size[chrom]


_GENOME = _Genome()


def set_genome(genome_file_path):
    """
    Make _Fasta object

    It is essential to execute this function for using this module.
    'genome_file_path' must be '.fa' format and there must be index file (.fai) in the same directory.
    If the 'genome_file_path' is wrong, AssertionError occur.
    """
    # Sanity check
    if not os.path.isfile(genome_file_path):
        raise AssertionError('ERROR: the genome file \'%s\' does not exist.' % genome_file_path)

    genome_idx_file_path = '%s.fai' % genome_file_path

    if not os.path.isfile(genome_idx_file_path):
        raise AssertionError('ERROR: the genome index file \'%s.fa\' does not exist' % genome_file_path)

    _GENOME.parse_file(genome_file_path)


def reverse_complement(seq):
    base_to_comp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
                    'a': 't', 'c': 'g', 'g': 'c', 't': 'a',
                    'N': 'N', 'n': 'n', '.': '.'}

    comp_seq = ''

    for base in seq:
        comp_seq += base_to_comp[base]

    return comp_seq[::-1]  # reverse


def get_seq(chrom, start, end, strand='+', upper=True):
    """
    :param chrom: chromosome ID (e.g. chr1, chr2, ...)
    :param start: a start position on the chromosome (0-based)
    :param end: an end position on the chromosome
    :param strand: '+' if you want top strand, else '-'
    :param upper: if True, return a sequence which characters are all capitalized.
    :return: the corresponding genome sequence
    """
    start = int(start)
    end = int(end)
    seq = _GENOME.fetch_seq(chrom, start, end)

    if strand == '-':
        seq = reverse_complement(seq)
    if upper:
        return seq.upper()
    else:
        return seq


def get_chr_size(chrom):
    """
    :param chrom: the chromosome ID (e.g. chr1, chr2, chr3, ..., chrX, chrY)
    :return: the size of the chromosome
    """
    return _GENOME.get_chr_size(chrom)


def bin_chr(chrom, bin_size, overlap=0):
    """
    :param chrom: the chromosome ID (e.g. chr1, chr2, chr3, ..., chrX, chrY)
    :param bin_size: (integer) the size of each bin
    :param overlap: a size of the overlap between two tandem bins
    :return: the list of tuples each of which represents the range (start, end) of each bin

             e.g. if a size of a chromosome is 2550 and the bin size is 1000,
                  then the list of tuples of the chromosome will be [(0, 1000), (1000, 2000), (2000, 2550)].
    """
    assert bin_size > overlap

    chr_size = _GENOME.get_chr_size(chrom)
    chr_bins = []

    bin_start = 0
    bin_end = bin_size

    while bin_start < chr_size:
        chr_bins.append((bin_start, bin_end))

        bin_start += (bin_size - overlap)
        bin_end = bin_start + bin_size

        if bin_end > chr_size:
            bin_end = chr_size

    return chr_bins
