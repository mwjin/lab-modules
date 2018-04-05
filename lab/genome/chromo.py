"""
Modules for the functions associated with the chromosome
1. get the size of the chromosome
2. bin the length of the chromosome
"""
from lab.genome.settings import GENOME_FILENAME


_CHR_TO_SIZE = {}  # chromosome ID to its size

with open('%s.fai' % GENOME_FILENAME, 'r') as GENOME_IDX_FILE:
    for line in GENOME_FILENAME.readlines():
        """
        Field 1: chromosome ID e.g. chr1
        Field 2: size of each chromosome without newline characters
        Field 3: size of each chromosome with newline characters
        Field 4: length of lines in each chromosome without newline characters
        Field 5: length of lines in each chromosome with newline characters
        """
        fields = line.split('\t')
        _CHR_TO_SIZE[fields[0]] = int(fields[1])


def get_chr_size(chrom):
    """
    :param chrom: the chromosome ID (e.g. chr1, chr2, chr3, ..., chrX, chrY)
    :return: the size of the chromosome
    """
    return _CHR_TO_SIZE[chrom]


def bin_chr(chrom, bin_size):
    """
    :param chrom: the chromosome ID (e.g. chr1, chr2, chr3, ..., chrX, chrY)
    :param bin_size: (integer) the size of each bin
    :return: the list of tuples each of which represents the range (start, end) of each bin

             e.g. if a size of a chromosome is 2550 and the bin size is 1000,
                  then the list of tuples of the chromosome will be [(0, 1000), (1000, 2000), (2000, 2550)].
    """
    chr_size = _CHR_TO_SIZE[chrom]
    bin_cnt = int(chr_size / bin_size) + 1
    chr_bins = []

    for i in range(bin_cnt):
        bin_start_idx = i * bin_size
        bin_end_idx = bin_start_idx + bin_size

        if bin_end_idx > chr_size:
            bin_end_idx = chr_size

        chr_bins.append((bin_start_idx, bin_end_idx))

    return chr_bins
