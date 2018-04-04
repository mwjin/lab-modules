"""
Genome utilities
1. make reverse complementary sequence for the input sequence
2. read partial sequence
3. get the sizes of chromosomes
4. set the range of bins of chromosomes
"""


import sys

from lab.genome.settings import GENOME_FILENAME
from lab.genome.fasta import Fasta

_GENOME = Fasta(GENOME_FILENAME)  # singleton


def reverse_complement(seq):
    return _GENOME.reverse_complement(seq)

# END: reverse_complement


def get_seq(chrom, start, end, upper=True):
    start = int(start)
    end = int(end)
    seq = _GENOME.fetch_seq(chrom, start, end)

    if upper:
        return seq.upper()
    else:
        return seq
# END: get_seq


def get_chr_sizes():
    """
    :return: a dictionary that maps chromosome ID to the size of this chromosome
    """
    chr_to_size = {}

    genome_idx_file = open('%s.fai' % GENOME_FILENAME, 'r')

    for line in genome_idx_file.readlines():
        """
        Column 1: chromosome ID e.g. chr1
        Column 2: size of each chromosome without newline characters
        Column 3: size of each chromosome with newline characters
        Column 4: length of lines in each chromosome without newline characters
        Column 5: length of lines in each chromosome with newline characters
        """
        fields = line.split('\t')

        chrom = fields[0]
        chr_size = int(fields[1])

        chr_to_size[chrom] = chr_size

    return chr_to_size


def get_chr_bin_ranges(bin_size):
    """
    :param bin_size: (integer) the size of each bin you want
    :return: a dictionary that maps chromosome ID to the list of tuples
             that represent the ranges of bins in each chromosome

             e.g. if a size of a chromosome is 2550 and the bin size is 1000,
                  then the list of tuples of the chromosome will be [(0, 1000), (1000, 2000), (2000, 2550)].
    """
    chr_to_size = get_chr_sizes()
    chr_to_bin_ranges = {}

    for chrom in chr_to_size:
        assert chrom not in chr_to_bin_ranges
        chr_to_bin_ranges[chrom] = []

        chr_size = chr_to_size[chrom]
        bin_cnt = int(chr_to_size[chrom] / bin_size) + 1

        for i in range(bin_cnt):
            bin_start_idx = i * bin_size
            bin_end_idx = bin_start_idx + bin_size

            if bin_end_idx > chr_size:
                bin_end_idx = chr_size

            chr_to_bin_ranges[chrom].append((bin_start_idx, bin_end_idx))

    return chr_to_bin_ranges


def main():
    pass
# END: main function

if __name__ == '__main__':
    if len(sys.argv) == 1:
        main()
    else:
        function_name = sys.argv[1]
        function_parameters = sys.argv[2:]
        if function_name in locals().keys():
            locals()[function_name](*function_parameters)
        else:
            sys.exit('ERROR: function_name=%s, parameters=%s' % (function_name, function_parameters))
