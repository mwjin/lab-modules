"""
Genome utilities
1. make reverse complementary sequence for the input sequence
2. read partial sequence
3. get the sizes of chromosomes
4. set the range of bins of chromosomes
"""


import sys

from lib_settings import GENOME_FILENAME
from genome.fasta import Fasta

GENOME = Fasta(GENOME_FILENAME)

def reverse_complement(seq):
    base_to_comp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
                    'a': 't', 'c': 'g', 'g': 'c', 't': 'a',
                    'N': 'N', '.':'.'}

    comp_seq = ''

    for base in seq:
        comp_seq += base_to_comp[base]

    return comp_seq[::-1]  # reverse
# END: reverse_complement


def read_partial_seq(chrID, start, end):
    start = int(start)
    end = int(end)
    seq = GENOME.fetch_seq(chrID, start, end)

    return seq.upper()
# END: read_partial_seq


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

        chrID = fields[0]
        chr_size = int(fields[1])

        chr_to_size[chrID] = chr_size

    return chr_to_size


def main():
    pass
# END: main function
    
if __name__ == '__main__':
    if len(sys.argv) == 1: main()
    else:
        function_name       = sys.argv[1]
        function_parameters = sys.argv[2:]
        if function_name in locals().keys(): locals()[function_name](*function_parameters)
        else: sys.exit('ERROR: function_name=%s, parameters=%s' % (function_name, function_parameters))
    #if END: len(sys.argv)
#if END: __name__
