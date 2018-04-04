"""
Module for the genome sequence
1. read partial sequence
2. make reverse complementary sequence for the input sequence
"""

from lab.genome.settings import GENOME_FILENAME
from lab.genome.fasta import Fasta

_GENOME = Fasta(GENOME_FILENAME)  # singleton


def get_seq(chrom, start, end, upper=True):
    """
    :param chrom: chromosome ID (e.g. chr1, chr2, ...)
    :param start: a start position on the chromosome (0-based)
    :param end: an end position on the chromosome
    :param upper: if True, return a sequence which characters are all capitalized.
    :return: the corresponding genome sequence
    """
    start = int(start)
    end = int(end)
    seq = _GENOME.fetch_seq(chrom, start, end)

    if upper:
        return seq.upper()
    else:
        return seq

# END: get_seq


def reverse_complement(seq):
    """
    :param seq: a string that consists of deoxy-nucleotides (A, T, G, C, ...)
    :return: reverse complementary sequence of that input
    """
    return _GENOME.reverse_complement(seq)

# END: reverse_complement
