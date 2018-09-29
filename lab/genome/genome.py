"""
Useful functions associated with human genome
1. Get a reverse complement sequence
2. Get a sequence of the arbitrary region
3. Get a size of the chromosome
4. Bin the chromosome

For using them except 1, we must call the function 'set_genome' first.
"""
import os
import re

__all__ = ['set_genome', 'get_seq', 'get_chr_size', 'bin_chrom',
           'reverse_complement', 'complementary_base', 'get_amino_acid']


class _Genome:
    """
    This class represents one genome.
    The methods in this class can be accessed by the outer functions in the same script.
    """
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
        if self.genome_file is not None and self.genome_file.name == genome_file_path:
            return

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
        :param start: a start position of the region (0-based)
        :param end: an end position of the region
        :return: the sequence of the input region
        """
        if self.genome_file is None:
            raise AssertionError('[ERROR] The genome file does not be entered. call \'parse_file\' first.')

        if chrom not in self.chroms:
            raise ValueError('[ERROR] Invalid chromosome ID \'%s\'' % chrom)

        chr_size = self.chrom_to_size[chrom]
        offset = self.chrom_to_offset[chrom]
        line_len = self.chrom_to_line_len[chrom]
        line_len_with_blk = self.chrom_to_line_len_with_blk[chrom]

        if start is None:
            start = 0

        if end is None:
            end = chr_size

        if not (0 <= start < end <= chr_size):
            raise AssertionError('[ERROR] Invalid region %s:%s-%s' % (chrom, start, end))

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


# Functions to access the methods of the class '_Genome'
def set_genome(genome_file_path):
    """
    Make _Genome object

    It is essential to execute this function for using this module.
    'genome_file_path' must be '.fa' format and there must be index file (.fai) in the same directory.
    If the 'genome_file_path' is wrong, AssertionError occur.
    """
    # Sanity check
    if not os.path.isfile(genome_file_path):
        raise AssertionError('[ERROR] the genome file \'%s\' does not exist.' % genome_file_path)

    genome_idx_file_path = '%s.fai' % genome_file_path

    if not os.path.isfile(genome_idx_file_path):
        raise AssertionError('[ERROR] the genome index file \'%s.fa\' does not exist' % genome_file_path)

    _GENOME.parse_file(genome_file_path)


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


def bin_chrom(chrom, bin_size, overlap=0):
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


# Functions not to access the class
_COMPLEMENTARY_BASE = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
                       'a': 't', 'c': 'g', 'g': 'c', 't': 'a',
                       'N': 'N', 'n': 'n', '.': '.'}

# For convenience, the key value is a DNA triplet.
# Types of amino acids
# F: Phenylalanine, L: Leucine, S: Serine, Y: Tyrosine, C: Cystein, W: Tryptophan, P: Proline, H: Histidine,
# Q: Glutamine, R: Arginine, I: Isoleucine, M: Methionine, T: Threonine, N: Asparagine, K: Lysine, V: Valine,
# A: Alanine, D: Aspartic acid, E: Glutamic acid, G: Glycine
_CODON_TABLE = {'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
                'TAT': 'Y', 'TAC': 'Y', 'TAA': 'STOP', 'TAG': 'STOP', 'TGT': 'C', 'TGC': 'C', 'TGA': 'STOP', 'TGG': 'W',
                'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
                'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
                'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
                'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
                'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
                'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'}


def reverse_complement(seq):
    """
    :param seq: a sequence consists of DNA nucleotides
    :return: the reverse complement sequence of that sequence
    """
    comp_seq = ''

    for base in seq:
        comp_seq += _COMPLEMENTARY_BASE[base]

    return comp_seq[::-1]  # reverse


def complementary_base(base):
    """
    :param base: a base of a nucleotide
    :return: the complementary base of that base. If the base is wrong, return None
    """
    return _COMPLEMENTARY_BASE.get(base)


def get_amino_acid(codon):
    """
    :param codon: a DNA triplet
    :return: an amino acid represented as a character or 'STOP'. If the codon is wrong, return None.
    """
    return _CODON_TABLE.get(codon)
