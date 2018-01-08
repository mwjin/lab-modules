#
#######################################################################
# This code reads a part of chromosome sequence only we want.
# Author: Jinman Park
# Modified by Minwoo Jeong
# last update: 2017. 9. 15.
#######################################################################
#

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
