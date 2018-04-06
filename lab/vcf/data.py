import sys
import os
import re
import gzip

from lab.utils import eprint
from lab.gene.anno import get_genic_region_val, parse_genic_region_val

__all__ = ['VCFData']


class VCFData:
    """ The object of this class represents one entry in VCF files (only consider SNV). """
    def __init__(self):
        self.chrom = ''
        self.pos = 0  # 1-base
        self.dbSNP_id = ''
        self.ref_nuc = ''
        self.alt_nuc = ''
        self.qual = 0.0
        self.filter = ''
        self.info = ''
        self.format = ''
        self.samples = []

        """ info from refFlat (in-house rep-isoforms) """
        self._genic_region_dict = {'+': {}, '-': {}}  # Mapping route: strand -> gene symbol -> ID -> genic region
        self._strand_to_region_val = {'+': 0, '-': 0}  # value: genic region value (see gene.utils)

    # END: __init__

    def get_assigned_strand(self):
        """
        * Notice: there is no strand for variants actually. However, to determine the transcript (+ or -)
                  this variant was more likely to be assigned, we set this virtual strand of this variant.

        :return: a strand where this variant was more likely to be assigned.
        """

        if self._strand_to_region_val['+'] >= self._strand_to_region_val['-']:
            return '+'
        else:
            return '-'
    # END: the function 'get_assigned_strand'

    def get_strand_region_val(self, strand):
        """
        :param strand: '+' or '-'
        :return: the genic region value on the input strand
        """
        return self._strand_to_region_val[strand]

    # END: the function 'get_strand_region_val'

    @staticmethod
    def parse_vcf_file(vcf_filename):

        if not os.path.isfile(vcf_filename):
            raise IOError('%s does not exist' % vcf_filename)

        # regex for VCF entries
        ptn_auto = re.compile('^[0-9]{1,2}$')
        ptn_sex = re.compile('^[XY]$')
        ptn_pos = re.compile('^[0-9]+$')

        var_list = []
        invalid_chr_to_cnt = {}

        if vcf_filename.endswith('.gz'):
            vcf_file = gzip.open(vcf_filename, 'rt')
        else:
            vcf_file = open(vcf_filename, 'r')

        for line in vcf_file:
            # File Format
            # Column Number:     | 0       | 1        | 2          | 3       | 4
            # Column Description:| chrom   | pos      | dbSNP_id   | ref_nuc | alt_nuc
            # Column Example:    | 13      | 32906558 | rs79483201 | T       | A
            # Column Number:     | 5       | 6        | 7          | 8              | 9./..
            # Column Description:| qual    | filter   | info       | format         | SampleIDs
            # Column Example:    | 5645.6  | PASS     | .          | GT:AD:DP:GQ:PL | Scores corresponding to format

            """
            Examples of the formats
            ##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths">
            ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth">
            ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
            ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
            ##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Ph
            """

            if line.startswith('#'):  # skip information headers
                continue

            fields = line.strip('\n').split('\t')

            chrom = fields[0]

            # filters for chrom
            if ptn_auto.match(chrom):
                if int(chrom) > 23:  # Wrong chromosome number
                    invalid_chrom = 'chr%s' % chrom

                    if invalid_chrom not in invalid_chr_to_cnt:
                        invalid_chr_to_cnt[invalid_chrom] = 0

                    invalid_chr_to_cnt[invalid_chrom] += 1

                    continue

            elif not ptn_sex.match(chrom):
                invalid_chrom = 'chr%s' % chrom

                if invalid_chrom not in invalid_chr_to_cnt:
                    invalid_chr_to_cnt[invalid_chrom] = 0

                invalid_chr_to_cnt[invalid_chrom] += 1

                continue

            variant = VCFData()
            variant.chrom = 'chr%s' % fields[0]

            if not ptn_pos.match(fields[1]):
                eprint('Invalid variant position \'%s\'' % fields[1])
                continue

            variant.pos = int(fields[1])  # 1-based
            variant.dbSNP_id = fields[2]
            variant.ref_nuc = fields[3]
            variant.alt_nuc = fields[4]
            variant.qual = float(fields[5]) if fields[5] != '.' else fields[5]
            variant.filter = fields[6]
            variant.info = fields[7]
            variant.format = fields[8]
            variant.samples = fields[9:]

            var_list.append(variant)
        # END: for loop 'line'

        vcf_file.close()

        # print the invalid chromosome ID
        eprint('\nInvalid chromosome ID of the variants in %s' % vcf_filename)

        for invalid_chrom in invalid_chr_to_cnt:
            eprint('%s: %d' % (invalid_chrom, invalid_chr_to_cnt[invalid_chrom]))

        eprint()

        return var_list

    # END: the function 'parse_vcf_file'

    def set_genic_region(self, genes_same_chr):
        """
        :param genes_same_chr: the list of 'RefFlat' object in the same chromosome with this variant
        """
        var_pos = self.pos - 1  # 1-base -> 0-base
        genes_same_chr.sort(key=lambda chr_gene: (chr_gene.tx_start, chr_gene.tx_end))

        for gene in genes_same_chr:
            if gene.chrom != self.chrom:
                eprint('Different chromosome ID')
                eprint('ChrID of %s (%s): %s' % gene.symbol, gene.id, gene.chrom)
                eprint('ChrID of the variant: %s' % self.chrom)
                sys.exit()

            if var_pos < gene.tx_start:
                break

            elif gene.tx_start <= var_pos < gene.tx_end:
                strand = gene.strand
                gene_sym = gene.symbol
                gene_id = gene.id
                genic_region = gene.find_genic_region(var_pos)

                if gene_sym not in self._genic_region_dict[strand]:
                    self._genic_region_dict[strand][gene_sym] = {}

                if gene_id in self._genic_region_dict[strand][gene_sym]:
                    eprint('There are RefFlat objects with same id.')
                    sys.exit()

                self._genic_region_dict[strand][gene_sym][gene_id] = genic_region

        # END: for loop 'gene'

        self._set_genic_region_value()

    # END: the function 'find_genic_region'

    def _set_genic_region_value(self):
        """
        makes up _strand_to_region_val
        """
        for strand in ['+', '-']:
            genic_region_to_bool = parse_genic_region_val(0)  # default (intergenic)

            for gene_sym in self._genic_region_dict[strand]:
                for gene_id in self._genic_region_dict[strand][gene_sym]:
                    genic_region = str(self._genic_region_dict[strand][gene_sym][gene_id])
                    genic_region_to_bool[genic_region] = True
                    genic_region_to_bool['intergenic'] = False

            region_val = get_genic_region_val(genic_region_to_bool)
            self._strand_to_region_val[strand] = region_val

    # END: the function '_set_genic_region_value'
# END: class 'VCFData'
