import sys
import os
import re
import gzip

from gene.utils import *


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
            sys.exit('File Not Found %s' % vcf_filename)

        # regex for VCF entries
        ptn_autosome = re.compile('^[0-9]{1,2}$')
        ptn_sexchr = re.compile('^[XY]$')
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

            # check if the variant passed the variant calling filters
            if 'PASS' not in line:
                continue

            fields = line.strip('\n').split('\t')

            chrom = fields[0]

            # filters for chrom
            if ptn_autosome.match(chrom):
                if int(chrom) > 23:  # Wrong chromosome number
                    invalid_chrom = 'chr%s' % chrom

                    if invalid_chrom not in invalid_chr_to_cnt:
                        invalid_chr_to_cnt[invalid_chrom] = 0

                    invalid_chr_to_cnt[invalid_chrom] += 1

                    continue

            elif not ptn_sexchr.match(chrom):
                invalid_chrom = 'chr%s' % chrom

                if invalid_chrom not in invalid_chr_to_cnt:
                    invalid_chr_to_cnt[invalid_chrom] = 0

                invalid_chr_to_cnt[invalid_chrom] += 1

                continue

            variant = VCFData()
            variant.chrom = 'chr%s' % fields[0]

            if not ptn_pos.match(fields[1]):
                print('Invalid variant position \'%s\'' % fields[1])
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
        print('\nInvalid chromosome ID of the variants in %s' % vcf_filename)

        for invalid_chrom in invalid_chr_to_cnt:
            print('%s: %d' % (invalid_chrom, invalid_chr_to_cnt[invalid_chrom]))

        print()

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
                print('Different chromosome ID')
                print('ChrID of %s (%s): %s' % gene.symbol, gene.id, gene.chrom)
                print('ChrID of the variant: %s' % self.chrom)
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
                    print('There are RefFlat objects with same id.')
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
            genic_region_to_bool = {genic_region: False for genic_region in GENIC_REGION_LIST}
            genic_region_to_bool['intergenic'] = True  # default

            for gene_sym in self._genic_region_dict[strand]:
                for gene_id in self._genic_region_dict[strand][gene_sym]:
                    genic_region = str(self._genic_region_dict[strand][gene_sym][gene_id])
                    genic_region_to_bool[genic_region] = True
                    genic_region_to_bool['intergenic'] = False

            region_val = get_genic_region_val(genic_region_to_bool)
            self._strand_to_region_val[strand] = region_val

    # END: the function '_set_genic_region_value'
# END: class 'VCFData'
