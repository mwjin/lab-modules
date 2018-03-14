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
        self.rep_gene_id = '.'  # the gene from which the representative genic region comes
        self.rep_strand = '.'  # the strand of the representative gene
        self.rep_genic_region = 'intergenic'  # default

        self._genic_region_dict = {'+': {}, '-': {}}  # Mapping route: strand -> gene symbol -> ID -> genic region
        self._strand_to_gene_ids = {'+': [], '-': []}
        self._strand_to_genic_region = {'+': 'intergenic', '-': 'intergenic'}
        self._strand_to_region_val = {'+': 0, '-': 0}

    # END: __init__

    def get_genic_region(self):
        """
        :return: representative genic region
        """
        return self.rep_genic_region

    # END: the function 'get_genic_region'

    def get_strand_genic_regions(self):
        """
        :return: a tuple with strand-specific genic regions (genic regions in +, genic regions in -)
        """
        return self._strand_to_genic_region['+'], self._strand_to_genic_region['-']

    # END: the function 'get_strand_genic_regions'

    def get_assigned_strand(self):
        """
        :return: a strand where this variant was more likely to be assigned.
        * Notice: there is no strand for variants actually. However, to determine the transcript (+ or -)
                  this variant was more likely to be assigned, we set this virtual strand of this variant.
        """
        return self.rep_strand

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

        self._set_rep_genic_region()

    # END: the function 'find_genic_region'

    def _set_rep_genic_region(self):
        genic_priority_dict = {'ORF': 1,
                               '5UTR': 2, '3UTR': 2, 'UTR': 2,  # UTR: both 5 and 3UTR
                               'ncRNA_exonic': 3,
                               'intronic': 4,
                               'ncRNA_intronic': 5,
                               'intergenic': 6}

        # set strand-specific representative genic region
        for strand in ['+', '-']:
            for gene_sym in self._genic_region_dict[strand]:
                for gene_id in self._genic_region_dict[strand][gene_sym]:
                    genic_region = str(self._genic_region_dict[strand][gene_sym][gene_id])

                    if genic_priority_dict[self._strand_to_genic_region[strand]] > genic_priority_dict[genic_region]:
                        self._strand_to_genic_region[strand] = genic_region
                        self._strand_to_gene_ids[strand] = [gene_id]

                    elif genic_priority_dict[self._strand_to_genic_region[strand]] == genic_priority_dict[genic_region]:
                        self._strand_to_gene_ids[strand].append(gene_id)
                        strand_genic_region = self._strand_to_genic_region[strand]

                        # deal with the case where one variable is annotated with both 5 and 3'UTR
                        if genic_priority_dict[strand_genic_region] == 2 and strand_genic_region != genic_region:
                            self._strand_to_genic_region[strand] = 'UTR'

                # END: for loop 'gene_id'
            # END: for loop 'gene_sym'
        # END: for loop 'strand'

        self._strand_to_gene_ids['+'].sort(key=lambda one_gene_id: int(one_gene_id[3:]))
        self._strand_to_gene_ids['-'].sort(key=lambda one_gene_id: int(one_gene_id[3:]))

        # determine the representative strand of the variant
        pos_genic_region = self._strand_to_genic_region['+']
        neg_genic_region = self._strand_to_genic_region['-']

        if genic_priority_dict[pos_genic_region] < genic_priority_dict[neg_genic_region]:
            self.rep_strand = '+'
        elif genic_priority_dict[pos_genic_region] > genic_priority_dict[neg_genic_region]:
            self.rep_strand = '-'
        else:
            pos_gene_cnt = len(self._strand_to_gene_ids['+'])
            neg_gene_cnt = len(self._strand_to_gene_ids['-'])

            if pos_gene_cnt > neg_gene_cnt:
                self.rep_strand = '+'
            elif pos_gene_cnt < neg_gene_cnt:
                self.rep_strand = '-'
            elif pos_gene_cnt != 0:
                pos_lowest_gene_id = self._strand_to_gene_ids['+'][0]
                neg_lowest_gene_id = self._strand_to_gene_ids['-'][0]

                if int(pos_lowest_gene_id[3:]) >= int(neg_lowest_gene_id[3:]):
                    self.rep_strand = '+'
                else:
                    self.rep_strand = '-'
            else:
                self.rep_strand = '+'

        # determine the id for representative genes where the variant is located
        rep_gene_id_list = self._strand_to_gene_ids[self.rep_strand]

        if not rep_gene_id_list:
            self.rep_gene_id = '.'
        else:
            self.rep_gene_id = ','.join(rep_gene_id_list)

        self.rep_genic_region = self._strand_to_genic_region[self.rep_strand]

    # END: the function '_set_rep_genic_region'

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
