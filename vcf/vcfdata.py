import sys
import os
import re
import gzip


class VCFData:
    """ The object of this class represents one entry in VCF files (only consider SNV). """
    def __init__(self):
        self.chrID = ''
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
        self.genic_region_dict = {}  # Mapping route: gene symbol -> ID -> (genic region, strand)
        self.rep_genic_region = 'intergenic'  # default
        self.rep_gene_id = '.'  # the gene from which the representative genic region comes
        self.rep_strand = '.'  # the strand of the representative gene

    # END: __init__

    def get_var_genic_region(self):
        """
        :return: representative genic region
        """
        return self.rep_genic_region

    @staticmethod
    def parse_vcf_file(vcf_filename):

        if not os.path.isfile(vcf_filename):
            sys.exit('File Not Found %s' % vcf_filename)

        # regex for VCF entries
        pat_autosome = re.compile('^[0-9]{1,2}$')
        pat_sexchr = re.compile('^[XY]$')
        pat_pos = re.compile('^[0-9]+$')

        var_list = []
        invalid_chr_to_cnt = {}

        if vcf_filename.endswith('.gz'):
            vcf_file = gzip.open(vcf_filename, 'rt')
        else:
            vcf_file = open(vcf_filename, 'r')

        for line in vcf_file:
            # File Format
            # Column Number:     | 0       | 1        | 2          | 3       | 4
            # Column Description:| chrID   | pos      | dbSNP_id   | ref_nuc | alt_nuc
            # Column Example:    | 13      | 32906558 | rs79483201 | T       | A
            # Column Number:     | 5       | 6        | 7          | 8              | 9./..
            # Column Description:| qual    | filter   | info       | format         | SampleIDs
            # Column Example:    | 5645.6  | PASS     | .          | GT:AD:DP:GQ:PL | Scores corresponding to format

            # Examples of the formats
            ##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
            ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
            ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
            ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
            ##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Ph

            if line.startswith('#'):  # skip information headers
                continue

            # check if the variant passed the variant calling filters
            if 'PASS' not in line:
                continue

            fields = line.strip('\n').split('\t')

            chrID = fields[0]

            # filters for chrID
            if pat_autosome.match(chrID):
                if int(chrID) > 23:  # Wrong chromosome number
                    invalid_chrID = 'chr%s' % chrID

                    if invalid_chrID not in invalid_chr_to_cnt:
                        invalid_chr_to_cnt[invalid_chrID] = 0

                    invalid_chr_to_cnt[invalid_chrID] += 1

                    continue

            elif not pat_sexchr.match(chrID):
                invalid_chrID = 'chr%s' % chrID

                if invalid_chrID not in invalid_chr_to_cnt:
                    invalid_chr_to_cnt[invalid_chrID] = 0

                invalid_chr_to_cnt[invalid_chrID] += 1

                continue

            variant = VCFData()
            variant.chrID = 'chr%s' % fields[0]

            if not pat_pos.match(fields[1]):
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

        for invalid_chrID in invalid_chr_to_cnt:
            print('%s: %d' % (invalid_chrID, invalid_chr_to_cnt[invalid_chrID]))

        print()

        return var_list
    # END: parse_vcf_file

    def set_genic_region(self, genes_same_chr):
        """
        :param genes_same_chr: the list of 'RefFlat' objects
                               this genes must be involved in the chromosome same with the variant.
        """
        var_pos = self.pos - 1  # 1-base -> 0-base
        genes_same_chr.sort(key=lambda gene: (gene.tx_start, gene.tx_end))

        for gene in genes_same_chr:
            if gene.chrID != self.chrID:
                print('Different chromosome ID')
                print('ChrID of %s (%s): %s' % gene.symbol, gene.id, gene.chrID)
                print('ChrID of the variant: %s' % self.chrID)
                sys.exit()

            if var_pos < gene.tx_start:
                break

            elif gene.tx_start <= var_pos < gene.tx_end:
                gene_sym = gene.symbol
                gene_id = gene.id
                genic_region = gene.find_genic_region(var_pos)
                strand = gene.strand

                if gene_sym not in self.genic_region_dict:
                    self.genic_region_dict[gene_sym] = {}

                if gene_id in self.genic_region_dict[gene_sym]:
                    print('There are RefFlat objects with same id.')
                    sys.exit()

                self.genic_region_dict[gene_sym][gene_id] = (genic_region, strand)

        # END: for loop 'gene'

        self._set_rep_genic_region()
    # END: find_genic_region

    def _set_rep_genic_region(self):
        genic_priority_dict = {'ORF': 1,
                               '5UTR': 2, '3UTR': 2, 'UTR': 2,  # UTR: both 5 and 3UTR
                               'ncRNA_exonic': 3,
                               'intronic': 4,
                               'ncRNA_intronic': 5,
                               'intergenic': 6}

        rep_gene_info_dict = {'+': [],  # mapping route: strand -> (gene_id, genic_region)
                              '-': []}

        for gene_sym in self.genic_region_dict:
            for gene_id in self.genic_region_dict[gene_sym]:
                genic_region, strand = self.genic_region_dict[gene_sym][gene_id]

                if genic_priority_dict[self.rep_genic_region] > genic_priority_dict[genic_region]:
                    self.rep_genic_region = genic_region
                    self.rep_strand = strand
                    rep_gene_info_dict[strand].append((gene_id, genic_region))

                elif genic_priority_dict[self.rep_genic_region] == genic_priority_dict[genic_region]:
                    rep_gene_info_dict[strand].append((gene_id, genic_region))

            # END: for loop 'gene_id'
        # END: for loop 'gene_sym'

        # deal with the case there are multiple genes with genic regions of the same priority
        pos_gene_cnt = len(rep_gene_info_dict['+'])
        neg_gene_cnt = len(rep_gene_info_dict['-'])

        # determine the representative strand of the variant
        if pos_gene_cnt > neg_gene_cnt:
            self.rep_strand = '+'
        elif pos_gene_cnt < neg_gene_cnt:
            self.rep_strand = '-'
        else:
            for strand in rep_gene_info_dict:
                rep_gene_info_dict[strand].sort(key=lambda gene_info: int(gene_info[0][:3]))

            pos_first_gene_id = rep_gene_info_dict['+'][0]
            neg_first_gene_id = rep_gene_info_dict['-'][0]

            if int(pos_first_gene_id[:3]) < int(neg_first_gene_id[:3]):
                self.rep_strand = '+'
            else:
                self.rep_strand = '-'

        self.rep_gene_id = ','.join([gene_id for gene_id, _ in rep_gene_info_dict[self.rep_strand]])

        if genic_priority_dict[self.rep_genic_region] == 2:
            for _, genic_region in rep_gene_info_dict[self.rep_strand]:
                if self.rep_genic_region != genic_region:
                    self.rep_genic_region = 'UTR'
                    break

    # END: _set_rep_genic_region
# END: class 'VCFData'
