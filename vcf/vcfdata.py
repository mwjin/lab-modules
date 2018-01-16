import sys, os, re
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
        self.genic_region_dict = {}  # Mapping route: gene symbol -> ID -> genic region
        self.rep_gene_sym = None
        self.rep_gene_id = None
        self.rep_genic_region = 'intergenic'  # default

    # END: __init__

    def get_var_genic_region(self):
        """
        :return: representative genic region
        """
        return self.rep_genic_region

    @staticmethod
    def parse_vcf_file(vcf_filename, isClustered=False):

        if not os.path.isfile(vcf_filename):
            sys.exit('File Not Found %s' % vcf_filename)

        # regex for VCF entries
        p_autosome = re.compile('^[0-9]{1,2}$')
        p_sexchr = re.compile('^[XY]$')
        p_pos = re.compile('^[0-9]+$')

        var_list = []

        if vcf_filename.endswith('.gz'):
            vcf_file = gzip.open(vcf_filename, 'r')
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

            ## check if the variant passed the variant calling filters
            if 'PASS' not in line:
                continue

            fields = line.strip('\n').split('\t')

            chrID = fields[0]

            ## filters for chrID
            if p_autosome.match(chrID):
                if int(chrID) > 23:  # Wrong chromosome number
                    print('Invalid chromosome ID in \'%s\'' % line)
                    continue

            elif not p_sexchr.match(chrID):
                print('Invalid chromosome ID in \'%s\'' % line)
                continue

            variant = VCFData()
            variant.chrID  = 'chr%s' % fields[0]

            if not p_pos.match(fields[1]):
                print('Invalid variant position in \'%s\'' % line)
                continue

            variant.pos = int(fields[1])
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

        return var_list
    # END: parse_vcf_file

    def set_genic_region(self, genes_same_chr):
        """
        :param genes_same_chr: the list of 'RefFlat' objects
                               this genes must be involved in the chromosome same with the variant.
        """
        var_pos = self.pos - 1  # 1-base -> 0-base
        genes_same_chr.sort(key=lambda gene: gene.tx_start)

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

                if gene_sym not in self.genic_region_dict:
                    self.genic_region_dict[gene_sym] = {}

                if gene_id in self.genic_region_dict[gene_sym]:
                    print('There are RefFlat objects with same id.')
                    sys.exit()

                self.genic_region_dict[gene_sym][gene_id] = genic_region
            # END: else
        # END: for loop 'gene'

        self._set_rep_genic_region()
    # END: find_genic_region

    def _set_rep_genic_region(self):
        genic_precedence_dict = {'ORF': 1,
                                 '3UTR': 2, '5UTR': 2,
                                 'ncRNA_exonic': 3,
                                 'intronic': 4,
                                 'ncRNA_intronic': 5,
                                 'intergenic': 6}

        for gene_sym in self.genic_region_dict:
            for gene_id in self.genic_region_dict[gene_sym]:
                genic_region = self.genic_region_dict[gene_sym][gene_id]

                if genic_precedence_dict[self.rep_genic_region] > genic_precedence_dict[genic_region]:
                    self.rep_gene_sym = gene_sym
                    self.rep_gene_id = gene_id
                    self.rep_genic_region = genic_region

                # set priority to the gene with lower ID
                elif genic_precedence_dict[self.rep_genic_region] == genic_precedence_dict[genic_region]:
                    if int(self.rep_gene_id[3:]) > int(gene_id[3:]):
                        self.rep_gene_sym = gene_sym
                        self.rep_gene_id = gene_id
                        self.rep_genic_region = genic_region

            # END: for loop 'gene_id'
        # END: for loop 'gene_sym'
    # END: _set_rep_genic_region
# END: class 'VCFData'
