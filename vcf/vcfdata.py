import sys, os

class VCFData:
    """ from annovar, only consider SNV """
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

        """ info from annovar """
        self.anno_gene_sym = ''
        self.anno_genic = ''
        self.anno_alt_genic = ''
        self.anno_SNV_type = ''

        """ info from refFlat (in-house rep-isoforms) """
        self.genic_region_dict = {}  # Mapping route: gene symbol -> ID -> genic region
        self.rep_gene_sym = None
        self.rep_gene_id = None
        self.rep_genic = 'intergenic'

    # END: __init__

    def set_genic_region(self, genes_same_chr):
        """
        :param genes_same_chr: the list of 'RefFlat' objects
                               this genes must be involved in the chromosome same with the variant.
        """
        var_pos = self.pos - 1  # 1-base -> 0-base
        genes_same_chr.sort(key=lambda gene: gene.tx_start)

        for gene in genes_same_chr:
            if gene.chrom != self.chrID:
                print('Different chromosome ID')
                print('ChrID of %s (%s): %s' % gene.symbol, gene.id, gene.chrom)
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
        genic_precedence_dict = {'exonic': 1,
                                 'ncRNA_intronic': 2, 'ncRNA_exonic': 2,
                                 '3UTR': 3, '5UTR': 3,
                                 'intronic': 4,
                                 'intergenic': 5}

        for gene_sym in self.genic_region_dict:
            for gene_id in self.genic_region_dict[gene_sym]:
                genic_region = self.genic_region_dict[gene_sym][gene_id]

                if genic_precedence_dict[self.rep_genic] > genic_precedence_dict[genic_region]:
                    self.rep_gene_sym = gene_sym
                    self.rep_gene_id = gene_id
                    self.rep_genic = genic_region
            # END: for loop 'gene_id'
        # END: for loop 'gene_sym'
    # END: _set_rep_genic_region

    @staticmethod
    def parse_vcf_file(vcf_filename, isClustered=False):
        alt_genic_dict = {'exonic': 'exonic', 'UTR3': '3UTR', 'UTR5': '5UTR', 'intronic': 'intronic',
                          'splicing': 'intronic', 'ncRNA_splicing': 'ncRNA_intronic', 'downstream': 'intergenic',
                          'ncRNA_exonic': 'ncRNA_exonic', 'intergenic': 'intergenic', 'ncRNA_UTR5': 'ncRNA_exonic',
                          'ncRNA_UTR3': 'ncRNA_exonic', 'upstream': 'intergenic', 'ncRNA_intronic': 'ncRNA_intronic',
                          '5UTR': '5UTR', '3UTR':'3UTR', 'other':'other'}

        if not os.path.isfile(vcf_filename):
            sys.exit('File Not Found %s' % vcf_filename)

        var_list = []
        vcf_file = open(vcf_filename, 'r')

        for line in vcf_file:
            # File Format
            # Column Number:     | 0       | 1        | 2          | 3       | 4
            # Column Description:| chrID   | pos      | dbSNP_id   | ref_nuc | alt_nuc
            # Column Example:    | chr13   | 32906558 | rs79483201 | T       | A
            # Column Number:     | 5       | 6        | 7          | 8              | 9./..
            # Column Description:| qual    | filter  | info        | format         | SampleIDs
            # Column Example:    | 5645.6  | PASS     | .          | GT:AD:DP:GQ:PL | Scores corresponding to format
            ##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
            ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
            ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
            ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
            ##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Ph

            if line.startswith('#'): continue  # SKIP information Headers
            fields = line.strip('\n').split('\t')

            ## filter for chrID
            if fields[0] == 'MT':
                continue
            if fields[0].startswith('<GL00'):
                continue
            if fields[0].startswith('.'):
                continue

            ## check if the variant passed the variant calling filters
            if 'PASS' not in line:
                continue

            variant = VCFData()
            variant.chrID  = 'chr%s' % fields[0]

            try: variant.pos = int(fields[1])
            except ValueError: continue

            variant.dbSNP_id = fields[2]
            variant.ref_nuc = fields[3]
            variant.alt_nuc = fields[4]
            variant.qual = float(fields[5]) if fields[5] != '.' else fields[5]
            variant.filter = fields[6]
            variant.info = fields[7]
            variant.format = fields[8]
            variant.list_sSamples = fields[9:]

            vcf_info_dict = dict(info.split('=') for info in variant.info.split(';')[:-1] if len(info.split('=')) == 2)

            variant.anno_gene_sym = vcf_info_dict['Gene.refGene'].split('\\x3b')[0]
            variant.anno_genic = vcf_info_dict['Func.refGene'].split('\\x3b')[0]
            variant.anno_alt_genic = alt_genic_dict[vcf_info_dict['Func.refGene'].split('\\x3b')[0]]
            variant.anno_SNV_type = vcf_info_dict['ExonicFunc.refGene'].split('\\x3b')[0]

            var_list.append(variant)
        # END: for loop 'line'

        vcf_file.close()

        return var_list
    # END: parse_vcf_file
# END: class 'VCFData'
