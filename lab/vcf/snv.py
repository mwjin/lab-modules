import sys
import os
import re
import gzip

from lab.utils import eprint, caller_file_and_line
from lab.gene.anno import get_anno_val

__all__ = ['SNV']


class SNV:
    """
    The object of this class represents one SNV entry in VCF files
    (SNV: Single nucleotide variant)
    """
    def __init__(self):
        """ Essential fields of a VCF file"""
        self.chrom = '.'
        self.pos = 0  # 0-base
        self.dbSNP_id = '.'
        self.ref_nuc = '.'
        self.alt_nuc = '.'
        self.qual = 0.0
        self.filter = '.'
        self.info = '.'

        """ INFO not in a VCF file """
        self._is_non_synonymous = False  # If this variant is non-synonymous for at least one gene, it is true.

        # Dictionary for the information of genes associated with this variant
        # Mapping route: strand -> gene symbol -> ID -> genic region
        self._gene_dict = {'+': {}, '-': {}}

        # For gene-based annotation of this variant
        # Annotation value: see lab.gene.anno.py
        self._strand_to_anno_val = {'+': 0, '-': 0}

        # Dictionary of mutation types of variant for all genes associated with this variant
        # Mapping route: strand -> gene symbol -> ID -> (is_non_synonymous, mutation type)
        self._mut_type_dict = {'+': {}, '-': {}}

    def __repr__(self):
        """
        Essential information for the variants which is represented as an dictionary
        This can be used to make a new 'SNV' object through the built-in 'eval' function
        """
        return "{chrom: '%s', pos: %d, ref_nuc: '%s', alt_nuc: '%s'}" % \
               (self.chrom, self.pos, self.ref_nuc, self.alt_nuc)

    def __str__(self):
        """
        Return an entry for a VCF file. It only contains essential fields.
        """
        return '%s\t%d\t%s\t%s\t%s\t%.1f\t%s\t%s' % (self.chrom, self.pos + 1, self.dbSNP_id, self.ref_nuc,
                                                     self.alt_nuc, self.qual, self.filter, self.info)

    def __eq__(self, other):
        """
        Compare the essential characteristics between two 'SNV' objects
        """
        self_info = (self.chrom, self.pos, self.ref_nuc, self.alt_nuc)
        other_info = (other.chrom, other.pos, other.ref_nuc, other.alt_nuc)

        return self_info == other_info

    def get_assigned_strand(self):
        """
        * Notice: there is no strand for variants actually. However, to determine the transcript (+ or -)
                  this variant was more likely to be assigned, we set this virtual strand of this variant.

        :return: a strand where this variant was more likely to be assigned.
        """

        if self._strand_to_anno_val['+'] >= self._strand_to_anno_val['-']:
            return '+'
        else:
            return '-'

    def get_assoc_genes(self, strand):
        """
        Return the dictionary that contains the information of the strand-specific genes associated with this variants
        :param strand: '+' or '-'
        :return: a dictionary (mapping route: gene symbol -> ID -> genic region)
        """
        return self._gene_dict[strand]

    def get_anno_val(self, strand):
        """
        :param strand: '+' or '-'
        :return: the genic region value on the input strand
        """
        return self._strand_to_anno_val[strand]

    def get_mut_type_dict(self, strand):
        """
        Return a dictionary that contains information of mutation types for each gene associated with this variant
        :param strand: '+' or '-'
        :return: a dictionary (mapping route: gene symbol -> ID -> (is_non_synonymous, mutation type))
        """
        return self._mut_type_dict[strand]

    def is_non_synonymous(self):
        """
        Return a boolean value that represents whether this SNV is non-synonymous or not
        """
        return self._is_non_synonymous

    def gene_based_anno(self, genes_same_chr):
        """
        Make up the '_gene_dict' attribute and set up the annotation values of this variant
        For the gene-based annotation, this method uses a list of genes.
        :param genes_same_chr: a list of 'RefFlat' objects in the same chromosome with this variant
        """
        genes_same_chr.sort(key=lambda chr_gene: (chr_gene.tx_start, chr_gene.tx_end))

        for gene in genes_same_chr:
            if gene.chrom != self.chrom:
                eprint('[ERROR] in %s' % caller_file_and_line())
                eprint('\tThere is a gene with the different chromosome ID')
                eprint('\tChromosome ID of %s (%s): %s' % gene.symbol, gene.id, gene.chrom)
                eprint('\tChromosome ID of the variant: %s' % self.chrom)
                sys.exit()

            if self.pos < gene.tx_start - 300:  # there will be no overlap in the next genes
                break

            elif self.pos < gene.tx_start:  # promoter
                strand = gene.strand
                gene_sym = gene.symbol
                gene_id = gene.id

                if gene_sym not in self._gene_dict[strand]:
                    self._gene_dict[strand][gene_sym] = {}

                if gene_sym not in self._mut_type_dict[strand]:
                    self._mut_type_dict[strand][gene_sym] = {}

                self._gene_dict[strand][gene_sym][gene_id] = 'promoter'
                self._mut_type_dict[strand][gene_sym][gene_id] = (False, None)

            elif gene.tx_start <= self.pos < gene.tx_end:
                strand = gene.strand
                gene_sym = gene.symbol
                gene_id = gene.id
                genic_region = gene.find_genic_region(self.pos)

                if gene_sym not in self._gene_dict[strand]:
                    self._gene_dict[strand][gene_sym] = {}

                if gene_sym not in self._mut_type_dict[strand]:
                    self._mut_type_dict[strand][gene_sym] = {}

                self._gene_dict[strand][gene_sym][gene_id] = genic_region

                # check the mutation type of this variant
                if genic_region == 'ORF':
                    is_non_synonymous, mut_type = gene.is_non_synonymous(self.pos, self.ref_nuc, self.alt_nuc)
                    self._mut_type_dict[strand][gene_sym][gene_id] = (is_non_synonymous, mut_type)

                    if is_non_synonymous:
                        self._is_non_synonymous = True

        self._set_anno_val()

    @staticmethod
    def parse_repr_str(repr_str):
        """
        Parse the representative string of the 'SNV' object and return a 'SNV' object
        There is less of information.
        """
        variant = SNV()
        attr_dict = eval(repr_str)

        for snv_attr in attr_dict:
            setattr(variant, snv_attr, attr_dict[snv_attr])

        return variant

    @staticmethod
    def parse_vcf_file(vcf_filename):

        if not os.path.isfile(vcf_filename):
            raise IOError('%s does not exist' % vcf_filename)

        # regex for VCF entries
        regex_auto = re.compile('^[0-9]{1,2}$')
        regex_sex = re.compile('^[XY]$')
        regex_pos = re.compile('^[0-9]+$')

        var_list = []
        invalid_chr_to_cnt = {}

        if vcf_filename.endswith('.gz'):
            vcf_file = gzip.open(vcf_filename, 'rt')
        else:
            vcf_file = open(vcf_filename, 'r')

        for line in vcf_file:
            """
            <File Format>
             Column Number:     | 0       | 1        | 2          | 3       | 4
             Column Description:| chrom   | pos      | dbSNP_id   | ref_nuc | alt_nuc
             Column Example:    | 13      | 32906558 | rs79483201 | T       | A
             Column Number:     | 5       | 6        | 7          | 8              | 9./..
             Column Description:| qual    | filter   | info       | format         | SampleIDs
             Column Example:    | 5645.6  | PASS     | .          | GT:AD:DP:GQ:PL | Scores corresponding to format

            <Examples of the formats>
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
            if regex_auto.match(chrom):
                if int(chrom) > 23:  # Wrong chromosome number
                    invalid_chrom = 'chr%s' % chrom

                    if invalid_chrom not in invalid_chr_to_cnt:
                        invalid_chr_to_cnt[invalid_chrom] = 0

                    invalid_chr_to_cnt[invalid_chrom] += 1

                    continue

            elif not regex_sex.match(chrom):
                invalid_chrom = 'chr%s' % chrom

                if invalid_chrom not in invalid_chr_to_cnt:
                    invalid_chr_to_cnt[invalid_chrom] = 0

                invalid_chr_to_cnt[invalid_chrom] += 1

                continue

            variant = SNV()
            variant.chrom = 'chr%s' % fields[0]

            if not regex_pos.match(fields[1]):
                eprint('[LOG] Invalid variant position \'%s\'' % fields[1])
                continue

            variant.pos = int(fields[1]) - 1  # 1-based -> 0-based
            variant.dbSNP_id = fields[2]
            variant.ref_nuc = fields[3]
            variant.alt_nuc = fields[4]
            variant.qual = float(fields[5]) if fields[5] != '.' else fields[5]
            variant.filter = fields[6]
            variant.info = fields[7]

            var_list.append(variant)

        vcf_file.close()

        # print the invalid chromosome ID
        if len(invalid_chr_to_cnt) != 0:
            eprint('\n[LOG] Invalid chromosome ID of the variants in %s' % vcf_filename)

            for invalid_chrom in invalid_chr_to_cnt:
                eprint('--- %s: %d' % (invalid_chrom, invalid_chr_to_cnt[invalid_chrom]))

            eprint()

        return var_list

    def _set_anno_val(self):
        """
        makes up _strand_to_anno_val
        """
        for strand in ['+', '-']:
            strand_anno_val = 0  # default (intergenic)

            for gene_sym in self._gene_dict[strand]:
                for gene_id in self._gene_dict[strand][gene_sym]:
                    genic_region = str(self._gene_dict[strand][gene_sym][gene_id])
                    anno_val = get_anno_val(genic_region)
                    region_bit = int(strand_anno_val / anno_val) % 2

                    if region_bit == 0:  # same genic region value is not added yet.
                        strand_anno_val += anno_val

            self._strand_to_anno_val[strand] = strand_anno_val
