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
        self.qual = '.'
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
        return '%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s' % (self.chrom, self.pos + 1, self.dbSNP_id, self.ref_nuc,
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
        :return: a dictionary 
                : mapping route: gene symbol -> ID -> (a boolean value, mutation type))
                : If the boolean is True, it means that the gene is protein-modified gene
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
                if genic_region == 'CDS':
                    prot_is_mod, mut_type = gene.point_mut_type(self.pos, self.ref_nuc, self.alt_nuc)
                    self._mut_type_dict[strand][gene_sym][gene_id] = (prot_is_mod, mut_type)

                    if prot_is_mod:  # a protein product of the gene is modified
                        self._is_non_synonymous = True

        self._set_anno_val()

    def parse_vcf_entry(self, vcf_entry):
        """
        By parsing a line of VCF file, get and store information to this VCF object.
        :param vcf_entry: a line of a VCF file (type: str)
        """
        """
        * VCF file format (tab seperated)
        1. CHROM
        2. POS (1-based)
        3. ID
        4. REF
        5. ALT
        6. QUAL
        7. FILTER
        8. INFO (Additional information, format depends of a variant caller)
        """
        fields = vcf_entry.strip('\n').split('\t')

        self.chrom = 'chr%s' % fields[0]
        self.pos = int(fields[1]) - 1  # 1-based -> 0-based
        self.dbSNP_id = fields[2]
        self.ref_nuc = fields[3]
        self.alt_nuc = fields[4]
        self.qual = fields[5]
        self.filter = fields[6]
        self.info = '\t'.join(fields[7:])

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
    def parse_vcf_file(vcf_path):
        """
        Parse a VCF file and return SNV objects
        This function only returns objects with a canonical chromosome ID
        :param vcf_path: a path of a VCF file
        :return: a list of SNV objects
        """
        if not os.path.isfile(vcf_path):
            raise IOError('%s does not exist' % vcf_path)

        # regex for VCF entries
        regex_chr = re.compile('^chr([0-9]{1,2}|[XY])$')
        variants = []

        if vcf_path.endswith('.gz'):
            vcf_file = gzip.open(vcf_path, 'rt')
        else:
            vcf_file = open(vcf_path, 'r')

        for line in vcf_file:
            if line.startswith('#'):  # skip information headers
                continue

            variant = SNV()
            variant.parse_vcf_entry(line)

            # filters for chrom
            if regex_chr.match(variant.chrom):
                variants.append(variant)

        vcf_file.close()

        return variants

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
