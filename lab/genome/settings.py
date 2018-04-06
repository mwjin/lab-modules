from lab.settings import ROOT_DIR

__all__ = ['GENOME_VER', 'GENOME_FILE_PATH']

# settings for genome
GENOME_VER = "hg19"
GENOME_FILE_PATH = "{0}/ref_genome/{1}/{1}.fa".format(ROOT_DIR, GENOME_VER)
