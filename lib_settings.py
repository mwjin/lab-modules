""" Settings for all the projects """

from os.path import dirname, abspath

# paths for all the projects
BASE_DIR = dirname(dirname(abspath(__file__)))  # the directory that contains all my projects
DATA_DIR = "%s/data" % BASE_DIR  # the directory that contains common data for my projects
ROOT_DIR = dirname(BASE_DIR)

# settings for genome
GENOME_VER = "hg19"
GENOME_FILENAME = "{0}/ref_genome/{1}/{1}.fa".format(ROOT_DIR, GENOME_VER)
