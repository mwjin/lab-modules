""" Definitions for my library """

from os.path import dirname, abspath

BASE_DIR = dirname(dirname(abspath(__file__)))
DATA_DIR = "%s/data" % BASE_DIR