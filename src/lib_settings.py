""" Settings for all the projects """

from os.path import dirname, abspath

# paths for all the projects
BASE_DIR = dirname(dirname(dirname(abspath(__file__))))  # the directory that contains all my projects
DATA_DIR = "%s/data" % BASE_DIR  # the directory that contains common data for my projects
ROOT_DIR = dirname(BASE_DIR)
