""" Definitions for my library """

from os.path import dirname, abspath

BASE_DIR = dirname(dirname(abspath(__file__)))  # the directory that contains all my projects
DATA_DIR = "%s/data" % BASE_DIR  # the directory that contains common data for my projects
ROOT_DIR = dirname(BASE_DIR)