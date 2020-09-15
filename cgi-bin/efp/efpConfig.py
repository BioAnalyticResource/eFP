#!/usr/bin/python3
"""
Created on Nov 24, 2009
@author: Robert Breit
"""
import os

# calculate top level directory for eFP browser code
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BASE_DIR = os.path.join(BASE_DIR, '../')

# Read-only database config
DB_HOST = ''  # hostname of MySQL database server
DB_USER = ''  # username for database access
DB_PASSWD = ''  # password for database access

DB_DATA_TABLE = 'sample_data'  # database table which contains expression data

# Species annotations
DB_ANNO = None  # database name for annotations lookup

# lookup table for gene annotation
DB_ANNO_TABLE = None  # SPECIES_annotation
DB_ANNO_GENEID_COL = None  # geneid

# lookup table for probeset id
DB_LOOKUP_TABLE = None  # SPECIES_lookup
DB_LOOKUP_GENEID_COL = None  # geneid

# Directory containing static files
STATIC_FILES = os.path.join(BASE_DIR, 'static')
STATIC_FILES_WEB = '../static'

# Directory to place output files
OUTPUT_FILES = os.path.join(BASE_DIR, 'output')
OUTPUT_FILES_WEB = '../output'

# lookup tables for ncbi ids
DB_NCBI_GENE_TABLE = None
DB_NCBI_PROT_TABLE = None
DB_NCBI_GENEID_COL = None

# Check Y-AXIS message
Y_AXIS = {'database_name': "Signal"}

# Check if lookup exists
LOOKUP = {'little_millet': "0"}

# initial gene ids when start page is loaded.
GENE_ID_DEFAULT1 = ''
GENE_ID_DEFAULT2 = ''

# the little graph on the tga image has a scale
# such that 1 unit is x pixels for different ranges on the x-axis of the graph
# the GRAPH_SCAL_UNIT consists of value pairs: upper end of range and scale unit
# so ((1000, 0.031), (10000, 0.003222), (1000000, 0.00031)) means:
# use 0.031 as scale unit for 0 < signal < 1000
# use 0.003222 as scale unit for 1000 < signal < 10000
# use 0.00031 as scale unit for 10000 < signal < 1000000
# see also efp.draw_image()
GRAPH_SCALE_UNIT = {'database_name': [(10, 3), (100, 0.3), (1000, 0.03)]}
# the default values are used if for the given data source no special values are defined
# GRAPH_SCALE_UNIT["default"] = ((0.06, 16) (0.006, 73), (0.0006, 130))

# define additional header text for individual data sources
# you can use key 'default' for each not individually defined
datasourceHeader = {'default': ''}

# define additional footer text for individual data sources
# you can use key 'default' for each not individually defined
datasourceFooter = {'default': ''}

# regular expression for check of gene id input
inputRegEx = r"^()$"

# default thresholds
minThreshold_Compare = 0.6  # Minimum color threshold for comparison is 0.6, giving us [-0.6, 0.6] on the color legend
minThreshold_Relative = 0.6  # Minimum color threshold for median is 0.6, giving us [-0.6, 0.6] on the color legend ~ 1.5-Fold
minThreshold_Absolute = 10  # Minimum color threshold for efp_max is 10, giving us [0, 10] on the color legend

# coordinates where to write gene id, probeset id and gene alias into image
GENE_ID1_POS = (5, 5)
GENE_ID2_POS = (5, 20)
GENE_PROBESET1_POS = (190, 5)
GENE_PROBESET2_POS = (190, 20)
GENE_ALIAS1_POS = (0, 0)
GENE_ALIAS2_POS = (0, 0)

GENE_ORTHO1_POS = (0, 285)
defaultDataSource = 'View_Name'
dataDir = os.path.join(BASE_DIR, 'data')
dataDirWeb = os.path.join('../data')

dbGroupDefault = 'group1'
# list of datasources in same group to find efp_max signal
groupDatasource = {'group1': ['View_Name']}

# mapping of xml files to show datasource name
groupDatasourceName = {"group1": {'View_Name': 'View Name'}}

species = 'SPECIES'
spec_names = {'SPECIES': 'Species'}
