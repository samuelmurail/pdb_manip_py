#!/usr/bin/env python3
# coding: utf-8

"""
#####################################
#########    PDB IN/OUT    ##########
#####################################
"""

# standard library
import os
import sys
import time
import urllib.request
import logging

# 3rd party packages
import numpy as np
from numpy.linalg import norm
from numpy import sin, cos
from scipy.spatial import distance_matrix
from scipy.spatial.transform import Rotation
from os_command_py import os_command


# Autorship information
__author__ = "Samuel Murail, Damien Espana, Pierre Tuffery"
__copyright__ = "Copyright 2020, RPBS"
__credits__ = ["Samuel Murail"]
__license__ = "GNU General Public License v2.0"
__version__ = "1.0.1"
__maintainer__ = "Samuel Murail"
__email__ = "samuel.murail@u-paris.fr"
__status__ = "Prototype"

# Logging
logger = logging.getLogger(__name__)


def show_log():
    """ To use only with Doctest !!!
    Redirect logger output to sys.stdout
    """
    # Delete all handlers
    logger.handlers = []
    # Set the logger level to INFO
    logger.setLevel(logging.INFO)
    # Add sys.sdout as handler
    logger.addHandler(logging.StreamHandler(sys.stdout))


# Test folder path
PDB_LIB_DIR = os.path.dirname(os.path.abspath(__file__))
TEST_PATH = os.path.abspath(os.path.join(PDB_LIB_DIR, "test/input/"))


# Global variables:

ATOM_MASS_DIST = {'H': 1,
                  'C': 6,
                  'N': 7,
                  'O': 8,
                  'P': 15,
                  'S': 16}


AA_DICT = {'GLY': 'G',
           'HIS': 'H',
           'HSE': 'H',
           'HSD': 'H',
           'HSP': 'H',
           'ARG': 'R',
           'LYS': 'K',
           'ASP': 'D',
           'GLU': 'E',
           'SER': 'S',
           'THR': 'T',
           'ASN': 'N',
           'GLN': 'Q',
           'CYS': 'C',
           'SEC': 'U',
           'PRO': 'P',
           'ALA': 'A',
           'ILE': 'I',
           'PHE': 'F',
           'TYR': 'Y',
           'TRP': 'W',
           'VAL': 'V',
           'LEU': 'L',
           'MET': 'M'}

PROTEIN_AA = AA_DICT.keys()

AA_1_TO_3_DICT = {'G': 'GLY',
                  'H': 'HIS',
                  'R': 'ARG',
                  'K': 'LYS',
                  'D': 'ASP',
                  'E': 'GLU',
                  'S': 'SER',
                  'T': 'THR',
                  'N': 'ASN',
                  'Q': 'GLN',
                  'C': 'CYS',
                  'U': 'SEC',
                  'P': 'PRO',
                  'A': 'ALA',
                  'I': 'ILE',
                  'F': 'PHE',
                  'Y': 'TYR',
                  'W': 'TRP',
                  'V': 'VAL',
                  'L': 'LEU',
                  'M': 'MET',
                  'X': 'ACE'}

# Atom names for each residues
BACK_ATOM = ['N', 'CA', 'C', 'O']

AA_ATOM_DICT = {'X': ['CH3', 'O', 'C'],  # X:ACE
                'G': BACK_ATOM,
                'A': BACK_ATOM + ['CB'],
                'S': BACK_ATOM + ['CB', 'OG'],
                'C': BACK_ATOM + ['CB', 'SG'],
                'T': BACK_ATOM + ['CB', 'OG1', 'CG2'],
                'V': BACK_ATOM + ['CB', 'CG1', 'CG2'],
                'I': BACK_ATOM + ['CB', 'CG1', 'CG2', 'CD'],
                'L': BACK_ATOM + ['CB', 'CG', 'CD1', 'CD2'],
                'N': BACK_ATOM + ['CB', 'CG', 'ND2', 'OD1'],
                'D': BACK_ATOM + ['CB', 'CG', 'OD1', 'OD2'],
                'M': BACK_ATOM + ['CB', 'CG', 'SD', 'CE'],
                'Q': BACK_ATOM + ['CB', 'CG', 'CD', 'NE2', 'OE1'],
                'E': BACK_ATOM + ['CB', 'CG', 'CD', 'OE1', 'OE2'],
                'K': BACK_ATOM + ['CB', 'CG', 'CD', 'CE', 'NZ'],
                'R': BACK_ATOM + ['CB', 'CG', 'CD', 'NE', 'CZ', 'NH1', 'NH2'],
                'F': BACK_ATOM + ['CB', 'CG', 'CD1', 'CE1', 'CZ', 'CD2',
                                  'CE2'],
                'Y': BACK_ATOM + ['CB', 'CG', 'CD1', 'CE1', 'CZ', 'CD2',
                                  'CE2', 'OH'],
                'H': BACK_ATOM + ['CB', 'CG', 'ND1', 'CE1', 'CD2', 'NE2'],
                'W': BACK_ATOM + ['CB', 'CG', 'CD1', 'NE1', 'CD2', 'CE2',
                                  'CE3', 'CZ3', 'CH2', 'CZ2'],
                'P': BACK_ATOM + ['CB', 'CG', 'CD']}

# Bond definition:
# Note that order is important
BACK_BOND = [['-C', 'N'], ['N', 'CA'], ['CA', 'C'], ['C', 'O']]
# X is for ACE special case
AA_BOND_DICT = {}
# Need to use a trick with unphysical bond
AA_BOND_DICT['X'] = [['CH3', 'O'], ['O', 'C']]
AA_BOND_DICT['G'] = BACK_BOND
AA_BOND_DICT['A'] = BACK_BOND + [['CA', 'CB']]
AA_BOND_DICT['S'] = BACK_BOND + [['CA', 'CB'], ['CB', 'OG']]
AA_BOND_DICT['C'] = BACK_BOND + [['CA', 'CB'], ['CB', 'SG']]
AA_BOND_DICT['T'] = BACK_BOND + [['CA', 'CB'], ['CB', 'OG1'], ['CB', 'CG2']]
AA_BOND_DICT['V'] = BACK_BOND + [['CA', 'CB'], ['CB', 'CG1'], ['CB', 'CG2']]
AA_BOND_DICT['I'] = BACK_BOND +\
    [['CA', 'CB'], ['CB', 'CG1'], ['CG1', 'CD'], ['CB', 'CG2']]
AA_BOND_DICT['L'] = BACK_BOND +\
    [['CA', 'CB'], ['CB', 'CG'], ['CG', 'CD1'], ['CG', 'CD2']]
AA_BOND_DICT['N'] = BACK_BOND +\
    [['CA', 'CB'], ['CB', 'CG'], ['CG', 'ND2'], ['CG', 'OD1']]
AA_BOND_DICT['D'] = BACK_BOND +\
    [['CA', 'CB'], ['CB', 'CG'], ['CG', 'OD1'], ['CG', 'OD2']]
AA_BOND_DICT['M'] = BACK_BOND +\
    [['CA', 'CB'], ['CB', 'CG'], ['CG', 'SD'], ['SD', 'CE']]
AA_BOND_DICT['Q'] = BACK_BOND +\
    [['CA', 'CB'], ['CB', 'CG'], ['CG', 'CD'], ['CD', 'NE2'],
     ['CD', 'OE1']]
AA_BOND_DICT['E'] = BACK_BOND +\
    [['CA', 'CB'], ['CB', 'CG'], ['CG', 'CD'], ['CD', 'OE1'],
     ['CD', 'OE2']]
AA_BOND_DICT['K'] = BACK_BOND +\
    [['CA', 'CB'], ['CB', 'CG'], ['CG', 'CD'], ['CD', 'CE'],
     ['CE', 'NZ']]
AA_BOND_DICT['R'] = BACK_BOND +\
    [['CA', 'CB'], ['CB', 'CG'], ['CG', 'CD'], ['CD', 'NE'],
     ['NE', 'CZ'], ['CZ', 'NH1'], ['CZ', 'NH2']]
AA_BOND_DICT['F'] = BACK_BOND +\
    [['CA', 'CB'], ['CB', 'CG'], ['CG', 'CD1'], ['CD1', 'CE1'],
     ['CE1', 'CZ'], ['CG', 'CD2'], ['CD2', 'CE2']]
AA_BOND_DICT['Y'] = BACK_BOND +\
    [['CA', 'CB'], ['CB', 'CG'], ['CG', 'CD1'], ['CD1', 'CE1'],
     ['CE1', 'CZ'], ['CZ', 'OH'], ['CG', 'CD2'], ['CD2', 'CE2']]
AA_BOND_DICT['H'] = BACK_BOND +\
    [['CA', 'CB'], ['CB', 'CG'], ['CG', 'ND1'], ['ND1', 'CE1'],
     ['CG', 'CD2'], ['CD2', 'NE2']]
AA_BOND_DICT['W'] = BACK_BOND +\
    [['CA', 'CB'], ['CB', 'CG'], ['CG', 'CD1'], ['CD1', 'NE1'],
     ['CG', 'CD2'], ['CD2', 'CE2'], ['CD2', 'CE3'], ['CE3', 'CZ3'],
     ['CZ3', 'CH2'], ['CH2', 'CZ2']]
AA_BOND_DICT['P'] = BACK_BOND +\
    [['CA', 'CB'], ['CB', 'CG'], ['CG', 'CD']]

# Distance, angle and dihedral angles parameters
BACK_DIST = [['N', 'CA', 1.46],
             ['CA', 'C', 1.52],
             ['C', 'O', 1.23],
             ['C', 'N', 1.29]]

BACK_ANGLE = [['N', 'CA', 'C', 110.9],
              ['CA', 'C', 'O', 122.0],
              ['CA', 'C', 'N', 110.9],
              ['CA', 'N', 'C', 121.3],
              ['N', 'C', 'O', 119.0]]  # Only for ACE-connexion

BACK_DIHE = [['N', 'CA', 'C', 'O', 0],
             ['N', 'CA', 'C', 'N', 180.0],
             ['CA', 'N', 'C', 'CA', 180.0],
             ['C', 'CA', 'N', 'C', -180.0],
             ['N', 'C', 'O', 'CH3', 180.0],  # Only for ACE-connexion
             ['CA', 'N', 'C', 'O', 0.0]]  # Only for ACE-connexion

DIST_DICT = {}
ANGLE_DICT = {}
DIHE_DICT = {}

# ACE X
DIST_DICT['X'] = [['CH3', 'O', 2.40], ['C', 'O', 1.23]]
ANGLE_DICT['X'] = [['CH3', 'O', 'C', 32.18]]
DIHE_DICT['X'] = BACK_DIHE


# Glycine
DIST_DICT['G'] = BACK_DIST
ANGLE_DICT['G'] = BACK_ANGLE
DIHE_DICT['G'] = BACK_DIHE

# Alanine
DIST_DICT['A'] = BACK_DIST + [['CA', 'CB', 1.52]]
ANGLE_DICT['A'] = BACK_ANGLE + [['CB', 'CA', 'N', 107.7]]
DIHE_DICT['A'] = BACK_DIHE + [['CB', 'CA', 'N', 'C', 56.3]]

# Serine
DIST_DICT['S'] = BACK_DIST + [['CA', 'CB', 1.52],
                              ['CB', 'OG', 1.42]]
ANGLE_DICT['S'] = BACK_ANGLE + [['CB', 'CA', 'N', 107.7],
                                ['CA', 'CB', 'OG', 110.8]]
DIHE_DICT['S'] = BACK_DIHE + [['CB', 'CA', 'N', 'C', 56.3],
                              ['N', 'CA', 'CB', 'OG', 69.4]]
# Cysteine
DIST_DICT['C'] = BACK_DIST + [['CA', 'CB', 1.52],
                              ['CB', 'SG', 1.81]]
ANGLE_DICT['C'] = BACK_ANGLE + [['CB', 'CA', 'N', 107.7],
                                ['CA', 'CB', 'SG', 110.8]]
DIHE_DICT['C'] = BACK_DIHE + [['CB', 'CA', 'N', 'C', 56.3],
                              ['N', 'CA', 'CB', 'SG', -173.8]]
# Threonine
DIST_DICT['T'] = BACK_DIST + [['CA', 'CB', 1.52],
                              ['CB', 'OG1', 1.42],
                              ['CB', 'CG2', 1.54]]
ANGLE_DICT['T'] = BACK_ANGLE + [['CB', 'CA', 'N', 107.7],
                                ['CA', 'CB', 'OG1', 110.6],
                                ['CA', 'CB', 'CG2', 116.3]]
DIHE_DICT['T'] = BACK_DIHE + [['CB', 'CA', 'N', 'C', 56.3],
                              ['N', 'CA', 'CB', 'OG1', -61.5],
                              ['N', 'CA', 'CB', 'CG2', 179.6]]

# Valine
DIST_DICT['V'] = BACK_DIST + [['CA', 'CB', 1.52],
                              ['CB', 'CG1', 1.54],
                              ['CB', 'CG2', 1.54]]
ANGLE_DICT['V'] = BACK_ANGLE + [['CB', 'CA', 'N', 107.7],
                                ['CA', 'CB', 'CG1', 110.6],
                                ['CA', 'CB', 'CG2', 116.3]]
DIHE_DICT['V'] = BACK_DIHE + [['CB', 'CA', 'N', 'C', 56.3],
                              ['N', 'CA', 'CB', 'CG1', -61.5],
                              ['N', 'CA', 'CB', 'CG2', 179.6]]

# Isoleucine
DIST_DICT['I'] = BACK_DIST + [['CA', 'CB', 1.52],
                              ['CB', 'CG1', 1.54],
                              ['CB', 'CG2', 1.54],
                              ['CG1', 'CD', 1.54]]
ANGLE_DICT['I'] = BACK_ANGLE + [['CB', 'CA', 'N', 107.7],
                                ['CA', 'CB', 'CG1', 110.6],
                                ['CA', 'CB', 'CG2', 116.3],
                                ['CB', 'CG1', 'CD', 116.3]]
DIHE_DICT['I'] = BACK_DIHE + [['CB', 'CA', 'N', 'C', 56.3],
                              ['N', 'CA', 'CB', 'CG1', -61.5],
                              ['N', 'CA', 'CB', 'CG2', 179.6],
                              ['CA', 'CB', 'CG1', 'CD', 179.6]]

# Isoleucine
DIST_DICT['L'] = BACK_DIST + [['CA', 'CB', 1.52],
                              ['CB', 'CG', 1.54],
                              ['CG', 'CD1', 1.54],
                              ['CG', 'CD2', 1.54]]
ANGLE_DICT['L'] = BACK_ANGLE + [['CB', 'CA', 'N', 107.7],
                                ['CA', 'CB', 'CG', 116.3],
                                ['CB', 'CG', 'CD1', 110.6],
                                ['CB', 'CG', 'CD2', 116.3]]
DIHE_DICT['L'] = BACK_DIHE + [['CB', 'CA', 'N', 'C', 56.3],
                              ['N', 'CA', 'CB', 'CG', -57.8],
                              ['CA', 'CB', 'CG', 'CD1', -61.5],
                              ['CA', 'CB', 'CG', 'CD2', 179.6]]

# Asparagine
DIST_DICT['N'] = BACK_DIST + [['CA', 'CB', 1.52],
                              ['CB', 'CG', 1.54],
                              ['CG', 'ND2', 1.29],
                              ['CG', 'OD1', 1.23]]
ANGLE_DICT['N'] = BACK_ANGLE + [['CB', 'CA', 'N', 107.7],
                                ['CA', 'CB', 'CG', 116.3],
                                ['CB', 'CG', 'ND2', 118.9],
                                ['CB', 'CG', 'OD1', 122.2]]
DIHE_DICT['N'] = BACK_DIHE + [['CB', 'CA', 'N', 'C', 56.3],
                              ['N', 'CA', 'CB', 'CG', -57.8],
                              ['CA', 'CB', 'CG', 'ND2', -78.2],
                              ['CA', 'CB', 'CG', 'OD1', 100.6]]

# Aspartic Acid
DIST_DICT['D'] = BACK_DIST + [['CA', 'CB', 1.52],
                              ['CB', 'CG', 1.54],
                              ['CG', 'OD1', 1.23],
                              ['CG', 'OD2', 1.23]]
ANGLE_DICT['D'] = BACK_ANGLE + [['CB', 'CA', 'N', 107.7],
                                ['CA', 'CB', 'CG', 116.3],
                                ['CB', 'CG', 'OD1', 118.9],
                                ['CB', 'CG', 'OD2', 122.2]]
DIHE_DICT['D'] = BACK_DIHE + [['CB', 'CA', 'N', 'C', 56.3],
                              ['N', 'CA', 'CB', 'CG', -177.0],
                              ['CA', 'CB', 'CG', 'OD1', 37.0],
                              ['CA', 'CB', 'CG', 'OD2', -140.7]]

# Methionine
DIST_DICT['M'] = BACK_DIST + [['CA', 'CB', 1.52],
                              ['CB', 'CG', 1.54],
                              ['CG', 'SD', 1.80],
                              ['SD', 'CE', 1.80]]
ANGLE_DICT['M'] = BACK_ANGLE + [['CB', 'CA', 'N', 107.7],
                                ['CA', 'CB', 'CG', 116.3],
                                ['CB', 'CG', 'SD', 118.9],
                                ['CG', 'SD', 'CE', 98.5]]
DIHE_DICT['M'] = BACK_DIHE + [['CB', 'CA', 'N', 'C', 56.3],
                              ['N', 'CA', 'CB', 'CG', -68.5],
                              ['CA', 'CB', 'CG', 'SD', -165.1],
                              ['CB', 'CG', 'SD', 'CE', -140.7]]

# Glutamine
DIST_DICT['Q'] = BACK_DIST + [['CA', 'CB', 1.52],
                              ['CB', 'CG', 1.54],
                              ['CG', 'CD', 1.54],
                              ['CD', 'NE2', 1.31],
                              ['CD', 'OE1', 1.22]]
ANGLE_DICT['Q'] = BACK_ANGLE + [['CB', 'CA', 'N', 107.7],
                                ['CA', 'CB', 'CG', 116.3],
                                ['CB', 'CG', 'CD', 118.9],
                                ['CG', 'CD', 'NE2', 121.1],
                                ['CG', 'CD', 'OE1', 120.0]]
DIHE_DICT['Q'] = BACK_DIHE + [['CB', 'CA', 'N', 'C', 56.3],
                              ['N', 'CA', 'CB', 'CG', -178.1],
                              ['CA', 'CB', 'CG', 'CD', -165.1],
                              ['CB', 'CG', 'CD', 'NE2', -0.9],
                              ['CB', 'CG', 'CD', 'OE1', 177.9]]

# Glutamic Acid
DIST_DICT['E'] = BACK_DIST + [['CA', 'CB', 1.52],
                              ['CB', 'CG', 1.54],
                              ['CG', 'CD', 1.54],
                              ['CD', 'OE1', 1.22],
                              ['CD', 'OE2', 1.22]]
ANGLE_DICT['E'] = BACK_ANGLE + [['CB', 'CA', 'N', 107.7],
                                ['CA', 'CB', 'CG', 116.3],
                                ['CB', 'CG', 'CD', 118.9],
                                ['CG', 'CD', 'OE1', 121.1],
                                ['CG', 'CD', 'OE2', 120.0]]
DIHE_DICT['E'] = BACK_DIHE + [['CB', 'CA', 'N', 'C', 56.3],
                              ['N', 'CA', 'CB', 'CG', -178.1],
                              ['CA', 'CB', 'CG', 'CD', -177.3],
                              ['CB', 'CG', 'CD', 'OE1', -0.9],
                              ['CB', 'CG', 'CD', 'OE2', 177.9]]

# Lysine
DIST_DICT['K'] = BACK_DIST + [['CA', 'CB', 1.52],
                              ['CB', 'CG', 1.54],
                              ['CG', 'CD', 1.54],
                              ['CD', 'CE', 1.54],
                              ['CE', 'NZ', 1.3]]
ANGLE_DICT['K'] = BACK_ANGLE + [['CB', 'CA', 'N', 107.7],
                                ['CA', 'CB', 'CG', 116.3],
                                ['CB', 'CG', 'CD', 118.9],
                                ['CG', 'CD', 'CE', 118.9],
                                ['CD', 'CE', 'NZ', 109.4]]
DIHE_DICT['K'] = BACK_DIHE + [['CB', 'CA', 'N', 'C', 56.3],
                              ['N', 'CA', 'CB', 'CG', -178.1],
                              ['CA', 'CB', 'CG', 'CD', -177.3],
                              ['CB', 'CG', 'CD', 'CE', 177.1],
                              ['CG', 'CD', 'CE', 'NZ', 177.9]]

# Arginine
DIST_DICT['R'] = BACK_DIST + [['CA', 'CB', 1.52],
                              ['CB', 'CG', 1.54],
                              ['CG', 'CD', 1.54],
                              ['CD', 'NE', 1.54],
                              ['NE', 'CZ', 1.3],
                              ['NH1', 'CZ', 1.3],
                              ['NH2', 'CZ', 1.3]]
ANGLE_DICT['R'] = BACK_ANGLE + [['CB', 'CA', 'N', 107.7],
                                ['CA', 'CB', 'CG', 116.3],
                                ['CB', 'CG', 'CD', 118.9],
                                ['CG', 'CD', 'NE', 118.9],
                                ['CD', 'NE', 'CZ', 125.3],
                                ['NE', 'CZ', 'NH1', 123.6],
                                ['NE', 'CZ', 'NH2', 123.6]]
DIHE_DICT['R'] = BACK_DIHE + [['CB', 'CA', 'N', 'C', 56.3],
                              ['N', 'CA', 'CB', 'CG', -178.1],
                              ['CA', 'CB', 'CG', 'CD', -177.3],
                              ['CB', 'CG', 'CD', 'NE', 177.1],
                              ['CG', 'CD', 'NE', 'CZ', 177.9],
                              ['CD', 'NE', 'CZ', 'NH1', 0.3],
                              ['CD', 'NE', 'CZ', 'NH2', -179.6]]

# Phenylalanine
DIST_DICT['F'] = BACK_DIST + [['CA', 'CB', 1.52],
                              ['CB', 'CG', 1.54],
                              ['CG', 'CD1', 1.4],
                              ['CG', 'CD2', 1.4],
                              ['CD1', 'CE1', 1.4],
                              ['CD2', 'CE2', 1.4],
                              ['CE1', 'CZ', 1.4]]
ANGLE_DICT['F'] = BACK_ANGLE + [['CB', 'CA', 'N', 107.7],
                                ['CA', 'CB', 'CG', 116.3],
                                ['CB', 'CG', 'CD1', 120.8],
                                ['CB', 'CG', 'CD2', 120.8],
                                ['CG', 'CD1', 'CE1', 120.4],
                                ['CG', 'CD2', 'CE2', 120.4],
                                ['CD1', 'CE1', 'CZ', 120.4]]
DIHE_DICT['F'] = BACK_DIHE + [['CB', 'CA', 'N', 'C', 56.3],
                              ['N', 'CA', 'CB', 'CG', -65.4],
                              ['CA', 'CB', 'CG', 'CD1', -78.1],
                              ['CA', 'CB', 'CG', 'CD2', 101.4],
                              ['CB', 'CG', 'CD1', 'CE1', 179.1],
                              ['CB', 'CG', 'CD2', 'CE2', 179.1],
                              ['CG', 'CD1', 'CE1', 'CZ', 0.1]]

# Tyrosine
DIST_DICT['Y'] = BACK_DIST + [['CA', 'CB', 1.52],
                              ['CB', 'CG', 1.54],
                              ['CG', 'CD1', 1.4],
                              ['CG', 'CD2', 1.4],
                              ['CD1', 'CE1', 1.4],
                              ['CD2', 'CE2', 1.4],
                              ['CE1', 'CZ', 1.4],
                              ['CZ', 'OH', 1.22]]
ANGLE_DICT['Y'] = BACK_ANGLE + [['CB', 'CA', 'N', 107.7],
                                ['CA', 'CB', 'CG', 116.3],
                                ['CB', 'CG', 'CD1', 120.8],
                                ['CB', 'CG', 'CD2', 120.8],
                                ['CG', 'CD1', 'CE1', 120.4],
                                ['CG', 'CD2', 'CE2', 120.4],
                                ['CD1', 'CE1', 'CZ', 120.4],
                                ['CE1', 'CZ', 'OH', 120.8]]
DIHE_DICT['Y'] = BACK_DIHE + [['CB', 'CA', 'N', 'C', 56.3],
                              ['N', 'CA', 'CB', 'CG', -65.4],
                              ['CA', 'CB', 'CG', 'CD1', -78.1],
                              ['CA', 'CB', 'CG', 'CD2', 101.4],
                              ['CB', 'CG', 'CD1', 'CE1', 179.1],
                              ['CB', 'CG', 'CD2', 'CE2', 179.1],
                              ['CG', 'CD1', 'CE1', 'CZ', 0.1],
                              ['CD1', 'CE1', 'CZ', 'OH', 179.1]]

# Histidine
DIST_DICT['H'] = BACK_DIST + [['CA', 'CB', 1.52],
                              ['CB', 'CG', 1.54],
                              ['CG', 'ND1', 1.4],
                              ['CG', 'CD2', 1.4],
                              ['ND1', 'CE1', 1.4],
                              ['CD2', 'NE2', 1.4]]
ANGLE_DICT['H'] = BACK_ANGLE + [['CB', 'CA', 'N', 107.7],
                                ['CA', 'CB', 'CG', 116.3],
                                ['CB', 'CG', 'ND1', 131.5],
                                ['CB', 'CG', 'CD2', 117.8],
                                ['CG', 'ND1', 'CE1', 105.4],
                                ['CG', 'CD2', 'NE2', 105.7]]
DIHE_DICT['H'] = BACK_DIHE + [['CB', 'CA', 'N', 'C', 56.3],
                              ['N', 'CA', 'CB', 'CG', -109.4],
                              ['CA', 'CB', 'CG', 'ND1', 107.6],
                              ['CA', 'CB', 'CG', 'CD2', -75.6],
                              ['CB', 'CG', 'ND1', 'CE1', 177.5],
                              ['CB', 'CG', 'CD2', 'NE2', 171.7]]

# Tryptophan
DIST_DICT['W'] = BACK_DIST + [['CA', 'CB', 1.52],
                              ['CB', 'CG', 1.54],
                              ['CG', 'CD1', 1.4],
                              ['CG', 'CD2', 1.4],
                              ['CD1', 'NE1', 1.4],
                              ['CD2', 'CE2', 1.4],
                              ['CD2', 'CE3', 1.4],
                              ['CE3', 'CZ3', 1.4],
                              ['CZ3', 'CH2', 1.4],
                              ['CH2', 'CZ2', 1.4]]
ANGLE_DICT['W'] = BACK_ANGLE + [['CB', 'CA', 'N', 107.7],
                                ['CA', 'CB', 'CG', 116.3],
                                ['CB', 'CG', 'CD1', 131.5],
                                ['CB', 'CG', 'CD2', 117.8],
                                ['CG', 'CD1', 'NE1', 105.4],
                                ['CG', 'CD2', 'CE2', 105.7],
                                ['CG', 'CD2', 'CE3', 131.5],
                                ['CD2', 'CE3', 'CZ3', 120.4],
                                ['CE3', 'CZ3', 'CH2', 120.4],
                                ['CZ3', 'CH2', 'CZ2', 120.4]]
DIHE_DICT['W'] = BACK_DIHE + [['CB', 'CA', 'N', 'C', 56.3],
                              ['N', 'CA', 'CB', 'CG', -109.4],
                              ['CA', 'CB', 'CG', 'CD1', 107.6],
                              ['CA', 'CB', 'CG', 'CD2', -75.6],
                              ['CB', 'CG', 'CD1', 'NE1', 177.5],
                              ['CB', 'CG', 'CD2', 'CE2', 180],
                              ['CB', 'CG', 'CD2', 'CE3', 0.0],
                              ['CG', 'CD2', 'CE3', 'CZ3', 180.0],
                              ['CD2', 'CE3', 'CZ3', 'CH2', 0.0],
                              ['CE3', 'CZ3', 'CH2', 'CZ2', 0]]

# Proline
DIST_DICT['P'] = BACK_DIST + [['CA', 'CB', 1.52],
                              ['CB', 'CG', 1.54],
                              ['CG', 'CD', 1.54]]
ANGLE_DICT['P'] = BACK_ANGLE + [['CB', 'CA', 'N', 101.9],
                                ['CA', 'CB', 'CG', 103.7],
                                ['CB', 'CG', 'CD', 103.3]]

DIHE_DICT['P'] = [['N', 'CA', 'C', 'O', 0],
                  ['N', 'CA', 'C', 'N', 180.0],
                  ['CA', 'N', 'C', 'CA', 180.0],
                  ['C', 'CA', 'N', 'C', -70.0],
                  ['N', 'C', 'O', 'CH3', 180.0],  # Only for ACE-connexion
                  ['CA', 'N', 'C', 'O', 0.0],
                  ['CB', 'CA', 'N', 'C', 168.6],
                  ['N', 'CA', 'CB', 'CG', 29.6],
                  ['CA', 'CB', 'CG', 'CD', -37.5]]

BLOSUM62 = {('W', 'F'): 1, ('L', 'R'): -2, ('S', 'P'): -1, ('V', 'T'): 0,
            ('Q', 'Q'): 5, ('N', 'A'): -2, ('Z', 'Y'): -2, ('W', 'R'): -3,
            ('Q', 'A'): -1, ('S', 'D'): 0, ('H', 'H'): 8, ('S', 'H'): -1,
            ('H', 'D'): -1, ('L', 'N'): -3, ('W', 'A'): -3, ('Y', 'M'): -1,
            ('G', 'R'): -2, ('Y', 'I'): -1, ('Y', 'E'): -2, ('B', 'Y'): -3,
            ('Y', 'A'): -2, ('V', 'D'): -3, ('B', 'S'): 0, ('Y', 'Y'): 7,
            ('G', 'N'): 0, ('E', 'C'): -4, ('Y', 'Q'): -1, ('Z', 'Z'): 4,
            ('V', 'A'): 0, ('C', 'C'): 9, ('M', 'R'): -1, ('V', 'E'): -2,
            ('T', 'N'): 0, ('P', 'P'): 7, ('V', 'I'): 3, ('V', 'S'): -2,
            ('Z', 'P'): -1, ('V', 'M'): 1, ('T', 'F'): -2, ('V', 'Q'): -2,
            ('K', 'K'): 5, ('P', 'D'): -1, ('I', 'H'): -3, ('I', 'D'): -3,
            ('T', 'R'): -1, ('P', 'L'): -3, ('K', 'G'): -2, ('M', 'N'): -2,
            ('P', 'H'): -2, ('F', 'Q'): -3, ('Z', 'G'): -2, ('X', 'L'): -1,
            ('T', 'M'): -1, ('Z', 'C'): -3, ('X', 'H'): -1, ('D', 'R'): -2,
            ('B', 'W'): -4, ('X', 'D'): -1, ('Z', 'K'): 1, ('F', 'A'): -2,
            ('Z', 'W'): -3, ('F', 'E'): -3, ('D', 'N'): 1, ('B', 'K'): 0,
            ('X', 'X'): -1, ('F', 'I'): 0, ('B', 'G'): -1, ('X', 'T'): 0,
            ('F', 'M'): 0, ('B', 'C'): -3, ('Z', 'I'): -3, ('Z', 'V'): -2,
            ('S', 'S'): 4, ('L', 'Q'): -2, ('W', 'E'): -3, ('Q', 'R'): 1,
            ('N', 'N'): 6, ('W', 'M'): -1, ('Q', 'C'): -3, ('W', 'I'): -3,
            ('S', 'C'): -1, ('L', 'A'): -1, ('S', 'G'): 0, ('L', 'E'): -3,
            ('W', 'Q'): -2, ('H', 'G'): -2, ('S', 'K'): 0, ('Q', 'N'): 0,
            ('N', 'R'): 0, ('H', 'C'): -3, ('Y', 'N'): -2, ('G', 'Q'): -2,
            ('Y', 'F'): 3, ('C', 'A'): 0, ('V', 'L'): 1, ('G', 'E'): -2,
            ('G', 'A'): 0, ('K', 'R'): 2, ('E', 'D'): 2, ('Y', 'R'): -2,
            ('M', 'Q'): 0, ('T', 'I'): -1, ('C', 'D'): -3, ('V', 'F'): -1,
            ('T', 'A'): 0, ('T', 'P'): -1, ('B', 'P'): -2, ('T', 'E'): -1,
            ('V', 'N'): -3, ('P', 'G'): -2, ('M', 'A'): -1, ('K', 'H'): -1,
            ('V', 'R'): -3, ('P', 'C'): -3, ('M', 'E'): -2, ('K', 'L'): -2,
            ('V', 'V'): 4, ('M', 'I'): 1, ('T', 'Q'): -1, ('I', 'G'): -4,
            ('P', 'K'): -1, ('M', 'M'): 5, ('K', 'D'): -1, ('I', 'C'): -1,
            ('Z', 'D'): 1, ('F', 'R'): -3, ('X', 'K'): -1, ('Q', 'D'): 0,
            ('X', 'G'): -1, ('Z', 'L'): -3, ('X', 'C'): -2, ('Z', 'H'): 0,
            ('B', 'L'): -4, ('B', 'H'): 0, ('F', 'F'): 6, ('X', 'W'): -2,
            ('B', 'D'): 4, ('D', 'A'): -2, ('S', 'L'): -2, ('X', 'S'): 0,
            ('F', 'N'): -3, ('S', 'R'): -1, ('W', 'D'): -4, ('V', 'Y'): -1,
            ('W', 'L'): -2, ('H', 'R'): 0, ('W', 'H'): -2, ('H', 'N'): 1,
            ('W', 'T'): -2, ('T', 'T'): 5, ('S', 'F'): -2, ('W', 'P'): -4,
            ('L', 'D'): -4, ('B', 'I'): -3, ('L', 'H'): -3, ('S', 'N'): 1,
            ('B', 'T'): -1, ('L', 'L'): 4, ('Y', 'K'): -2, ('E', 'Q'): 2,
            ('Y', 'G'): -3, ('Z', 'S'): 0, ('Y', 'C'): -2, ('G', 'D'): -1,
            ('B', 'V'): -3, ('E', 'A'): -1, ('Y', 'W'): 2, ('E', 'E'): 5,
            ('Y', 'S'): -2, ('C', 'N'): -3, ('V', 'C'): -1, ('T', 'H'): -2,
            ('P', 'R'): -2, ('V', 'G'): -3, ('T', 'L'): -1, ('V', 'K'): -2,
            ('K', 'Q'): 1, ('R', 'A'): -1, ('I', 'R'): -3, ('T', 'D'): -1,
            ('P', 'F'): -4, ('I', 'N'): -3, ('K', 'I'): -3, ('M', 'D'): -3,
            ('V', 'W'): -3, ('W', 'W'): 11, ('M', 'H'): -2, ('P', 'N'): -2,
            ('K', 'A'): -1, ('M', 'L'): 2, ('K', 'E'): 1, ('Z', 'E'): 4,
            ('X', 'N'): -1, ('Z', 'A'): -1, ('Z', 'M'): -1, ('X', 'F'): -1,
            ('K', 'C'): -3, ('B', 'Q'): 0, ('X', 'B'): -1, ('B', 'M'): -3,
            ('F', 'C'): -2, ('Z', 'Q'): 3, ('X', 'Z'): -1, ('F', 'G'): -3,
            ('B', 'E'): 1, ('X', 'V'): -1, ('F', 'K'): -3, ('B', 'A'): -2,
            ('X', 'R'): -1, ('D', 'D'): 6, ('W', 'G'): -2, ('Z', 'F'): -3,
            ('S', 'Q'): 0, ('W', 'C'): -2, ('W', 'K'): -3, ('H', 'Q'): 0,
            ('L', 'C'): -1, ('W', 'N'): -4, ('S', 'A'): 1, ('L', 'G'): -4,
            ('W', 'S'): -3, ('S', 'E'): 0, ('H', 'E'): 0, ('S', 'I'): -2,
            ('H', 'A'): -2, ('S', 'M'): -1, ('Y', 'L'): -1, ('Y', 'H'): 2,
            ('Y', 'D'): -3, ('E', 'R'): 0, ('X', 'P'): -2, ('G', 'G'): 6,
            ('G', 'C'): -3, ('E', 'N'): 0, ('Y', 'T'): -2, ('Y', 'P'): -3,
            ('T', 'K'): -1, ('A', 'A'): 4, ('P', 'Q'): -1, ('T', 'C'): -1,
            ('V', 'H'): -3, ('T', 'G'): -2, ('I', 'Q'): -3, ('Z', 'T'): -1,
            ('C', 'R'): -3, ('V', 'P'): -2, ('P', 'E'): -1, ('M', 'C'): -1,
            ('K', 'N'): 0, ('I', 'I'): 4, ('P', 'A'): -1, ('M', 'G'): -3,
            ('T', 'S'): 1, ('I', 'E'): -3, ('P', 'M'): -2, ('M', 'K'): -1,
            ('I', 'A'): -1, ('P', 'I'): -3, ('R', 'R'): 5, ('X', 'M'): -1,
            ('L', 'I'): 2, ('X', 'I'): -1, ('Z', 'B'): 1, ('X', 'E'): -1,
            ('Z', 'N'): 0, ('X', 'A'): 0, ('B', 'R'): -1, ('B', 'N'): 3,
            ('F', 'D'): -3, ('X', 'Y'): -1, ('Z', 'R'): 0, ('F', 'H'): -1,
            ('B', 'F'): -3, ('F', 'L'): 0, ('X', 'Q'): -1, ('B', 'B'): 4}


class Coor:
    """ Topologie base on coordinates like pdb or gro.

    The coor object containt a dictionnary of atoms indexed
    on the atom num and the crystal packing info.


    :param atom_dict: dictionnary of atom
    :type atom_dict: dict

    :param crystal_pack: crystal packing
    :type crystal_pack: str

    **Atom dictionnary parameters**

    :param field: pdb field
    :type field: str

    :param num: atom number
    :type num: int

    :param name: atom name
    :type name: str

    :param alter_loc: atom number
    :type alter_loc: str

    :param res_name: residue name (3 letters)
    :type res_name: str

    :param chain: chain ID
    :type chain: str

    :param res_num: residue number (based on pdb file)
    :type res_num: int

    :param uniq_resid: unique residue number
    :type uniq_resid: int

    :param insert_res: atom number
    :type insert_res: str

    :param xyz: coordinate
    :type x: numpy array

    :param occ: occupation
    :type occ: float

    :param beta: beta flactor
    :type beta: float


    .. note::
        The atom num index in the dictionnary, is not the same as the
        ``atom_num`` field of the dictionnary.

    .. note::
        Files necessary for testing : ../test/input/1y0m.pdb,
        ../test/input/1rxz.pdb and ../test/input/4n1m.pdb.
        To do the unitary test, execute pdb_mani.py (-v for verbose mode)

    .. todo::
        Add an atom class ?

    """

    def __init__(self, coor_in=None, pdb_lines=None):
        self.atom_dict = dict()
        self.crystal_pack = None
        self._num = None
        self.atom_dict = {}
        self.title = None

        if coor_in is not None:
            self.read_file(coor_in)
        elif pdb_lines is not None:
            self.parse_pdb_lines(pdb_lines)

    @property
    def num(self):
        return len(self.atom_dict)

    @property
    def view(self):
        """ Return a `nglview` object to be view in
        a jupyter notebook.

        Example:

        >>> TEST_OUT = str(getfixture('tmpdir'))
        >>> prot_coor = Coor()
        >>> prot_coor.get_PDB('3EAM', os.path.join(TEST_OUT, '3eam.pdb'))
        Succeed to read file ...3eam.pdb ,  13505 atoms found
        >>> view = prot_coor.view #doctest: +SKIP
        >>> view #doctest: +SKIP
        """
        import nglview as nv

        struct_str = nv.TextStructure(self.get_structure_string())
        view = nv.NGLWidget(struct_str)

        return view

    def get_PDB(self, pdb_ID, out_file=None, check_file_out=True):
        """Get a pdb file from the PDB using its ID
        and return a Coor object.

        :param pdb_ID: Protein Data Bank structure ID
        :type pdb_ID: str

        :param out_file: path of the pdb file to save
        :type out_file: str, optional, default=None

        :param check_file_out: flag to check or not if
            file has already been created.
            If the file is present then the command break.
        :type check_file_out: bool, optional, default=True


        :Example:

        >>> show_log()
        >>> TEST_OUT = str(getfixture('tmpdir'))
        >>> prot_coor = Coor()
        >>> prot_coor.get_PDB('3EAM', os.path.join(TEST_OUT, '3eam.pdb'))\
        #doctest: +ELLIPSIS
        Succeed to read file ...3eam.pdb ,  13505 atoms found
        """

        # Define output file:
        if out_file is None:
            out_file = '{}.pdb'.format(pdb_ID)

        if check_file_out and os_command.check_file_and_create_path(out_file):
            logger.info("PDB file {} already exist, file not saved".format(
                out_file))
            self.read_pdb(out_file)
            return

        # Get the pdb file from the PDB:
        urllib.request.urlretrieve(
            'http://files.rcsb.org/download/{}.pdb'.format(pdb_ID), out_file)

        self.read_pdb(out_file)

    def read_file(self, file_in):
        """Read a pdb file and return atom informations as a dictionnary
        indexed on the atom num. The fonction can also read pqr files if
        specified with ``pqr_format = True``,
        it will only change the column format of beta and occ factors.

        :param pdb_in: path of the pdb file to read
        :type pdb_in: str

        :param pqr_format: Flag for .pqr file format reading.
        :type pqr_format: bool, default=False

        :Example:

        >>> prot_coor = Coor()
        >>> prot_coor.read_file(os.path.join(TEST_PATH, '1y0m.pdb'))\
        #doctest: +ELLIPSIS
        Succeed to read file ...1y0m.pdb ,  648 atoms found
        >>> prot_coor.read_file(os.path.join(TEST_PATH, '1y0m.gro'))\
        #doctest: +ELLIPSIS
        Succeed to read file ...1y0m.gro ,  648 atoms found

        """

        file_lines = open(file_in)
        lines = file_lines.readlines()
        if str(file_in).endswith('.gro'):
            self.parse_gro_lines(lines)
        elif str(file_in).endswith('.pqr'):
            self.parse_pdb_lines(lines, pqr_format=True)
        elif str(file_in).endswith('.pdb'):
            self.parse_pdb_lines(lines, pqr_format=False)
        else:
            logger.warning('File name doesn\'t finish with .pdb'
                           'read it as .pdb anyway')
            self.parse_pdb_lines(lines, pqr_format=False)

        logger.info("Succeed to read file {} ,  {} atoms found".format(
            os.path.relpath(file_in), self.num))

    def read_pdb(self, pdb_in, pqr_format=False):
        """Read a pdb file and return atom informations as a dictionnary
        indexed on the atom num. The fonction can also read pqr files if
        specified with ``pqr_format = True``,
        it will only change the column format of beta and occ factors.

        :param pdb_in: path of the pdb file to read
        :type pdb_in: str

        :param pqr_format: Flag for .pqr file format reading.
        :type pqr_format: bool, default=False

        :Example:

        >>> prot_coor = Coor()
        >>> prot_coor.read_pdb(os.path.join(TEST_PATH, '1y0m.pdb'))\
        #doctest: +ELLIPSIS
        Succeed to read file ...1y0m.pdb ,  648 atoms found

        """

        file_in = open(pdb_in)
        lines = file_in.readlines()
        self.parse_pdb_lines(lines, pqr_format=pqr_format)

        logger.info("Succeed to read file {} ,  {} atoms found".format(
            os.path.relpath(pdb_in), self.num))

    def parse_pdb_lines(self, pdb_lines, pqr_format=False):
        """Parse the pdb lines and return atom informations as a dictionnary
        indexed on the atom num. The fonction can also read pqr files if
        specified with ``pqr_format = True``,
        it will only change the column format of beta and occ factors.

        :param pdb_lines: lines to parse
        :type pdb_lines: list of str

        :Example:

        >>> prot_coor = Coor()
        >>> f = open(os.path.join(TEST_PATH, '1y0m.pdb'))
        >>> lines = f.readlines()
        >>> prot_coor.parse_pdb_lines(lines)
        >>> prot_coor.num
        648

        """

        atom_index = 0
        uniq_resid = -1
        old_res_num = -1

        for line in pdb_lines:
            # print(line)
            if line.startswith("CRYST1"):
                self.crystal_pack = line
            if line.startswith('ATOM') or line.startswith("HETATM"):

                field = line[:6].strip()
                atom_num = int(line[6:11])
                atom_name = line[12:16].strip()

                res_name = line[17:20].strip()
                chain = line[21]
                res_num = int(line[22:26])
                insert_res = line[26:27]
                xyz = np.array([float(line[30:38]),
                                float(line[38:46]),
                                float(line[46:54])])
                if pqr_format:
                    alter_loc = ""
                    res_name = line[16:20].strip()
                    occ, beta = line[54:62].strip(), line[62:70].strip()
                    elem_symbol = ''
                else:
                    alter_loc = line[16:17]
                    res_name = line[17:20].strip()
                    occ, beta = line[54:60].strip(), line[60:66].strip()
                    elem_symbol = line[76:78]

                if occ == "":
                    occ = 0.0
                else:
                    occ = float(occ)

                if beta == "":
                    beta = 0.0
                else:
                    beta = float(beta)

                if res_num != old_res_num:
                    uniq_resid += 1
                    old_res_num = res_num

                atom = {"field": field,
                        "num": atom_num,
                        "name": atom_name,
                        "alter_loc": alter_loc,
                        "res_name": res_name,
                        "chain": chain,
                        "res_num": res_num,
                        "uniq_resid": uniq_resid,
                        "insert_res": insert_res,
                        "xyz": xyz,
                        "occ": occ,
                        "beta": beta,
                        "elem": elem_symbol}

                self.atom_dict[atom_index] = atom
                atom_index += 1

        logger.debug("Succeeded to parse lines. %d atoms found" % atom_index)
        return

    def parse_gro_lines(self, gro_lines):
        """Parse a gro file and return atom informations as a dictionnary
        indexed on the atom num. 

        :param gro_in: lines to parse
        :type gro_in: list of str

        :Example:

        >>> prot_coor = Coor()
        >>> f = open(os.path.join(TEST_PATH, '1y0m.gro'))
        >>> lines = f.readlines()
        >>> prot_coor.parse_gro_lines(lines)
        >>> prot_coor.num
        648

        """

        line_len = len(gro_lines)

        atom_index = 0
        uniq_resid = -1
        old_res_num = -1

        for i, line in enumerate(gro_lines):
            #print(line)
            if i == 0:
                self.title = line
            elif i ==1:
                num = int(line.strip())
            elif i == line_len - 1:
                self.crystal_pack = line
            elif i>= 2:
                # "%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f"
                res_num = int(line[:5])
                res_name = line[5:10].strip()
                atom_name = line[10:15].strip()
                atom_num = int(line[15:20])
                xyz = np.array([float(line[20:28]) * 10,
                                float(line[28:36]) * 10,
                                float(line[36:44]) * 10])
                occ = 0.0
                beta = 0.0
                field = 'ATOM'
                alter_loc = ''
                chain = ''
                insert_res = ''
                elem_symbol = ''

                if res_num != old_res_num:
                    uniq_resid += 1
                    old_res_num = res_num

                atom = {"field": field,
                        "num": atom_num,
                        "name": atom_name,
                        "alter_loc": alter_loc,
                        "res_name": res_name,
                        "chain": chain,
                        "res_num": res_num,
                        "uniq_resid": uniq_resid,
                        "insert_res": insert_res,
                        "xyz": xyz,
                        "occ": occ,
                        "beta": beta,
                        "elem": elem_symbol}

                self.atom_dict[atom_index] = atom
                atom_index += 1

        if num != len(self.atom_dict):
            logger.warning('Mismatch with atom number in gro file')

        logger.debug(f"Succeeded to parse lines. {atom_index} atoms found")
        return

    def get_structure_string(self):
        """Return a coor object as a pdb string.

        :Example:

        >>> prot_coor = Coor()
        >>> prot_coor.read_pdb(os.path.join(TEST_PATH, '1y0m.pdb'))\
        #doctest: +ELLIPSIS
        Succeed to read file ...1y0m.pdb ,  648 atoms found
        >>> pdb_str = prot_coor.get_structure_string()
        >>> print('Number of caracters: {}'.format(len(pdb_str)))
        Number of caracters: 51264

        """

        str_out = ""
        if self.crystal_pack is not None:
            str_out += self.cryst_convert(format_out='pdb')

        for atom_num, atom in sorted(self.atom_dict.items()):
            # Atom name should start a column 14, with the type of atom ex:
            #   - with atom type 'C': ' CH3'
            # for 2 letters atom type, it should start at coulumn 13 ex:
            #   - with atom type 'FE': 'FE1'
            name = atom["name"]
            if len(name) <= 3 and name[0] in ['C', 'H', 'O', 'N', 'S', 'P']:
                name = " " + name

            str_out += "{:6s}{:5d} {:4s}{:1s}{:3s} {:1s}{:4d}{:1s}"\
                       "   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}"\
                       "          {:2s}\n".format(
                            atom["field"],
                            atom["num"],
                            name,
                            atom["alter_loc"],
                            atom["res_name"],
                            atom["chain"],
                            atom["res_num"],
                            atom["insert_res"],
                            atom["xyz"][0],
                            atom["xyz"][1],
                            atom["xyz"][2],
                            atom["occ"],
                            atom["beta"],
                            atom['elem'])

        return str_out

    def cryst_convert(self, format_out='pdb'):
        """
        PDB format:
        https://www.wwpdb.org/documentation/file-format-content/format33/sect8.html

        Gro to pdb:
        https://mailman-1.sys.kth.se/pipermail/gromacs.org_gmx-users/2008-May/033944.html

        https://en.wikipedia.org/wiki/Fractional_coordinates

        >>> prot_coor = Coor()
        >>> prot_coor.read_file(os.path.join(TEST_PATH, '1y0m.gro'))\
        #doctest: +ELLIPSIS
        Succeed to read file ...1y0m.gro ,  648 atoms found
        >>> prot_coor.cryst_convert(format_out='pdb')
        'CRYST1   28.748   30.978   29.753  90.00  92.12  90.00 P 1           1\\n'
        >>> prot_coor = Coor()
        >>> prot_coor.read_file(os.path.join(TEST_PATH, '1y0m.pdb'))\
        #doctest: +ELLIPSIS
        Succeed to read file ...1y0m.pdb ,  648 atoms found
        >>> prot_coor.cryst_convert(format_out='gro')
        '   2.87480   3.09780   2.97326   0.00000   0.00000   0.00000   0.00000  -0.11006   0.00000\\n'
        """
        line = self.crystal_pack
        if line.startswith("CRYST1"):
            format_in = 'pdb'
            a = float(line[6:15])
            b = float(line[15:24])
            c = float(line[24:33])
            alpha = float(line[33:40])
            beta = float(line[40:47])
            gamma = float(line[47:54])
            sGroup = line[56:66]
            z = int(line[67:70])
        else:
            format_in = 'gro'
            line_split = line.split()
            #  v1(x) v2(y) v3(z) v1(y) v1(z) v2(x) v2(z) v3(x) v3(y)
            if len(line_split) == 3:
                v1 = np.array([float(line_split[0]), 0., 0.])
                v2 = np.array([0., float(line_split[1]), 0.])
                v3 = np.array([0., 0., float(line_split[2])])
            elif len(line_split) == 9:
                v1 = np.array([float(line_split[0]), float(line_split[3]), float(line_split[4])])
                v2 = np.array([float(line_split[5]), float(line_split[1]), float(line_split[6])])
                v3 = np.array([float(line_split[7]), float(line_split[8]), float(line_split[2])])

        # Convert:
        if format_out == 'pdb':
            if format_in == 'gro':
                a = sum(v1**2)**0.5 * 10
                b = sum(v2**2)**0.5 * 10
                c = sum(v3**2)**0.5 * 10
                alpha = np.rad2deg(Coor.angle_vec(v2, v3))
                beta = np.rad2deg(Coor.angle_vec(v1, v3))
                gamma = np.rad2deg(Coor.angle_vec(v1, v2))
                # Following is wrong, to check !!!
                sGroup = '1'
                z = 1
            new_line = "CRYST1{:9.3f}{:9.3f}{:9.3f}{:7.2f}{:7.2f}{:7.2f} P {:9} {:3d}\n".format(
                    a, b, c, alpha, beta, gamma, sGroup, z)
        elif format_out == 'gro':
            if format_in == 'pdb':
                alpha = np.deg2rad(alpha)
                beta = np.deg2rad(beta)
                gamma = np.deg2rad(gamma)
                v1 = [a / 10, 0., 0.]
                v2 = [b * cos(gamma) / 10, b * sin(gamma) / 10, 0.]
                v = (1.0 - cos(alpha)**2 - cos(beta)**2 - cos(gamma)**2 +
                    2.0 * cos(alpha) * cos(beta) * cos(gamma))**0.5 *\
                    a * b * c
                v3 = [c * cos(beta) / 10,
                      (c / sin(gamma)) * (cos(alpha) - cos(beta) * cos(gamma)) / 10,
                      v / (a * b * sin(gamma)) /10]

            new_line = "{:10.5f}{:10.5f}{:10.5f}{:10.5f}{:10.5f}{:10.5f}{:10.5f}{:10.5f}{:10.5f}\n".format(
                    v1[0],v2[1],v3[2],v1[1],v1[2],v2[0],v2[2],v3[0],v3[1])

        return(new_line)

    def get_gro_structure_string(self):
        """Return a coor object as a pdb string.

        :Example:

        >>> prot_coor = Coor()
        >>> prot_coor.read_pdb(os.path.join(TEST_PATH, '1y0m.pdb'))\
        #doctest: +ELLIPSIS
        Succeed to read file ...1y0m.pdb ,  648 atoms found
        >>> gro_str = prot_coor.get_gro_structure_string()
        >>> print('Number of caracters: {}'.format(len(gro_str)))
        Number of caracters: 29283

        """

        str_out = ""
        if self.title is not None:
            str_out += self.title
        else:
            str_out += "Create with pdb_manip_py\n"
        str_out += f"{self.num:6}\n"
        for atom_num, atom in sorted(self.atom_dict.items()):

            str_out += "{:5d}{:5s}{:>5s}{:5d}"\
                       "{:8.3f}{:8.3f}{:8.3f}\n".format(
                            atom["res_num"],
                            atom["res_name"],
                            atom["name"],
                            atom["num"],
                            atom["xyz"][0]/10,
                            atom["xyz"][1]/10,
                            atom["xyz"][2]/10)
        if self.crystal_pack is not None:
            str_out += self.cryst_convert(format_out='gro')
        #print(str_out)
        return str_out

    def write_pdb(self, pdb_out, check_file_out=True):
        """Write a pdb file.

        :param pdb_out: path of the pdb file to write
        :type pdb_out: str

        :param check_file_out: flag to check or not if
            file has already been created.
            If the file is present then the command break.
        :type check_file_out: bool, optional, default=True

        :Example:

        >>> TEST_OUT = str(getfixture('tmpdir'))
        >>> prot_coor = Coor(os.path.join(TEST_PATH, '1y0m.pdb'))\
        #doctest: +ELLIPSIS
        Succeed to read file ...1y0m.pdb ,  648 atoms found
        >>> prot_coor.write_pdb(os.path.join(TEST_OUT, 'tmp.pdb'))\
        #doctest: +ELLIPSIS
        Succeed to save file ...tmp.pdb

        """

        if check_file_out and os_command.check_file_and_create_path(pdb_out):
            logger.info("PDB file {} already exist, file not saved".format(
                pdb_out))
            return

        filout = open(pdb_out, 'w')
        filout.write(self.get_structure_string())

        logger.info("Succeed to save file %s" % os.path.relpath(pdb_out))
        return

    def get_aa_seq(self):
        """Get the amino acid sequence from a coor object.

        :return: dictionnary of sequence indexed by the chain ID
        :rtype: dict

        :Example:

        >>> prot_coor = Coor(os.path.join(TEST_PATH, '1y0m.pdb'))\
        #doctest: +ELLIPSIS
        Succeed to read file ...1y0m.pdb ,  648 atoms found
        >>> prot_coor.get_aa_seq()
        {'A': 'TFKSAVKALFDYKAQREDELTFTKSAIIQNVEKQDGGWWRGDYGGKKQLWFPSNYVEEMIN'}

        .. warning::
            If atom chains are not arranged sequentialy (A,A,A,B,B,A,A,A ...),
            the first atom seq will be overwritten by the last one.

        """

        # Get CA atoms
        CA_index_list = self.get_index_selection({"name": ["CA"]})

        seq = ""
        seq_dict = {}
        chain_first = self.atom_dict[CA_index_list[0]]['chain']

        for index in sorted(CA_index_list):
            loop_atom = self.atom_dict[index]

            if loop_atom['chain'] != chain_first:
                seq_dict[chain_first] = seq
                seq = ""
                chain_first = loop_atom['chain']

            seq = seq + AA_DICT[loop_atom['res_name']]

        seq_dict[chain_first] = seq

        return seq_dict

    def get_aa_num(self):
        """Get the amino acid number of a coor object.

        :return: Number of residues
        :rtype: int

        :Example:

        >>> prot_coor = Coor(os.path.join(TEST_PATH, '1y0m.pdb'))\
        #doctest: +ELLIPSIS
        Succeed to read file ...1y0m.pdb ,  648 atoms found
        >>> prot_coor.get_aa_num()
        61

        .. note::
            Only count Ca atoms, this may not be the best choice ?

        """

        CA_index_list = self.get_index_selection({"name": ["CA"]})

        return len(CA_index_list)

    def get_array(self, field='xyz', index_list=None):
        """ Convert atom dict as a numpy array.

        :param field: field to extract
        :type field: str (Default='xyz')

        :param index_list: list of index to extract
        :type index_list: list (Default=None)

        :return: coordinates array
        :rtype: np.array

        :Example:

        >>> prot_coor = Coor(os.path.join(TEST_PATH, '1y0m.pdb'))\
        #doctest: +ELLIPSIS
        Succeed to read file ...1y0m.pdb ,  648 atoms found
        >>> coor_array = prot_coor.get_array()
        >>> print(coor_array[1:4])
        [[-0.971  9.213 12.734]
         [-0.185  7.901 12.631]
         [ 0.456  7.524 13.61 ]]
        >>> coor_array = prot_coor.get_array(field='name')
        >>> print(coor_array[1:4])
        ['CA' 'C' 'O']
        """

        if index_list is None:
            return(np.array([atom[field] for key, atom in sorted(
                self.atom_dict.items())]))
        else:
            return(np.array([self.atom_dict[index][field] for index in index_list]))

    def change_pdb_field(self, change_dict):
        """Change all atom field of a coor object,
        the change is based on the change_dict dictionnary.

        :param change_dict: change ditionnay eg. {"chain" : "A"}
        :type change_dict: dict

        :Example:

        >>> prot_coor = Coor(os.path.join(TEST_PATH, '1y0m.pdb'))\
        #doctest: +ELLIPSIS
        Succeed to read file ...1y0m.pdb ,  648 atoms found
        >>> prot_coor.get_aa_seq()
        {'A': 'TFKSAVKALFDYKAQREDELTFTKSAIIQNVEKQDGGWWRGDYGGKKQLWFPSNYVEEMIN'}
        >>> prot_coor.change_pdb_field(change_dict = {"chain" : "B"})\
        #doctest: +ELLIPSIS
        <...Coor object at ...>
        >>> prot_coor.get_aa_seq()
        {'B': 'TFKSAVKALFDYKAQREDELTFTKSAIIQNVEKQDGGWWRGDYGGKKQLWFPSNYVEEMIN'}

        """

        for atom_num, atom in self.atom_dict.items():
            for change, val in change_dict.items():
                atom[change] = val

        return self

    def change_index_pdb_field(self, index_list, change_dict):
        """Change all atom field of a part of coor object defined by ``index``,
        the change is based on the change_dict dictionnary.

        :param index_list: list of atom index to change
        :type index_list: list

        :param change_dict: change ditionnay eg. {"chain" : "A"}
        :type change_dict: dict

        :Example:

        >>> prot_coor = Coor(os.path.join(TEST_PATH, '1y0m.pdb'))\
        #doctest: +ELLIPSIS
        Succeed to read file ...1y0m.pdb ,  648 atoms found
        >>> res_826_852 = prot_coor.get_index_selection({'res_num' :
        ... range(826,852)})
        >>> prot_coor.change_index_pdb_field(index_list = res_826_852,
        ... change_dict = {"chain" : "B"}) #doctest: +ELLIPSIS
        <...Coor object at ...>
        >>> prot_seq = prot_coor.get_aa_seq()
        >>> prot_seq == {'A': 'TFKSAVKALFDYKAQREDELTFTKSAIIQNVEKQD',
        ... 'B': 'GGWWRGDYGGKKQLWFPSNYVEEMIN'}
        True

        """

        for atom_num in index_list:
            for change, val in change_dict.items():
                self.atom_dict[atom_num][change] = val

        return self

    def select_part_dict(self, selec_dict):
        """Select atom of a coor object defined,
        the selection is based on the change_dict dictionnary.
        Return a new coor object.

        :param selec_dict: change ditionnay eg. {"chain" : ["A","G"]}
        :type selec_dict: dict

        :return: a new coor object
        :rtype: coor

        :Example:

        >>> prot_coor = Coor(os.path.join(TEST_PATH, '1y0m.pdb'))\
        #doctest: +ELLIPSIS
        Succeed to read file ...1y0m.pdb ,  648 atoms found
        >>> prot_coor.get_aa_num()
        61
        >>> prot_20_coor = prot_coor.select_part_dict(selec_dict = \
{'res_num': list(range(791,800))})
        >>> prot_20_coor.get_aa_seq()
        {'A': 'TFKSAVKAL'}
        >>> prot_20_coor.get_aa_num()
        9
        >>> prot_N_atom = prot_coor.select_part_dict(selec_dict = {'name'\
: ['ZN']})
        >>> # WARNING using selec_dict = {'name' : 'ZN'} will
        >>> # give you 61 residues !!
        >>> print(prot_N_atom.num)
        0
        >>> # Select only protein atoms
        >>> print(prot_coor.num)
        648
        >>> prot_only = prot_coor.select_part_dict(selec_dict = \
{'res_name': PROTEIN_AA})
        >>> print(prot_only.num)
        526
        """

        coor_out = Coor()

        for atom_num, atom in self.atom_dict.items():
            selected = True
            for selection in selec_dict.keys():
                # print("select",selection, selec_dict[selection],".")
                # print("atom",atom[selection],".")
                if atom[selection] not in selec_dict[selection]:
                    selected = False
                    break
            if selected:
                coor_out.atom_dict[atom_num] = atom

        return coor_out

    def get_index_selection(self, selec_dict):
        """Select atom of a coor object based on the change_dict dictionnary.
        Return the list of index of selected atoms.

        :param selec_dict: select ditionnay eg. {"chain" : ["A","G"]}
        :type selec_dict: dict

        :return: list of atom index
        :rtype: list of int

        :Example:

        >>> prot_coor = Coor(os.path.join(TEST_PATH, '1y0m.pdb'))\
        #doctest: +ELLIPSIS
        Succeed to read file ...1y0m.pdb ,  648 atoms found
        >>> prot_coor.get_index_selection({'res_num' : [826,827]})
        [297, 298, 299, 300, 301, 302, 303, 304]

        """

        index_list = []

        for atom_num, atom in self.atom_dict.items():
            selected = True
            # print("atom_num:",atom_num,"atom:",atom)
            for selection in selec_dict.keys():
                # print("select",selection, selec_dict[selection],".")
                # print("selection=",selection)
                # print("atom:",atom)
                # print("atom",atom[selection],".")
                if atom[selection] not in selec_dict[selection]:
                    selected = False
                    break
            if selected:
                # print(atom)
                index_list.append(atom_num)

        return index_list

    def select_from_index(self, index_list):
        """Select atom of a coor object from an atom index
        list.
        Return a new coor object.

        :param index_list: list of index eg. [0, 4]
        :type index_list: list

        :return: a new coor object
        :rtype: coor

        :Example:

        >>> prot_coor = Coor(os.path.join(TEST_PATH, '1y0m.pdb'))\
        #doctest: +ELLIPSIS
        Succeed to read file ...1y0m.pdb ,  648 atoms found
        >>> index_list = prot_coor.get_index_selection({'res_num' : [826,827]})
        >>> index_list
        [297, 298, 299, 300, 301, 302, 303, 304]
        >>> prot_sel = prot_coor.select_from_index(index_list)
        >>> len(prot_sel.atom_dict) == len(index_list)
        True

        """

        coor_out = Coor()

        for atom_num, atom in self.atom_dict.items():
            if atom_num in index_list:
                coor_out.atom_dict[atom_num] = atom

        return coor_out

    def get_attribute_selection(self, selec_dict={}, attribute='uniq_resid',
                                index_list=None):
        """Select atom of a coor object based on the change_dict dictionnary.
        Return the list of unique attribute of the selected atoms.

        :param selec_dict: select ditionnay eg. {"chain" : ["A","G"]}
        :type selec_dict: dict

        :return: list of atom index
        :rtype: list of int

        :Example:

        >>> prot_coor = Coor(os.path.join(TEST_PATH, '1y0m.pdb'))\
        #doctest: +ELLIPSIS
        Succeed to read file ...1y0m.pdb ,  648 atoms found
        >>> prot_coor.get_attribute_selection({'res_num' : [826,827]},\
attribute='uniq_resid')
        [35, 36]

        """

        attr_list = []

        if index_list is None:
            index_list = self.atom_dict.keys()

        # for atom_num, atom in sorted(self.atom_dict.items()):
        for atom_num in index_list:
            atom = self.atom_dict[atom_num]
            selected = True
            # print("atom_num:",atom_num,"atom:",atom)
            for selection in selec_dict.keys():
                # print("select",selection, selec_dict[selection],".")
                # print("selection=",selection)
                # print("atom:",atom)
                # print("atom",atom[selection],".")
                if atom[selection] not in selec_dict[selection]:
                    selected = False
                    break
            if selected:
                # print(atom)
                attr_list.append(atom[attribute])

        return list(set(attr_list))

    def del_atom_index(self, index_list):
        """Delete atoms of a coor object defined by their ``index``.

        :param index_list: list of atom index to delete
        :type index_list: list

        :Example:

        >>> prot_coor = Coor(os.path.join(TEST_PATH, '1y0m.pdb'))\
        #doctest: +ELLIPSIS
        Succeed to read file ...1y0m.pdb ,  648 atoms found
        >>> prot_coor.get_aa_seq()
        {'A': 'TFKSAVKALFDYKAQREDELTFTKSAIIQNVEKQDGGWWRGDYGGKKQLWFPSNYVEEMIN'}
        >>> res_810_852 = prot_coor.get_index_selection({'res_num' :\
range(810,852)})
        >>> prot_coor.del_atom_index(index_list = res_810_852)\
        #doctest: +ELLIPSIS
        <...Coor object at ...>
        >>> prot_coor.get_aa_seq()
        {'A': 'TFKSAVKALFDYKAQREDE'}

        """

        for index in index_list:
            # print(index, self.atom_dict[index])
            del self.atom_dict[index]

        return self

    def correct_chain(self, Ca_cutoff=4.5):
        """Correct the chain ID's of a coor object, by checking consecutive
        Calphas atoms distance. If the distance is higher than ``Ca_cutoff``
        , the former atoms are considered as in a different chain.

        :param Ca_cutoff: cutoff for distances between Calphas atoms (X)
        :type Ca_cutoff: float, default=4.5

        :Example:

        >>> prot_coor = Coor(os.path.join(TEST_PATH, '1y0m.pdb'))\
        #doctest: +ELLIPSIS
        Succeed to read file ...1y0m.pdb ,  648 atoms found
        >>> res_810 = prot_coor.get_index_selection({'res_num' : [810]})
        >>> prot_coor = prot_coor.del_atom_index(index_list = res_810)
        >>> prot_coor.get_aa_seq()
        {'A': 'TFKSAVKALFDYKAQREDETFTKSAIIQNVEKQDGGWWRGDYGGKKQLWFPSNYVEEMIN'}
        >>> prot_coor.correct_chain() #doctest: +ELLIPSIS
        Chain: A  Residue: 0 to 18
        Chain: B  Residue: 20 to 60
        <...Coor object at ...>
        >>> # As a residue is missing, Calphas after residue 18 is no
        >>> # more consecutive


        .. note::
            This is specially usefull for pdb2gmx which cut the protein
            chains based on the chain ID's.

        """

        Ca_atom = self.select_part_dict({"name": ["CA"]})
        first_flag = True
        chain_res_list = []
        res_list = []
        old_atom = None

        # Identify Chain uniq_resid
        # Need to use sorted to be sure to check consecutive residues (atoms)
        for key, atom in sorted(Ca_atom.atom_dict.items()):
            if first_flag:
                first_flag = False
            else:
                distance = Coor.atom_dist(atom, old_atom)
                # print(distance)
                if distance < Ca_cutoff:
                    old_atom = atom
                else:
                    # print("New chain")
                    chain_res_list.append(res_list)
                    res_list = []
            res_list.append(atom['uniq_resid'])
            old_atom = atom
        chain_res_list.append(res_list)

        # Change chain ID :
        for i, chain_res in enumerate(chain_res_list):
            logger.info("Chain: {}  Residue: {} to {}".format(
                chr(65 + i), chain_res[0], chain_res[-1]))
            chain_index = self.get_index_selection({'uniq_resid': chain_res})
            # print(chain_index)
            self.change_index_pdb_field(chain_index, {"chain": chr(65 + i)})
            # pdb_dict_out.update(chain_dict)

        # print(pdb_dict_out)

        return self

    def correct_his_name(self):
        """ Get his protonation state from pdb2pqr and replace HIS resname
        with HSE, HSD, HSP resname.
        To do after pdb2pqr, in order that protonation is recognize by pdb2gmx.

        :Example:

        >>> try: #doctest: +ELLIPSIS
        ...   print("Start import")
        ...   from . import pdb2pqr
        ... except ImportError:
        ...   import pdb2pqr
        Start import...
        >>> TEST_OUT = str(getfixture('tmpdir'))
        >>> # Compute protonation with pdb2pqr:
        >>> pdb2pqr.compute_pdb2pqr(os.path.join(TEST_PATH, '4n1m.pdb'),\
os.path.join(TEST_OUT, '4n1m.pqr')) #doctest: +ELLIPSIS
        Succeed to read file ...4n1m.pdb ,  2530 atoms found
        Succeed to save file ...tmp_pdb2pqr.pdb
        pdb2pqr... --ff CHARMM --ffout CHARMM --chain \
--ph-calc-method=propka ...tmp_pdb2pqr.pdb ...4n1m.pqr
        0
        >>> prot_coor = Coor(os.path.join(TEST_OUT, '4n1m.pqr')) \
#doctest: +ELLIPSIS
        Succeed to read file ...4n1m.pqr ,  2549 atoms found
        >>> HSD_index = prot_coor.get_index_selection({'res_name' : ['HSD'],\
'name':['CA']})
        >>> print(len(HSD_index))
        4
        >>> HSE_index = prot_coor.get_index_selection({'res_name' : ['HSE'],\
'name':['CA']})
        >>> print(len(HSE_index))
        0
        >>> HSP_index = prot_coor.get_index_selection({'res_name' : ['HSP'],\
'name':['CA']})
        >>> print(len(HSP_index))
        1
        >>> prot_coor.correct_his_name() #doctest: +ELLIPSIS
        <...Coor object at 0x...
        >>> HIS_index = prot_coor.get_index_selection({'res_name' : ['HIS'],\
'name':['CA']})
        >>> print(len(HIS_index))
        0

        .. note::
            This function seems useless. Since last version of pdb2pqr
            residue name seems correct.
        """

        # FIND HISTIDINE res

        # HSD:
        hsd_uniq_res = self.get_attribute_selection({"res_name": ["HIS"],
                                                     "name": ["HD1"]},
                                                    attribute='uniq_resid')
        # HSE:
        hse_uniq_res = self.get_attribute_selection({"res_name": ["HIS"],
                                                     "name": ["HE2"]},
                                                    attribute='uniq_resid')
        # HSP: find res in common with both hsd and hse
        hsp_uniq_res = [res for res in hsd_uniq_res if res in hse_uniq_res]
        # remove HSP res from HSE HSD list
        if hsp_uniq_res:
            for res in hsp_uniq_res:
                hsd_uniq_res.remove(res)
                hse_uniq_res.remove(res)

        # Replace HIS resname :
        all_his_uniq_res = hsd_uniq_res + hse_uniq_res + hsp_uniq_res

        for atom_num, atom in self.atom_dict.items():
            if atom["uniq_resid"] in all_his_uniq_res:
                if atom["uniq_resid"] in hsd_uniq_res:
                    atom["res_name"] = "HSD"
                elif atom["uniq_resid"] in hse_uniq_res:
                    atom["res_name"] = "HSE"
                else:
                    atom["res_name"] = "HSP"

        return self

    def add_zinc_finger(self, ZN_pdb, cutoff=3.2):
        """ Change protonation state of cysteins and histidine coordinating
        Zinc atoms.
        To do after `correct_his_name` and `correct_cys_name`,
        in order that protonation is recognize by pdb2gmx.

        :Example:

        >>> try: #doctest: +ELLIPSIS
        ...   print("Start import")
        ...   from . import pdb2pqr
        ... except ImportError:
        ...   import pdb2pqr
        Start import...
        >>> TEST_OUT = str(getfixture('tmpdir'))
        >>> show_log()
        >>> # Read the pdb 1jd4 and keep only chain A
        >>> input_pdb = Coor(os.path.join(TEST_PATH, '1jd4.pdb'))\
        #doctest: +ELLIPSIS
        Succeed to read file ...1jd4.pdb ,  1586 atoms found
        >>> chain_A = input_pdb.select_part_dict(selec_dict =\
{'chain' : ['A']})
        >>> chain_A.write_pdb(os.path.join(TEST_OUT, '1jd4_A.pdb'))\
        #doctest: +ELLIPSIS
        Succeed to save file ...1jd4_A.pdb
        >>>
        >>> # Compute protonation with pdb2pqr:
        >>> pdb2pqr.compute_pdb2pqr(os.path.join(TEST_OUT, '1jd4_A.pdb'),\
os.path.join(TEST_OUT, '1jd4.pqr')) #doctest: +ELLIPSIS
        Succeed to read file ...1jd4_A.pdb ,  793 atoms found
        Succeed to save file ...tmp_pdb2pqr.pdb
        pdb2pqr... --ff CHARMM --ffout CHARMM --chain --ph-calc-method=propka \
...tmp_pdb2pqr.pdb ...1jd4.pqr
        0
        >>> prot_coor = Coor(os.path.join(TEST_OUT, '1jd4.pqr'))
        Succeed to read file ...1jd4.pqr ,  1549 atoms found
        >>> prot_coor.correct_cys_name() #doctest: +ELLIPSIS
        <...Coor object at 0x...
        >>> prot_coor.correct_his_name() #doctest: +ELLIPSIS
        <...Coor object at 0x...
        >>> prot_coor.correct_chain() #doctest: +ELLIPSIS
        Chain: A  Residue: 0 to 95
        <...Coor object at 0x...
        >>> ZN_index = prot_coor.get_index_selection({'name':['ZN']})
        >>> print(len(ZN_index))
        0
        >>> prot_coor.add_zinc_finger(os.path.join(TEST_OUT, '1jd4_A.pdb'))\
        #doctest: +ELLIPSIS
        Succeed to read file ...1jd4_A.pdb ,  793 atoms found
        Presence of 1 Zinc detected
        change cystein residue(s) : [48, 51, 75]
        change histidine residue(s) : [68]
        True
        >>> ZN_index = prot_coor.get_index_selection({'name':['ZN']})
        >>> print(len(ZN_index))
        1

        .. note::
            This function seems useless. Since last version of pdb2pqr
            residue name seems correct.
        """

        # Check the number of ZN atoms:
        coor_pre_pqr = Coor()
        coor_pre_pqr.read_pdb(ZN_pdb)
        Zinc_sel = coor_pre_pqr.select_part_dict(selec_dict={'name': ['ZN']})
        Zinc_num = Zinc_sel.num

        if Zinc_num == 0:
            return False
        else:
            logger.info("Presence of {} Zinc detected".format(Zinc_num))

        # Add the Zinc atoms:
        for key, val in Zinc_sel.atom_dict.items():
            atom = val
            atom['chain'] = 'Z'
            self.atom_dict[self.num] = val

        # Check cystein and histidine atoms close to ZN:
        close_atom = self.dist_under_index(Zinc_sel, cutoff=cutoff)
        cys_uniq_res_list = []
        his_uniq_res_list = []
        for atom in close_atom:
            local_atom = self.atom_dict[atom]
            if local_atom['res_name'] in ['CYS']:
                cys_uniq_res_list.append(local_atom['uniq_resid'])
            if local_atom['res_name'] in ['HIS', 'HSD', 'HSE', 'HSP']:
                his_uniq_res_list.append(local_atom['uniq_resid'])

        cys_uniq_res_list = list(set(cys_uniq_res_list))
        his_uniq_res_list = list(set(his_uniq_res_list))

        # Change CYS to CYN:
        logger.info("change cystein residue(s) : {}".format(
            sorted(cys_uniq_res_list)))
        to_change_cys = self.get_index_selection(
            {'uniq_resid': cys_uniq_res_list})
        self.change_index_pdb_field(to_change_cys, {'res_name': 'CYN'})

        # Change Histidine to HSD or HSE
        # the non protonated nitrogen have to be the closest to the ZN atom:
        logger.info("change histidine residue(s) : {}".format(
            sorted(his_uniq_res_list)))
        for his_uniq_res in his_uniq_res_list:
            # NE2 ND1
            epsilon_his_index = self.get_index_selection(
                {'uniq_resid': [his_uniq_res], 'name': ['NE2']})[0]
            delta_his_index = self.get_index_selection(
                {'uniq_resid': [his_uniq_res], 'name': ['ND1']})[0]
            for key, val in Zinc_sel.atom_dict.items():

                epsilon_dist = Coor.atom_dist(
                    val, self.atom_dict[epsilon_his_index])
                delta_dist = Coor.atom_dist(
                    val, self.atom_dict[delta_his_index])

                if (epsilon_dist < cutoff or delta_dist < cutoff):
                    to_change_his = self.get_index_selection(
                        {'uniq_resid': [his_uniq_res]})
                    if (epsilon_dist < delta_dist):
                        self.change_index_pdb_field(
                            to_change_his, {'res_name': 'HSD'})
                    else:
                        self.change_index_pdb_field(
                            to_change_his, {'res_name': 'HSE'})

        return True

    def correct_cys_name(self):
        """ Correct the CYS resname from pdb2pqr

        :Example:

        >>> try: #doctest: +ELLIPSIS
        ...   print("Start import")
        ...   from . import pdb2pqr
        ... except ImportError:
        ...   import pdb2pqr
        Start import...
        >>> TEST_OUT = str(getfixture('tmpdir'))
        >>> # Compute protonation with pdb2pqr:
        >>> pdb2pqr.compute_pdb2pqr(os.path.join(TEST_PATH, '1dpx.pdb'),\
os.path.join(TEST_OUT, '1dpx.pqr')) #doctest: +ELLIPSIS
        Succeed to read file ...1dpx.pdb ,  1192 atoms found
        Succeed to save file ...tmp_pdb2pqr.pdb
        pdb2pqr... --ff CHARMM --ffout CHARMM --chain --ph-calc-method\
=propka ...tmp_pdb2pqr.pdb ...1dpx.pqr
        0
        >>> prot_coor = Coor(os.path.join(TEST_OUT, '1dpx.pqr'))
        Succeed to read file ...1dpx.pqr ,  1961 atoms found
        >>> Isu_index = prot_coor.get_index_selection({'res_name' : ['DISU']})
        >>> print(len(Isu_index))
        16
        >>> prot_coor.correct_cys_name() #doctest: +ELLIPSIS
        <...Coor object at 0x...
        >>> Isu_index = prot_coor.get_index_selection({'res_name' : ['DISU']})
        >>> print(len(Isu_index))
        0
        """

        # FIND ISU res
        isu_index_list = self.get_index_selection({"res_name": ["DISU"]})
        if not isu_index_list:
            # print("Nothing to Fix")
            return self

        # Replace CYS resname :

        for atom_num in isu_index_list:
            self.atom_dict[atom_num]["res_name"] = "CYS"
            if self.atom_dict[atom_num]["name"] == "1CB":
                self.atom_dict[atom_num]["name"] = "CB"
            if self.atom_dict[atom_num]["name"] == "1SG":
                self.atom_dict[atom_num]["name"] = "SG"

        return self

    def correct_ion_octa(self, ion_name):
        """ For specified ion, create an octahedral dummy model described
        by a set of 6 cationic dummy atoms connected around a central metal
        atom.
        From Duarte et al. J Phys Chem B 2014.

        :param ion_name: name of metal present in .pdb to transform
            in octahedral dummy model
        :type ion_name: str

        :Example:

        >>> TEST_OUT = str(getfixture('tmpdir'))
        >>> prot_coor = Coor(os.path.join(TEST_PATH, '1jd4.pdb'))\
        #doctest: +ELLIPSIS
        Succeed to read file ...1jd4.pdb ,  1586 atoms found
        >>> ion_index = prot_coor.get_index_selection({'res_name' : ['ZN']})
        >>> print(len(ion_index))
        2
        >>> prot_coor.correct_ion_octa('ZN') #doctest: +ELLIPSIS
        <...Coor object at 0x...
        >>> ion_index = prot_coor.get_index_selection({'res_name' : ['ZN']})
        >>> print(len(ion_index))
        14

        """

        index = 0
        new_atom_dict = dict()

        for atom_num, atom in sorted(self.atom_dict.items()):
            if ion_name in atom["name"]:

                ion_atom = atom.copy()
                ion_atom["num"] = index
                new_atom_dict[index] = ion_atom
                index += 1

                ion_atom = atom.copy()
                ion_atom["name"] = "D1"
                ion_atom["num"] = index
                ion_atom["xyz"] = np.array([ion_atom["xyz"][0] + float(0.9),
                                            ion_atom["xyz"][1],
                                            ion_atom["xyz"][2]])
                new_atom_dict[index] = ion_atom
                index += 1

                ion_atom = atom.copy()
                ion_atom["name"] = "D2"
                ion_atom["num"] = index
                ion_atom["xyz"] = np.array([ion_atom["xyz"][0] - float(0.9),
                                            ion_atom["xyz"][1],
                                            ion_atom["xyz"][2]])
                new_atom_dict[index] = ion_atom
                index += 1

                ion_atom = atom.copy()
                ion_atom["name"] = "D3"
                ion_atom["num"] = index
                ion_atom["xyz"] = np.array([ion_atom["xyz"][0],
                                            ion_atom["xyz"][1] + float(0.9),
                                            ion_atom["xyz"][2]])
                new_atom_dict[index] = ion_atom
                index += 1

                ion_atom = atom.copy()
                ion_atom["name"] = "D4"
                ion_atom["num"] = index
                ion_atom["xyz"] = np.array([ion_atom["xyz"][0],
                                            ion_atom["xyz"][1] - float(0.9),
                                            ion_atom["xyz"][2]])
                new_atom_dict[index] = ion_atom
                index += 1

                ion_atom = atom.copy()
                ion_atom["name"] = "D5"
                ion_atom["num"] = index
                ion_atom["xyz"] = np.array([ion_atom["xyz"][0],
                                            ion_atom["xyz"][1],
                                            ion_atom["xyz"][2] + float(0.9)])
                new_atom_dict[index] = ion_atom
                index += 1

                ion_atom = atom.copy()
                ion_atom["name"] = "D6"
                ion_atom["num"] = index
                ion_atom["xyz"] = np.array([ion_atom["xyz"][0],
                                            ion_atom["xyz"][1],
                                            ion_atom["xyz"][2] - float(0.9)])
                new_atom_dict[index] = ion_atom
                index += 1

            else:
                atom["num"] = index
                new_atom_dict[index] = atom
                index += 1

        self.atom_dict = new_atom_dict

        return self

    def water_to_ATOM(self):
        """ Change `HETATM` field of water to `ATOM`, as pdb2pqr only use ATOM field.

        :Example:

        >>> try: #doctest: +ELLIPSIS
        ...   print("Start import")
        ...   from . import pdb2pqr
        ... except ImportError:
        ...   import pdb2pqr
        Start import...
        >>> prot_coor = Coor(os.path.join(TEST_PATH, '1dpx.pdb'))\
        #doctest: +ELLIPSIS
        Succeed to read file ...1dpx.pdb ,  1192 atoms found
        >>> hetatm_index = prot_coor.get_index_selection({'field':['HETATM']})
        >>> print(len(hetatm_index))
        179
        >>> prot_coor.water_to_ATOM() #doctest: +ELLIPSIS
        <...Coor object at 0x...
        >>> hetatm_index = prot_coor.get_index_selection({'field':['HETATM']})
        >>> print(len(hetatm_index))
        2
        >>> water_index = prot_coor.get_index_selection({'res_name':['HOH']})
        >>> print(len(water_index))
        177
        """

        # FIND Water res
        water_index_list = self.get_index_selection(
            {'res_name': ['HOH'], 'field': ['HETATM']})
        if not water_index_list:
            return self
        else:
            self.change_index_pdb_field(water_index_list, {'field': 'ATOM'})
            return self

    def correct_water_name(self):
        """ Correct the water resname from pdb2pqr

        :Example:

        >>> try: #doctest: +ELLIPSIS
        ...   print("Start import")
        ...   from . import pdb2pqr
        ... except ImportError:
        ...   import pdb2pqr
        Start import...
        >>> TEST_OUT = str(getfixture('tmpdir'))
        >>> prot_coor = Coor(os.path.join(TEST_PATH, '1dpx.pdb'))\
        #doctest: +ELLIPSIS
        Succeed to read file ...1dpx.pdb ,  1192 atoms found
        >>> prot_coor.water_to_ATOM() #doctest: +ELLIPSIS
        <...Coor object at 0x...
        >>> prot_coor.write_pdb(os.path.join(TEST_OUT, '1dpx_water.pdb'))\
        #doctest: +ELLIPSIS
        Succeed to save file ...1dpx_water.pdb
        >>> # Compute protonation with pdb2pqr:
        >>> pdb2pqr.compute_pdb2pqr(
        ... os.path.join(TEST_OUT, '1dpx_water.pdb'),
        ... os.path.join(TEST_OUT, '1dpx_water.pqr')) #doctest: +ELLIPSIS
        Succeed to read file ...1dpx_water.pdb ,  1192 atoms found
        Succeed to save file ...tmp_pdb2pqr.pdb
        pdb2pqr... --ff CHARMM --ffout CHARMM --chain --ph-calc-method=propka \
...tmp_pdb2pqr.pdb ...1dpx_water.pqr
        0
        >>> prot_coor = Coor(os.path.join(
        ... TEST_OUT, '1dpx_water.pqr')) #doctest: +ELLIPSIS
        Succeed to read file ...1dpx_water.pqr ,  2492 atoms found
        >>> water_index = prot_coor.get_index_selection(
        ... {'res_name':['TP3M'], 'name':['OH2']})
        >>> print(len(water_index))
        177
        """

        # FIND Water res
        water_index_list = self.get_index_selection({"res_name": ["TP3M"]})
        if not water_index_list:
            # print("Nothing no water to fix")
            return self

        # Replace SOL resname :

        for atom_num in water_index_list:
            self.atom_dict[atom_num]["res_name"] = "SOL"
            if self.atom_dict[atom_num]["name"] == "OH2":
                self.atom_dict[atom_num]["name"] = "OW"
            if self.atom_dict[atom_num]["name"] == "H1":
                self.atom_dict[atom_num]["name"] = "HW1"
            if self.atom_dict[atom_num]["name"] == "H2":
                self.atom_dict[atom_num]["name"] = "HW2"

        return self

    def insert_mol(self, pdb_out, out_folder, mol_chain, mol_num,
                   check_file_out=True, prot_atom_name=['CA'],
                   sol_res_name=['SOL'], sol_atom_name=['OW']):
        """
        Insert molecules defined by chain ID ``mol_chain`` in a water solvant.
        Check which water molecules are within ``cutoff_prot_off=12.0`` X
        and ``cutoff_prot_in=15.0`` X of protein and peptide C alpha atoms.
        Move the molecules to be inserted at the position of water molecules.
        Then delete all water molecules within ``cutoff_water_clash=1.2`` X of
        the inserted molecule atoms.

        :param pdb_out: name of output pdb file
        :type pdb_out: str

        :param out_folder: path of the ouput directory
        :type out_folder: str

        :param mol_chain: chain ID of the molecule to be inserted,
        :type mol_chain: str

        :param mol_num: Number of molecule to be inserted,
        :type mol_num: int

        :param check_file_out: flag to check or not if
            file has already been created.
            If the file is present then the command break.
        :type check_file_out: bool, optional, default=True

        .. warning::
            self.atom_dict file must contain alredy a concatenated
            system with a ligand (chain: ``mol_chain``) and a
            solvated system.
            Molecules to insert should have different residue number
            or at least non consecutive.
        """

        # Create the out_folder:
        pdb_out = os.path.join(out_folder, pdb_out)
        os_command.create_dir(out_folder)

        # Parameters for molecule insertion:
        cutoff_water_clash = 1.2
        cutoff_prot_off = 12.0
        cutoff_prot_in = 15.0
        cutoff_mol_off = 5.0

        logger.info("Insert mol in system")

        # Check if output files exist:
        if check_file_out and os.path.isfile(pdb_out):
            logger.info("Insert Mol %s already exist" % pdb_out)
            return None

        # Select protein, water and molecule atoms :
        prot_insert_CA = self.select_part_dict(selec_dict={'name': prot_atom_name})
        water = self.select_part_dict(selec_dict={'res_name': sol_res_name})
        water_O = self.select_part_dict(
            selec_dict={'res_name': sol_res_name, 'name': sol_atom_name})
        insert = self.select_part_dict(selec_dict={'chain': [mol_chain]})
        #insert_index_list = self.get_index_selection(
        #  selec_dict={'chain': [mol_chain]})
        #lig_res_name = insert.atom_dict[insert_index_list[0]]['res_name']
        #lig_atom_name = insert.atom_dict[insert_index_list[0]]['name']
        # print(insert.atom_dict[insert_index_list[0]])
        #lig_res_num = insert.atom_dict[insert_index_list[0]]['res_num']

        #insert_first_atom = self.select_part_dict(
        #    selec_dict={'chain': [mol_chain], 'name': [lig_atom_name],
        #                'res_name': [lig_res_name],
        #                'res_num': [lig_res_num]})

        #mol_num = insert_first_atom.num
        #print(insert.atom_dict)
        res_insert_list = list(set(insert.get_attribute_selection(
            attribute='uniq_resid')))
        # Need to sort the resid, to have consecutive residues
        res_insert_list.sort()
        # print(res_insert_list)
        logger.info("Residue list = {}".format(res_insert_list))
        mol_len = int(len(res_insert_list) / mol_num)

        logger.info("Insert {} mol of {:d} residues each".format(
            mol_num, mol_len))
        start_time = time.time()
        # Insert one molecule at a time:
        for i in range(mol_num):
            start_time = time.time()
            # Get water not too close and too far from the protein
            water_prot_index = water_O.get_index_dist_between(
                prot_insert_CA, cutoff_max=cutoff_prot_in,
                cutoff_min=cutoff_prot_off)
            # Get water not close from ligand:
            water_ligand_index = water_O.get_index_dist_between(
                insert, cutoff_max=1e10,
                cutoff_min=cutoff_mol_off)
            water_good_index = [value for value in water_prot_index
                                if value in water_ligand_index]
            logger.info('insert mol {:3}, water mol {:5}, time={:.2f}'.format(
                i + 1, len(water_good_index), time.time() - start_time))
            insert_unique = insert.select_part_dict(
                selec_dict={'chain': [mol_chain],
                            'uniq_resid': res_insert_list[
                            (mol_len * i):(mol_len * (i + 1))]})
            com_insert = insert_unique.center_of_mass()
            # print(self.atom_dict[water_good_index[0]]['uniq_resid'])
            trans_vector = self.atom_dict[water_good_index[0]]['xyz'] -\
                com_insert

            insert_unique.translate(trans_vector)

        # Delete water residues in which at leat one atom is close
        # enough to peptides
        water_to_del_index = water.get_index_dist_between(
            insert, cutoff_max=cutoff_water_clash)
        water_res_to_del = water.get_attribute_selection(
            selec_dict={}, attribute='uniq_resid',
            index_list=water_to_del_index)
        water_to_del_index = water.get_index_selection(
            selec_dict={'uniq_resid': water_res_to_del})
        logger.info("Delete {} overlapping water atoms".format(
            len(water_to_del_index)))
        self.del_atom_index(index_list=water_to_del_index)

        self.write_pdb(pdb_out)
        return self

    def translate(self, vector):
        """ Translate all atoms of a coor object by a given ``vector``

        :param vector: 3d translation vector
        :type vector: list

        :Example:

        >>> prot_coor = Coor(os.path.join(TEST_PATH, '1y0m.pdb'))\
        #doctest: +ELLIPSIS
        Succeed to read file ...1y0m.pdb ,  648 atoms found
        >>> com_1y0m = prot_coor.center_of_mass()
        >>> print("x:{:.2f} y:{:.2f} z:{:.2f}".format(*com_1y0m))
        x:16.01 y:0.45 z:8.57
        >>> prot_coor.translate(-com_1y0m)
        >>> com_1y0m = prot_coor.center_of_mass()
        >>> print("x:{:.2f} y:{:.2f} z:{:.2f}".format(*com_1y0m))
        x:-0.00 y:0.00 z:0.00

        """
        for atom_num, atom in self.atom_dict.items():
            atom['xyz'] += vector

        return

    def get_mass_array(self):
        """ Extract mass of each `atom_dict` and return it as an numpy array
        Avoid using atoms with 2 letters atom name like NA Cl ...

        :return: mass array
        :rtype: np.array


        :Example:

        >>> prot_coor = Coor(os.path.join(TEST_PATH, '1y0m.pdb'))\
        #doctest: +ELLIPSIS
        Succeed to read file ...1y0m.pdb ,  648 atoms found
        >>> mass_1y0m = prot_coor.get_mass_array()
        >>> print("Mass of 10 first atoms: {}".format(mass_1y0m[:10]))
        Mass of 10 first atoms: [7 6 6 8 6 8 6 7 6 6]
        >>> prot_coor_ca = prot_coor.select_part_dict({'name':['CA']})
        >>> mass_1y0m_ca = prot_coor_ca.get_mass_array()
        >>> print("Mass of 10 first atoms: {}".format(mass_1y0m_ca[:10]))
        Mass of 10 first atoms: [6 6 6 6 6 6 6 6 6 6]
        """

        name_array = self.get_array(field='name')
        mass_list = []
        for name in name_array:
            if name[0] in ATOM_MASS_DIST:
                mass_list.append(ATOM_MASS_DIST[name[0]])
            else:
                logger.warning('Warning atom {} mass could not be'
                               ' founded'.format(name))
                mass_list.append(0)
        mass_array = np.array(mass_list)

        return mass_array

    def center_of_mass(self, selec_dict={}):
        """ Compute the center of mass of a selection
        Avoid using atoms with 2 letters atom name like NA Cl ...
        If selection is empy, take all atoms.

        :param selec_dict: selection dictionnary
        :type selec_dict: dict, default={}

        :Example:

        >>> prot_coor = Coor(os.path.join(TEST_PATH, '1y0m.pdb'))\
        #doctest: +ELLIPSIS
        Succeed to read file ...1y0m.pdb ,  648 atoms found
        >>> com_1y0m = prot_coor.center_of_mass()
        >>> print("x:{:.2f} y:{:.2f} z:{:.2f}".format(*com_1y0m))
        x:16.01 y:0.45 z:8.57
        >>> com_1y0m_ca = prot_coor.center_of_mass({'name':['CA']})
        >>> print("x:{:.2f} y:{:.2f} z:{:.2f}".format(*com_1y0m_ca))
        x:15.95 y:0.72 z:8.96

        .. warning::
            Atom name must start with its type letter (H, C, N, O, P, S).
        """

        local_select = self.select_part_dict(selec_dict=selec_dict)
        coor_array = local_select.get_array()
        mass_array = local_select.get_mass_array()

        return (coor_array.T * mass_array).sum(axis=1) / sum(mass_array)

    def centroid(self, selec_dict={}):
        """ Compute the centroid of a selection
        If selection is empy, take all atoms.

        :param selec_dict: selection dictionnary
        :type selec_dict: dict, default={}

        :Example:

        >>> prot_coor = Coor(os.path.join(TEST_PATH, '1y0m.pdb'))\
        #doctest: +ELLIPSIS
        Succeed to read file ...1y0m.pdb ,  648 atoms found
        >>> com_1y0m = prot_coor.centroid()
        >>> print("x:{:.2f} y:{:.2f} z:{:.2f}".format(*com_1y0m))
        x:16.03 y:0.44 z:8.57
        >>> com_1y0m_ca = prot_coor.centroid(selec_dict={'name':['CA']})
        >>> print("x:{:.2f} y:{:.2f} z:{:.2f}".format(*com_1y0m_ca))
        x:15.95 y:0.72 z:8.96


        """

        coor_array = self.select_part_dict(selec_dict=selec_dict).get_array()

        return coor_array.mean(axis=0)

    def get_box_dim(self, selec_dict={}):
        """ Compute the x, y, z dimension of a selection

        :param selec_dict: selection dictionnary
        :type selec_dict: dict, default={}

        :Example:

        >>> prot_coor = Coor(os.path.join(TEST_PATH, '1y0m.pdb'))\
        #doctest: +ELLIPSIS
        Succeed to read file ...1y0m.pdb ,  648 atoms found
        >>> box_1y0m = prot_coor.get_box_dim()
        >>> print("x:{:.2f} y:{:.2f} z:{:.2f}".format(*box_1y0m))
        x:39.14 y:29.36 z:31.19

        .. warning::
            Atom name must start with its type letter (H, C, N, O, P, S).
        """

        if len(selec_dict) > 0:
            sel_coor = self.select_part_dict(selec_dict=selec_dict)
            coor_array = sel_coor.get_array()
        else:
            coor_array = self.get_array()

        min_val = np.amin(coor_array, axis=0)
        max_val = np.amax(coor_array, axis=0)

        return max_val - min_val

    def get_index_dist_between(self, atom_sel_2, cutoff_min=0, cutoff_max=10):
        """ Check is distance between atoms of self.atom_dict is under cutoff
        with the atoms of group 1.
        Then return list of index of atoms of self.coor under cutoff ditance.

        :param atom_sel_1: atom dictionnary
        :type atom_sel_1: dict

        :param atom_sel_2: atom dictionnary
        :type atom_sel_2: dict

        :param cutoff_min: maximum distance cutoff
        :type cutoff_min: float, default=0.0

        :param cutoff_max: minimum distance cutoff
        :type cutoff_max: float, default=10.0

        :Example:

        >>> prot_coor = Coor(os.path.join(TEST_PATH, '1y0m.pdb'))\
        #doctest: +ELLIPSIS
        Succeed to read file ...1y0m.pdb ,  648 atoms found
        >>> res_810 = prot_coor.select_part_dict({'res_num' : [810]})
        >>> close_r810 = prot_coor.get_index_dist_between(
        ... res_810, cutoff_min = 3, cutoff_max = 5)
        >>> print(len(close_r810))
        65

        """

        coor_array = self.get_array()
        index_array = np.array(
            [key for key, atom in self.atom_dict.items()])

        coor_array_2 = atom_sel_2.get_array()

        # Compute distance matrix
        dist_mat = distance_matrix(coor_array, coor_array_2)

        # Compute index of matrix column under cutoff_max and over cutoff_min:
        dist_mat_good = np.where(
            (dist_mat.min(1) < cutoff_max) & (dist_mat.min(1) > cutoff_min))[0]

        return index_array[dist_mat_good]

    def compute_rmsd_to(self, atom_sel_2, selec_dict={'name': ['CA']},
                        index_list=None):
        """ Compute RMSD between two atom_dict
        Then return the RMSD value.

        :param atom_sel_1: atom dictionnary
        :type atom_sel_1: dict

        :param atom_sel_2: atom dictionnary
        :type atom_sel_2: dict

        :param selec_dict: selection dictionnary
        :type selec_dict: dict, default={}

        :return: RMSD
        :rtype: float

        :Example:

        >>> prot_coor = Coor()

        """

        if index_list is None:
            sel_1_coor = self.select_part_dict(selec_dict=selec_dict)
            sel_2_coor = atom_sel_2.select_part_dict(selec_dict=selec_dict)
            coor_array_1 = sel_1_coor.get_array()
            coor_array_2 = sel_2_coor.get_array()
        else:
            coor_array_1 = self.get_array(index_list=index_list[0])
            coor_array_2 = atom_sel_2.get_array(index_list=index_list[1])

        diff = coor_array_1 - coor_array_2
        N = len(coor_array_1)

        rmsd = np.sqrt((diff * diff).sum() / N)

        return rmsd

    def remove_alter_position(self):
        """ Remove alternative position.

        :Example:

        >>> prot_coor = Coor(os.path.join(TEST_PATH, '4n1m.pdb'))\
        #doctest: +ELLIPSIS
        Succeed to read file ...4n1m.pdb ,  2530 atoms found
        >>> print('Atom num = {}'.format(prot_coor.num))
        Atom num = 2530
        >>> prot_coor.remove_alter_position()
        >>> print('Atom num = {}'.format(prot_coor.num))
        Atom num = 2475
        """

        # Remove alter_loc B, C, D
        alter_loc_bcd = self.get_index_selection(
            {'alter_loc': ['B', 'C', 'D']})
        self.del_atom_index(index_list=alter_loc_bcd)
        self.change_pdb_field(change_dict={"alter_loc": ""})
        return

    def align_to(self, atom_sel_2, selec_dict={'name': ['CA']},
                 index_list=None, rot_kabsch=True):
        """ Align structure to an Coor object.

        :param atom_sel_1: atom dictionnary
        :type atom_sel_1: dict

        :param atom_sel_2: atom dictionnary
        :type atom_sel_2: dict

        :param selec_dict: selection dictionnary
        :type selec_dict: dict, default={'name': ['CA']}

        :param rot_kabsch: method for rotation kabsh, if not quaternion
        :type rot_kabsch: bool, default=True

        :Example:

        >>> prot_coor = Coor(os.path.join(TEST_PATH, '1jd4.pdb'))\
        #doctest: +ELLIPSIS
        Succeed to read file ...1jd4.pdb ,  1586 atoms found
        >>> chain_A = prot_coor.select_part_dict(selec_dict={'chain': ['A']})
        >>> chain_B = prot_coor.select_part_dict(selec_dict={'chain': ['B']})
        >>> rmsd = chain_A.compute_rmsd_to(chain_B)
        >>> print('RMSD before alignement is {:.2f} Å'.format(rmsd))
        RMSD before alignement is 37.47 Å
        >>> chain_A.align_to(chain_B)
        >>> rmsd = chain_A.compute_rmsd_to(chain_B)
        >>> print('RMSD after alignement is {:.2f} Å'.format(rmsd))
        RMSD after alignement is 0.06 Å
        >>> chain_A.align_to(chain_B)
        >>> rmsd = chain_A.compute_rmsd_to(chain_B)
        >>> print('RMSD after 2nd alignement is {:.2f} Å'.format(rmsd))
        RMSD after 2nd alignement is 0.06 Å

        """

        if index_list is None:
            sel_1_coor = self.select_part_dict(selec_dict=selec_dict)
            sel_2_coor = atom_sel_2.select_part_dict(selec_dict=selec_dict)

            coor_array_1 = sel_1_coor.get_array()
            coor_array_2 = sel_2_coor.get_array()
        else:
            # print(index_list[0], index_list[1])
            coor_array_1 = self.get_array(index_list=index_list[0])
            coor_array_2 = atom_sel_2.get_array(index_list=index_list[1])

        all_coor_array_1 = self.get_array()

        centroid_1 = coor_array_1.mean(axis=0)
        centroid_2 = coor_array_2.mean(axis=0)
        coor_array_1 -= centroid_1
        coor_array_2 -= centroid_2

        if rot_kabsch:
            rot_mat = Coor.kabsch(coor_array_1, coor_array_2)
        else:
            rot_mat = Coor.quaternion_rotate(coor_array_1, coor_array_2)

        all_coor_array_1 -= centroid_1
        all_coor_array_1 = np.dot(all_coor_array_1, rot_mat)
        all_coor_array_1 += centroid_2

        for i, atom_num in enumerate(sorted(self.atom_dict)):
            self.atom_dict[atom_num]['xyz'] = all_coor_array_1[i]

        return

    def align_seq_coor_to(self, atom_sel_2, chain_1=['A'], chain_2=['A']):
        """ Align 2 strucures, using a sequence alignement to determine
        which residue to align.
        Compute RMSD between two atom_dict
        Then return the RMSD value.


        :param atom_sel_2: atom dictionnary
        :type atom_sel_2: dict

        :param chain_1: list of chain
        :type chain_1: list

        :param chain_2: list of chain
        :type chain_2: list

        :return: rmsd and alignement index
        :rtype: float and list

        :Example:

        >>> prot_1_coor = Coor(os.path.join(TEST_PATH, '1jd4.pdb'))\
        #doctest: +ELLIPSIS
        Succeed to read file ...1jd4.pdb ,  1586 atoms found
        >>> prot_2_coor = Coor(os.path.join(TEST_PATH, '1dpx.pdb'))\
        #doctest: +ELLIPSIS
        Succeed to read file ...1dpx.pdb ,  1192 atoms found
        >>> rmsd, align_sel = prot_1_coor.align_seq_coor_to(prot_2_coor)\
        #doctest: +NORMALIZE_WHITESPACE
        ----------------------------------------NYFPQYPEYAIETARLRTFEAWPRNLKQKP--HQLAEAGF
                                                |   |  | | | | *|  | *  *  | \
*  ||*||
        KVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPC
        <BLANKLINE>
        FYTGVGDRVRCFSCGGGLMDWNDNDEPWEQHALWLSQCRFVKLMKGQLYIDTVAAKPV
             |* |   |*|  ||| ||| | *   | * ||*| | |   * |*    |  |
        SALLSSDITASVNCAKKIVSDGNGMNAW---VAWRNRCKGTDV---QAWI---RGCRL
        <BLANKLINE>
        >>> print('RMSD = {:.2f} Å'.format(rmsd))
        RMSD = 12.81 Å
        """

        sel_1_CA = self.select_part_dict(
            selec_dict={'name': ['CA'], 'chain': chain_1,
                        'res_name': PROTEIN_AA})
        sel_2_CA = atom_sel_2.select_part_dict(
            selec_dict={'name': ['CA'], 'chain': chain_2,
                        'res_name': PROTEIN_AA})

        sel_1_seq = sel_1_CA.get_aa_seq()
        sel_2_seq = sel_2_CA.get_aa_seq()

        seq_1 = ''
        for chain in chain_1:
            seq_1 += sel_1_seq[chain]
        seq_2 = ''
        for chain in chain_2:
            seq_2 += sel_2_seq[chain]

        align_seq_1, align_seq_2 = Coor.align_seq(seq_1, seq_2)
        Coor.print_align_seq(align_seq_1, align_seq_2)

        sel_index_1 = np.array(
            [key for key, atom in sorted(sel_1_CA.atom_dict.items())])
        sel_index_2 = np.array(
            [key for key, atom in sorted(sel_2_CA.atom_dict.items())])
        # print(len(sel_index_1), len(sel_index_2),
        # len(seq_1), len(seq_2), sel_index_1, sel_index_2)

        align_sel_1 = []
        align_sel_2 = []
        index_sel_1 = 0
        index_sel_2 = 0

        for i in range(len(align_seq_1)):
            # print(i, index_sel_1, index_sel_2)
            if align_seq_1[i] != '-' and align_seq_2[i] != '-':
                align_sel_1.append(sel_index_1[index_sel_1])
                align_sel_2.append(sel_index_2[index_sel_2])
                index_sel_1 += 1
                index_sel_2 += 1
            elif align_seq_1[i] != '-':
                index_sel_1 += 1
            else:
                index_sel_2 += 1

        self.align_to(atom_sel_2, index_list=[align_sel_1, align_sel_2])
        rmsd = self.compute_rmsd_to(
            atom_sel_2, index_list=[align_sel_1, align_sel_2])

        return rmsd, [align_sel_1, align_sel_2]

    @staticmethod
    def align_seq(seq_1, seq_2, gap_cost=-8, gap_extension=-2):
        """ Align two amino acid sequences using the Waterman - Smith Algorithm.

        :param seq_1: amino acid sequence 1
        :type seq_1: str

        :param seq_2: amino acid sequence 2
        :type seq_2: str

        :param gap_cost: Gap cost
        :type gap_cost: int (Default -8)

        :param gap_extension: Gap extension cost
        :type gap_extension: int (Default -2)

        :return: the two aligned sequences
        :rtype: str, str

        :Example:


        >>> seq_1 = 'AQDMVSPPPPIADEPLTVNTGIYLIECYSLDDKAETFKVNAFLSLSWKDRRLAFDPV\
RSGVRVKTYEPEAIWIPEIRFVNVENARDADVVDISVSPDGTVQYLERFSARVLSPLDFRRYPFDSQTLHIYLIVRSV\
DTRNIVLAVDLEKVGKNDDVFLTGWDIESFTAVVKPANFALEDRLESKLDYQLRISRQYFSYIPNIILPMLFILFISW\
TAFWSTSYEANVTLVVSTLIAHIAFNILVETNLPKTPYMTYTGAIIFMIYLFYFVAVIEVTVQHYLKVESQPARAASI\
TRASRIAFPVVFLLANIILAFLFFGF'
        >>> seq_2 = 'APSEFLDKLMGKVSGYDARIRPNFKGPPVNVTCNIFINSFGSIAETTMDYRVNIFLR\
QQWNDPRLAYSEYPDDSLDLDPSMLDSIWKPDLFFANEKGANFHEVTTDNKLLRISKNGNVLYSIRITLVLACPMDLK\
NFPMDVQTCIMQLESFGYTMNDLIFEWDEKGAVQVADGLTLPQFILKEEKDLRYCTKHYNTGKFTCIEARFHLERQMG\
YYLIQMYIPSLLIVILSWVSFWINMDAAPARVGLGITTVLTMTTQSSGSRASLPKVSYVKAIDIWMAVCLLFVFSALL\
EYAAVNFIARAGTKLFISRAKRIDTVSRVAFPLVFLIFNIFYWITY\
KLVPR'
        >>> align_seq_1, align_seq_2 = Coor.align_seq(seq_1, seq_2)
        >>> Coor.print_align_seq(align_seq_1, align_seq_2)\
        #doctest: +NORMALIZE_WHITESPACE
        -------------AQDMVSPPPPIADEPLTVNTGIYLIECYSLDDKAETFKVNAFLSLSWKDRRLAFDPVRS-GVRVKTY
                     |   |   * |   *||*| |*|| |  *| | |  ||** **  |*|* ***|\
|   | || |
        APSEFLDKLMGKVSGYDARIRPNFKGPPVNVTCNIFINSFGSIAETTMDYRVNIFLRQQWNDPRLAYSEYPDDSLDLDPS
        <BLANKLINE>
        EPEAIWIPEIRFVNVENARDADVV----DISVSPDGTVQYLERFSARVLSPLDFRRYPFDSQTLHIYLIVRSVDTRNIVL
          ||** *|| *|* ||*|  |*|     | |* |*|* *  *||  |  *|*||||*|* **  | \
*   |   |||||
        MLDSIWKPDLFFANEKGANFHEVTTDNKLLRISKNGNVLYSIRITLVLACPMDLKNFPMDVQTCIMQLESFGYTMNDLIF
        <BLANKLINE>
        AVDLEKVGKNDDVFLTGWDIESFTAVV-KPANFALEDRLESK---LDYQLRISRQYFSYIPNIILPMLFILFISWTAFW-
          * ** *  |     *  | |*     *  ||| |    |*   || |||||**   *| || |* *\
|*|||**||**
        EWD-EK-GAVQ--VADGLTLPQFILKEEKDLRYCTKHYNTGKFTCIEARFHLERQMGYYLIQMYIPSLLIVILSWVSFWI
        <BLANKLINE>
        -STSYEANVTLVVSTLIAHIAFNILVETNLPKTPYMTYTGAIIFMIYLFYFVAVIEVTVQHYL-KVESQP--ARAASITR
           |  *|* * ||*|||  | |   |||***| *|      | |  ** * *||* || ||| || |\
|   |**  *
        NMDAAPARVGLGITTVLTMTTQSSGSRASLPKVSYVKAIDIWMAVCLLFVFSALLEYAAVNFIARAGTKLFISRAKRIDT
        <BLANKLINE>
        ASRIAFPVVFLLANII--LAFLFFGF
        |**|***|***| **|  ||| |
        VSRVAFPLVFLIFNIFYWITYKLVPR
        <BLANKLINE>
        """

        seq_1_len = len(seq_1)
        seq_2_len = len(seq_2)

        direction_matrix = np.zeros((seq_1_len + 1, seq_2_len + 1), int)
        score_matrix = np.zeros((seq_1_len + 1, seq_2_len + 1), int)

        score_matrix[0] = 0
        score_matrix[:, 0] = 0

        # left = 1, top = 2, top-left = 3
        direction_matrix[0] = 1
        direction_matrix[:, 0] = 2
        direction_matrix[0, 0] = 0

        # Compute the score and direction matrix
        for i in range(1, seq_1_len + 1):
            for j in range(1, seq_2_len + 1):
                gap_top, gap_left = gap_cost, gap_cost
                # Check if previous gap is also a gap, then do an extension
                if direction_matrix[i, j - 1] == 1:
                    gap_top = gap_extension
                if direction_matrix[i - 1, j] == 2:
                    gap_left = gap_extension

                top = score_matrix[i, j - 1] + gap_top
                left = score_matrix[i - 1, j] + gap_left

                if (seq_1[i - 1], seq_2[j - 1]) in BLOSUM62:
                    diag = score_matrix[i - 1, j - 1] +\
                        BLOSUM62[seq_1[i - 1], seq_2[j - 1]]
                else:
                    diag = score_matrix[i - 1, j - 1] +\
                        BLOSUM62[seq_2[j - 1], seq_1[i - 1]]

                score_matrix[i, j] = max(top, left, diag, 0)

                # Save direction
                if score_matrix[i, j] == top:
                    direction_matrix[i, j] = 1
                elif score_matrix[i, j] == left:
                    direction_matrix[i, j] = 2
                else:
                    direction_matrix[i, j] = 3

        # Compute the sequence alignement from the direction matrix.
        i = seq_1_len
        j = seq_2_len

        seq_1_align = ''
        seq_2_align = ''

        while i != 0 or j != 0:
            if direction_matrix[i, j] == 3:
                seq_1_align = seq_1[i - 1] + seq_1_align
                seq_2_align = seq_2[j - 1] + seq_2_align
                i -= 1
                j -= 1
            elif direction_matrix[i, j] == 1:
                seq_1_align = '-' + seq_1_align
                seq_2_align = seq_2[j - 1] + seq_2_align
                j -= 1
            elif direction_matrix[i, j] == 2:
                seq_1_align = seq_1[i - 1] + seq_1_align
                seq_2_align = '-' + seq_2_align
                i -= 1

        return seq_1_align, seq_2_align

    @staticmethod
    def print_align_seq(seq_1, seq_2, line_len=80):

        sim_seq = ''
        for i in range(len(seq_1)):

            if seq_1[i] == seq_2[i]:
                sim_seq += '*'
                continue
            elif seq_1[i] != '-' and seq_2[i] != '-':
                if (seq_1[i], seq_2[i]) in BLOSUM62:
                    mut_score = BLOSUM62[seq_1[i], seq_2[i]]
                else:
                    mut_score = BLOSUM62[seq_2[i], seq_1[i]]
                if mut_score >= 0:
                    sim_seq += '|'
                    continue
            sim_seq += ' '

        for i in range(1 + len(seq_1) // line_len):
            print(seq_1[i * line_len: (i + 1) * line_len])
            print(sim_seq[i * line_len: (i + 1) * line_len])
            print(seq_2[i * line_len: (i + 1) * line_len])
            print()
        return

    @staticmethod
    def kabsch(coor_1, coor_2):
        """ Source: https://github.com/charnley/rmsd/blob/master/rmsd/\
        calculate_rmsd.py
        Using the Kabsch algorithm with two sets of paired point P and Q, \
        centered around the centroid. Each vector set is represented as an NxD
        matrix, where D is the the dimension of the space.

        The algorithm works in three steps:
        - a centroid translation of P and Q (assumed done before this
        function call)
        - the computation of a covariance matrix C
        - computation of the optimal rotation matrix U
        For more info see http://en.wikipedia.org/wiki/Kabsch_algorithm

        :param coor_1: coordinates array of size (N, D),
            where N is points and D is dimension.
        :type coor_1: np.array

        :param coor_2: coordinates array of size (N, D),
            where N is points and D is dimension.
        :type coor_2: np.array

        :return: rotation matrix
        :rtype: np.array of size (D, D)
        """

        # Computation of the covariance matrix
        C = np.dot(np.transpose(coor_1), coor_2)

        # Computation of the optimal rotation matrix
        # This can be done using singular value decomposition (SVD)
        # Getting the sign of the det(V)*(W) to decide
        # whether we need to correct our rotation matrix to ensure a
        # right-handed coordinate system.
        # And finally calculating the optimal rotation matrix U
        # see http://en.wikipedia.org/wiki/Kabsch_algorithm

        V, S, W = np.linalg.svd(C)
        d = (np.linalg.det(V) * np.linalg.det(W)) < 0.0

        if d:
            S[-1] = -S[-1]
            V[:, -1] = -V[:, -1]

        # Create Rotation matrix U
        rot_mat = np.dot(V, W)

        return rot_mat

    @staticmethod
    def quaternion_transform(r):
        """
        Source: https://github.com/charnley/rmsd/blob/master/rmsd/\
        calculate_rmsd.py
        Get optimal rotation
        note: translation will be zero when the centroids of each
        molecule are the same.
        """
        Wt_r = Coor.makeW(*r).T
        Q_r = Coor.makeQ(*r)
        rot = Wt_r.dot(Q_r)[:3, :3]
        return rot

    @staticmethod
    def makeW(r1, r2, r3, r4=0):
        """
        Source: https://github.com/charnley/rmsd/blob/master/rmsd/\
        calculate_rmsd.py
        matrix involved in quaternion rotation
        """
        W = np.asarray([
            [r4, r3, -r2, r1],
            [-r3, r4, r1, r2],
            [r2, -r1, r4, r3],
            [-r1, -r2, -r3, r4]])
        return W

    @staticmethod
    def makeQ(r1, r2, r3, r4=0):
        """
        Source: https://github.com/charnley/rmsd/blob/master/rmsd/\
        calculate_rmsd.py
        matrix involved in quaternion rotation
        """
        Q = np.asarray([
            [r4, -r3, r2, r1],
            [r3, r4, -r1, r2],
            [-r2, r1, r4, r3],
            [-r1, -r2, -r3, r4]])
        return Q

    @staticmethod
    def quaternion_rotate(X, Y):
        """
        Source: https://github.com/charnley/rmsd/blob/master/rmsd/\
        calculate_rmsd.py
        Calculate the rotation

        :param coor_1: coordinates array of size (N, D),\
            where N is points and D is dimension.
        :type coor_1: np.array

        :param coor_2: coordinates array of size (N, D),\
            where N is points and D is dimension.
        :type coor_2: np.array

        :return: rotation matrix
        :rtype: np.array of size (D, D)
        """
        N = X.shape[0]
        W = np.asarray([Coor.makeW(*Y[k]) for k in range(N)])
        Q = np.asarray([Coor.makeQ(*X[k]) for k in range(N)])
        Qt_dot_W = np.asarray([np.dot(Q[k].T, W[k]) for k in range(N)])
        # NOTE UNUSED W_minus_Q = np.asarray([W[k] - Q[k] for k in range(N)])
        A = np.sum(Qt_dot_W, axis=0)
        eigen = np.linalg.eigh(A)
        r = eigen[1][:, eigen[0].argmax()]
        rot = Coor.quaternion_transform(r)
        return rot

    def rotation_angle(self, tau_x, tau_y, tau_z):
        """ Compute coordinates of a system after a rotation on x, y and z axis.

        :param tau_x: angle of rotation (degrees) on the x axis
        :type tau_x: float

        :param tau_y: angle of rotation (degrees) on the y axis
        :type tau_y: float

        :param tau_z: angle of rotation (degrees) on the z axis
        :type tau_z: float

        :Example:

        >>> prot_coor = Coor(os.path.join(TEST_PATH, '1y0m.pdb'))\
        #doctest: +ELLIPSIS
        Succeed to read file ...1y0m.pdb ,  648 atoms found
        >>> com_1y0m = prot_coor.center_of_mass()
        >>> print("x:{:.2f} y:{:.2f} z:{:.2f}".format(*com_1y0m))
        x:16.01 y:0.45 z:8.57
        >>> prot_coor.rotation_angle(90, 90, 90)
        >>> com_1y0m = prot_coor.center_of_mass()
        >>> print("x:{:.2f} y:{:.2f} z:{:.2f}".format(*com_1y0m))\
        #doctest: +ELLIPSIS
        x:9.98 y:-4.03 z:-14.63

        """

        coor_array = self.get_array()

        tau_x = np.degrees(tau_x)
        tau_y = np.degrees(tau_y)
        tau_z = np.degrees(tau_z)

        x_rot_mat = np.array([[1, 0, 0],
                              [0, np.cos(tau_x), -np.sin(tau_x)],
                              [0, np.sin(tau_x), np.cos(tau_x)]])

        y_rot_mat = np.array([[np.cos(tau_y), 0, np.sin(tau_y)],
                              [0, 1, 0],
                              [-np.sin(tau_y), 0, np.cos(tau_y)]])

        z_rot_mat = np.array([[np.cos(tau_z), -np.sin(tau_z), 0],
                              [np.sin(tau_z), np.cos(tau_z), 0],
                              [0, 0, 1]])

        rotation_matrix = np.dot(np.dot(x_rot_mat, y_rot_mat), z_rot_mat)

        coor_array = np.dot(coor_array, rotation_matrix)

        for i, atom_num in enumerate(sorted(self.atom_dict)):
            self.atom_dict[atom_num]['xyz'] = coor_array[i]

        return

    def moment_inertia(self):
        """ Tensor moment of inertia relative to center of mass.

        Taken from:
        https://github.com/MDAnalysis/mdanalysis/blob/develop/package/MDAnalysis/core/topologyattrs.py

        :return: tensor matrix 3*3
        :rtype: np.array

        :Example:

        >>> prot_coor = Coor(os.path.join(TEST_PATH, '1y0m.pdb'))\
        #doctest: +ELLIPSIS
        Succeed to read file ...1y0m.pdb ,  648 atoms found
        >>> tensor = prot_coor.moment_inertia()
        >>> print(tensor)
        [[343695.37614973  31289.34303697  40416.71166176]
         [ 31289.34303697 478561.60329129  -9443.32698326]
         [ 40416.71166176  -9443.32698326 413840.08602987]]
        """

        coor_array = self.get_array()

        # Convert to local coordinates
        com = self.center_of_mass()

        coor_array -= com

        masses = self.get_mass_array()

        # Create the inertia tensor
        # m_i = mass of atom i
        # (x_i, y_i, z_i) = pos of atom i
        # Ixx = sum(m_i*(y_i^2+z_i^2));
        # Iyy = sum(m_i*(x_i^2+z_i^2));
        # Izz = sum(m_i*(x_i^2+y_i^2))
        # Ixy = Iyx = -1*sum(m_i*x_i*y_i)
        # Ixz = Izx = -1*sum(m_i*x_i*z_i)
        # Iyz = Izy = -1*sum(m_i*y_i*z_i)
        tens = np.zeros((3, 3), dtype=np.float64)
        # xx
        tens[0][0] = (
            masses * (coor_array[:, 1] ** 2 + coor_array[:, 2] ** 2)).sum()
        # xy & yx
        tens[0][1] = tens[1][0] = - (
            masses * coor_array[:, 0] * coor_array[:, 1]).sum()
        # xz & zx
        tens[0][2] = tens[2][0] = - (
            masses * coor_array[:, 0] * coor_array[:, 2]).sum()
        # yy
        tens[1][1] = (
            masses * (coor_array[:, 0] ** 2 + coor_array[:, 2] ** 2)).sum()
        # yz + zy
        tens[1][2] = tens[2][1] = - (
            masses * coor_array[:, 1] * coor_array[:, 2]).sum()
        # zz
        tens[2][2] = (
            masses * (coor_array[:, 0] ** 2 + coor_array[:, 1] ** 2)).sum()

        return tens

    def principal_axis(self):
        """ Calculate the principal axes from the moment of inertia.

        Taken from:
        https://github.com/MDAnalysis/mdanalysis/blob/develop/\
        package/MDAnalysis/core/topologyattrs.py

        :return: 3 principal axis
        :rtype: list of np.array


        :Example:

        >>> prot_coor = Coor(os.path.join(TEST_PATH, '1y0m.pdb'))\
        #doctest: +ELLIPSIS
        Succeed to read file ...1y0m.pdb ,  648 atoms found
        >>> princ_axis = prot_coor.principal_axis()
        >>> print(princ_axis)
        [[-0.21318784 -0.97697417  0.00850932]
         [ 0.39228202 -0.07761763  0.91656441]
         [-0.89479929  0.19873844  0.39979654]]
        """

        e_val, e_vec = np.linalg.eig(self.moment_inertia())

        # Sort
        indices = np.argsort(e_val)[::-1]
        # Return transposed in more logical form. See Issue 33.
        return e_vec[:, indices].T

    def align_principal_axis(self, axis=2,
                             vector=[0, 0, 1],
                             selec_dict={}):
        """Align principal axis with index `axis` with `vector`.

        Taken from:
        https://github.com/MDAnalysis/mdanalysis/blob/develop/\
package/MDAnalysis/core/topologyattrs.py


        :param axis: principal axis to align (0, 1, or 2)
        :type axis: int (Default 2)

        :param vector: vector to align
        :type vector: list (Default z axis)

        :param selec_dict: selection dictionnary
        :type selec_dict: dict, default={'name': ['CA']}

        :Example:

        >>> show_log()
        >>> prot_coor = Coor(os.path.join(TEST_PATH, '1y0m.pdb'))\
        #doctest: +ELLIPSIS
        Succeed to read file ...1y0m.pdb ,  648 atoms found
        >>> prot_coor.align_principal_axis()
        Do a rotation of 66.43°
        >>> prot_coor.align_principal_axis()
        Do a rotation of 0.00°
        """

        sele_dict = self.select_part_dict(selec_dict)

        p = sele_dict.principal_axis()[axis]
        # print(p, vector)
        angle = np.degrees(Coor.angle_vec(p, vector))
        logger.info('Do a rotation of {:.2f}°'.format(angle))

        r, rmsd = Rotation.align_vectors(
            p.reshape((1, 3)), np.array(vector).reshape((1, 3)))
        # print(r.as_matrix())

        coor_array = self.get_array()
        coor_array -= self.center_of_mass()
        coor_array = np.dot(coor_array, r.as_matrix())

        for i, atom_num in enumerate(sorted(self.atom_dict)):
            self.atom_dict[atom_num]['xyz'] = coor_array[i]

        return

    def get_max_size(self):
        """Get maximum size of a molecule.

        :Example:

        >>> prot_coor = Coor(os.path.join(TEST_PATH, '1y0m.pdb'))\
        #doctest: +ELLIPSIS
        Succeed to read file ...1y0m.pdb ,  648 atoms found
        >>> print('Maximum size without alignment is {:.2f} Å'.format(
        ... np.ceil(prot_coor.get_box_dim()).max()))
        Maximum size without alignment is 40.00 Å
        >>> max_size = prot_coor.get_max_size()
        Do a rotation of 66.43°
        >>> print('Maximum size is {:.2f} Å'.format(max_size))
        Maximum size is 43.00 Å

        """

        # Save initial coordinates
        coor_array = self.get_array()

        self.align_principal_axis(selec_dict={})
        max_dim = np.ceil(self.get_box_dim()).max()

        # Put back intial coordinates
        for i, atom_num in enumerate(sorted(self.atom_dict)):
            self.atom_dict[atom_num]['xyz'] = coor_array[i]

        return(max_dim)

    @staticmethod
    def angle_vec(vec_a, vec_b):
        """ Compute angle between two vectors.

        :param vec_a: vector
        :type vec_a: list

        :param vec_b: vector
        :type vec_b: list

        :return: angle in radian
        :rtype: float

        :Example:

        >>> angle = Coor.angle_vec([1, 0, 0], [0, 1, 0])
        >>> print('angle = {:.2f}'.format(np.degrees(angle)))
        angle = 90.00
        >>> angle = Coor.angle_vec([1, 0, 0], [1, 0, 0])
        >>> print('angle = {:.2f}'.format(np.degrees(angle)))
        angle = 0.00
        >>> angle = Coor.angle_vec([1, 0, 0], [1, 1, 0])
        >>> print('angle = {:.2f}'.format(np.degrees(angle)))
        angle = 45.00
        >>> angle = Coor.angle_vec([1, 0, 0], [-1, 0, 0])
        >>> print('angle = {:.2f}'.format(np.degrees(angle)))
        angle = 180.00

        """

        unit_vec_a = vec_a / np.linalg.norm(vec_a)
        unit_vec_b = vec_b / np.linalg.norm(vec_b)

        dot_product = np.dot(unit_vec_a, unit_vec_b)

        angle = np.arccos(dot_product)

        return(angle)

    def dist_under_index(self, atom_sel_2, cutoff=10.0):
        """ Check is distance between atoms of self.coor is under cutoff with
        atoms of group 1.
        Then return list of index of atoms of self.coor under ctuoff ditance.

        :param atom_sel_1: atom dictionnary
        :type atom_sel_1: dict

        :param atom_sel_2: atom dictionnary
        :type atom_sel_2: dict

        :param cutoff: distance cutoff
        :type cutoff: float, default=10.0

        :return: array of index
        :rtype: np.array

        """

        coor_array = self.get_array()
        index_array = np.array(
            [key for key, atom in sorted(self.atom_dict.items())])

        coor_array_2 = atom_sel_2.get_array()

        # Compute distance matrix
        dist_mat = distance_matrix(coor_array, coor_array_2)

        # Compute index of matrix column under cutoff_max and over cutoff_min:
        dist_mat_good = np.where((dist_mat.min(1) < cutoff))[0]

        return index_array[dist_mat_good]

    def make_peptide(self, sequence, pdb_out, check_file_out=True):
        """
        Create a linear peptide structure.

        :param sequence: peptide sequence
        :type sequence: str

        :param pdb_out: name of output pdb file
        :type pdb_out: str

        :param check_file_out: flag to check or not if file has
            already been created. If the file is present then the
            command break.
        :type check_file_out: bool, optional, default=True

        """

        # Create and go in out_folder:
        # This is necessary for the topologie creation
        out_folder = os.path.dirname(pdb_out)
        # print(out_folder)
        os_command.create_dir(out_folder)

        logger.info("-Make peptide: {}".format(sequence))

        # Check if output files exist:
        if check_file_out and os_command.check_file_and_create_path(pdb_out):
            logger.info("make_peptide %s already exist" % pdb_out)
            return os.path.abspath(pdb_out)

        pep = Coor()
        seq = 'X' + sequence
        atom_num = 0
        uniq_resid = 0
        connect_dict = {}
        prev_res_name_index = {}

        # Initialize atom_dict:
        for res_name in seq:
            logger.info("residue name:{}".format(res_name))
            res_name_index = {}

            for atom_name in AA_ATOM_DICT[res_name]:
                # print("\tatom name:{}".format(atom_name))

                if atom_num == 0:
                    xyz = np.zeros(3)
                else:
                    # Look for the previous bonded atom:
                    for dist in AA_BOND_DICT[res_name]:
                        if atom_name == dist[0]:
                            if dist[1][0] != '-':
                                connect_name = dist[1]
                                connect_index = res_name_index[connect_name]
                            else:
                                connect_name = dist[1][1:]
                                connect_index =\
                                    prev_res_name_index[connect_name]
                            break
                        elif atom_name == dist[1]:
                            if dist[0][0] != '-':
                                connect_name = dist[0]
                                connect_index = res_name_index[connect_name]
                            else:
                                connect_name = dist[0][1:]
                                connect_index =\
                                    prev_res_name_index[connect_name]
                            break
                    # print("{} connect to {} for {}".format(
                    #    atom_name, connect_name, res_name))
                    bond_len = Coor.find_dist(res_name,
                                              atom_name,
                                              connect_name)
                    # print("Bond : {}-{} = {} X".format(
                    #    atom_name, connect_name, bond_len))

                    if atom_num == 1:
                        xyz = pep.atom_dict[connect_index]['xyz'] +\
                            [bond_len, 0, 0]
                    connect_dict[atom_num] = connect_index

                    if atom_num > 1:
                        connect_2_index = connect_dict[connect_index]
                        connect_2_name = pep.atom_dict[connect_2_index]['name']
                        angle = Coor.find_angle(res_name, atom_name,
                                                connect_name, connect_2_name)
                        angle_rad = np.deg2rad(angle)
                        # print("Angle: {}-{}-{} = {}°".format(
                        # atom_name, connect_name, connect_2_name, angle))
                        if atom_num == 2:
                            xyz = pep.atom_dict[connect_index]['xyz'] + [
                                -bond_len * np.cos(np.deg2rad(angle)),
                                -bond_len * np.sin(np.deg2rad(angle)),
                                0]

                    if atom_num > 2:
                        if connect_2_index not in connect_dict:
                            xyz = pep.atom_dict[connect_index]['xyz'] +\
                                [-bond_len * np.cos(np.deg2rad(angle)),
                                 -bond_len * np.sin(np.deg2rad(angle)), 0]
                        else:
                            connect_3_index = connect_dict[connect_2_index]
                            connect_3_name =\
                                pep.atom_dict[connect_3_index]['name']
                            dihe = Coor.find_dihe_angle(res_name, atom_name,
                                                        connect_name,
                                                        connect_2_name,
                                                        connect_3_name)
                            dihe_rad = np.deg2rad(dihe)
                            # print("Dihedral Angle: {}-{}-{}-{} = {}°".format(
                            # atom_name, connect_name, connect_2_name,
                            # connect_3_name, dihe))
                            # From https://github.com/ben-albrecht/qcl/blob/\
                            # master/qcl/ccdata_xyz.py#L208
                            vec_1 = pep.atom_dict[connect_index]['xyz'] -\
                                pep.atom_dict[connect_2_index]['xyz']
                            vec_2 = pep.atom_dict[connect_index]['xyz'] -\
                                pep.atom_dict[connect_3_index]['xyz']

                            vec_n = np.cross(vec_1, vec_2)
                            vec_nn = np.cross(vec_1, vec_n)

                            vec_n /= norm(vec_n)
                            vec_nn /= norm(vec_nn)

                            vec_n *= -sin(dihe_rad)
                            vec_nn *= cos(dihe_rad)

                            vec_3 = vec_n + vec_nn
                            vec_3 /= norm(vec_3)
                            vec_3 *= bond_len * sin(angle_rad)

                            vec_1 /= norm(vec_1)
                            vec_1 *= bond_len * cos(angle_rad)

                            xyz = pep.atom_dict[connect_index]['xyz'] +\
                                vec_3 - vec_1
                            # print(xyz)

                atom = {"field": 'ATOM',
                        "num": atom_num,
                        "name": atom_name,
                        "alter_loc": "",
                        "res_name": AA_1_TO_3_DICT[res_name],
                        "chain": "P",
                        "res_num": uniq_resid,
                        "uniq_resid": uniq_resid,
                        "insert_res": "",
                        "xyz": xyz,
                        "occ": 0.0,
                        "beta": 0.0,
                        "elem": ''}
                res_name_index[atom_name] = atom_num
                # print(atom)
                pep.atom_dict[atom_num] = atom
                atom_num += 1
            prev_res_name_index = res_name_index
            uniq_resid += 1

        pep.write_pdb(pdb_out)
        pdb_out = os.path.abspath(pdb_out)

        return pdb_out

    @staticmethod
    def find_dist(aa_name, name_a, name_b):
        for dist in DIST_DICT[aa_name]:
            if dist[:2] == [name_a, name_b] or dist[:2] == [name_b, name_a]:
                return dist[2]
        raise ValueError('Distance param {}-{} for {} not found !!'.format(
            name_a, name_b, aa_name))

    @staticmethod
    def find_angle(aa_name, name_a, name_b, name_c):
        for angle in ANGLE_DICT[aa_name]:
            if (angle[:3] == [name_a, name_b, name_c] or
                    angle[:3] == [name_c, name_b, name_a]):
                return angle[3]
        raise ValueError('Angle param {}-{}-{} for {} not found !!'.format(
            name_a, name_b, name_c, aa_name))

    @staticmethod
    def find_dihe_angle(aa_name, name_a, name_b, name_c, name_d):
        for angle in DIHE_DICT[aa_name]:
            if (angle[:4] == [name_a, name_b, name_c, name_d] or
                    angle[:4] == [name_d, name_c, name_b, name_a]):
                return angle[4]
        raise ValueError(
            'Angle param {}-{}-{}-{} for {} not found !!'.format(
                name_a, name_b, name_c, name_d, aa_name))

    @staticmethod
    def atom_dist(atom_a, atom_b):
        """Compute the distance between 2 atoms.

        :param atom_a: atom dictionnary
        :type atom_a: dict

        :param atom_b: atom dictionnary
        :type atom_b: dict

        :return: distance
        :rtype: float

        :Example:

        >>> atom_1 = {'xyz': np.array([0.0, 0.0, 0.0])}
        >>> atom_2 = {'xyz': np.array([0.0, 1.0, 0.0])}
        >>> atom_3 = {'xyz': np.array([1.0, 1.0, 1.0])}
        >>> Coor.atom_dist(atom_1, atom_2)
        1.0
        >>> Coor.atom_dist(atom_1, atom_3)
        1.7320508075688772
        """

        distance = np.linalg.norm(atom_a['xyz'] - atom_b['xyz'])
        return distance

    @staticmethod
    def atom_angle(atom_a, atom_b, atom_c):
        """Compute the anlge between 3 atoms.

        :param atom_a: atom dictionnary
        :type atom_a: dict

        :param atom_b: atom dictionnary
        :type atom_b: dict

        :param atom_c: atom dictionnary
        :type atom_c: dict

        :return: angle (degrees)
        :rtype: float

        :Example:

        >>> atom_1 = {'xyz': np.array([0.0, 0.0, 0.0])}
        >>> atom_2 = {'xyz': np.array([0.0, 1.0, 0.0])}
        >>> atom_3 = {'xyz': np.array([1.0, 1.0, 1.0])}
        >>> Coor.atom_angle(atom_1, atom_2, atom_3)
        90.0
        >>> print('{:.3f}'.format(Coor.atom_angle(atom_1, atom_3, atom_2)))
        35.264
        """

        ba = atom_a['xyz'] - atom_b['xyz']
        bc = atom_c['xyz'] - atom_b['xyz']
        angle = Coor.angle_vec(ba, bc)

        return np.degrees(angle)

    @staticmethod
    def atom_dihed_angle(atom_a, atom_b, atom_c, atom_d):
        """Compute the dihedral anlge using 4 atoms.

        :param atom_a: atom dictionnary
        :type atom_a: dict

        :param atom_b: atom dictionnary
        :type atom_b: dict

        :param atom_c: atom dictionnary
        :type atom_c: dict

        :param atom_d: atom dictionnary
        :type atom_d: dict

        :return: dihedral angle
        :rtype: float

        :Example:

        >>> atom_1 = {'xyz': np.array([0.0, -1.0, 0.0])}
        >>> atom_2 = {'xyz': np.array([0.0, 0.0, 0.0])}
        >>> atom_3 = {'xyz': np.array([1.0, 0.0, 0.0])}
        >>> atom_4 = {'xyz': np.array([1.0, 1.0, 0.0])}
        >>> atom_5 = {'xyz': np.array([1.0, -1.0, 0.0])}
        >>> atom_6 = {'xyz': np.array([1.0, -1.0, 1.0])}
        >>> angle_1 = Coor.atom_dihed_angle(atom_1, atom_2, atom_3, atom_4)
        >>> print('{:.3f}'.format(angle_1))
        180.000
        >>> angle_2 = Coor.atom_dihed_angle(atom_1, atom_2, atom_3, atom_5)
        >>> print('{:.3f}'.format(angle_2))
        0.000
        >>> angle_3 = Coor.atom_dihed_angle(atom_1, atom_2, atom_3, atom_6)
        >>> print('{:.3f}'.format(angle_3))
        -45.000
        """

        ab = -1 * (atom_b['xyz'] - atom_a['xyz'])
        bc = atom_c['xyz'] - atom_b['xyz']
        cd = atom_d['xyz'] - atom_c['xyz']

        v1 = np.cross(ab, bc)
        v2 = np.cross(cd, bc)
        v1_x_v2 = np.cross(v1, v2)

        y = np.dot(v1_x_v2, bc)*(1.0/np.linalg.norm(bc))
        x = np.dot(v1, v2)
        angle = np.arctan2(y, x)

        return np.degrees(angle)

    @staticmethod
    def concat_pdb(*pdb_in_files, pdb_out):
        """Concat a list of pdb files in one.

        :param pdb_in_files: list of pdb files
        :type pdb_in_files: list

        :param pdb_out: atom dictionnary
        :type pdb_out: dict

        :Example:

        >>> TEST_OUT = str(getfixture('tmpdir'))
        >>> Coor.concat_pdb(os.path.join(TEST_PATH, '1y0m.pdb'),
        ...                 os.path.join(TEST_PATH, '1rxz.pdb'),
        ...                 pdb_out = os.path.join(TEST_OUT, 'tmp_2.pdb'))\
         #doctest: +ELLIPSIS
        Succeed to save concat file: ...tmp_2.pdb
        """

        if os_command.check_file_and_create_path(pdb_out):
            logger.info("File  %s  already exist" % pdb_out)
            return

        filout = open(pdb_out, 'w')
        count = 0

        for pdb_in in pdb_in_files:
            # print("Concat :", os.path.relpath(pdb_in))
            with open(pdb_in) as pdbfile:
                for line in pdbfile:
                    if ((count == 0 and line[:6] == "CRYST1") or
                            line[:4] == 'ATOM' or line[:6] == "HETATM"):
                        filout.write(line)
            count += 1

        filout.close()
        logger.info("Succeed to save concat file:  %s" % pdb_out)


class Multi_Coor:
    """ Topologie base on coordinates like pdb or gro.

    The coor object containt a a list of dictionnary of atoms indexed
    on the atom num and the crystal packing info.


    :param coor_list: list of dictionnary of atom
    :type coor_listt: list

    :param crystal_pack: crystal packing
    :type crystal_pack: str

    .. note::

        Files necessary for testing :
        * ../test/input/1y0m.pdb, ../test/input/1rxz.pdb\
        and ../test/input/4n1m.pdb.
        To do the unitary test, execute pdb_mani.py (-v for verbose mode)

    """

    def __init__(self, pdb_in=None, pqr_format=False):
        self.coor_list = []
        self.crystal_pack = None

        if pdb_in is not None:
            self.read_pdb(pdb_in, pqr_format)

    def read_pdb(self, pdb_in, pqr_format=False):
        """Read a pdb file and return atom informations as a dictionnary
        indexed on the atom num.
        The fonction can also read pqr files if specified with
        ``pqr_format = True``, it will only change the column
        format of beta and occ factors.

        :param pdb_in: path of the pdb file to read
        :type pdb_in: str

        :param pqr_format: Flag for .pqr file format reading.
        :type pqr_format: bool, default=False

        :Example:

        >>> VIP_coor = Multi_Coor(os.path.join(TEST_PATH, '2rri.pdb'))\
        #doctest: +ELLIPSIS
        Read 20 Model(s)
        Succeed to read file ...2rri.pdb, 479 atoms found
        """

        atom_index = 0
        uniq_resid = -1
        old_res_num = -1
        model_num = 1

        model_coor = Coor()

        with open(pdb_in) as pdbfile:
            for line in pdbfile:
                if line.startswith("CRYST1"):
                    self.crystal_pack = line
                if line.startswith("MODEL"):
                    # print('Read Model {}'.format(model_num))
                    model_num += 1
                if line.startswith("ENDMDL"):
                    if model_coor.num != 0:
                        self.coor_list.append(model_coor)
                        model_coor = Coor()
                        atom_index = 0
                        uniq_resid = -1
                        old_res_num = -1

                if line.startswith('ATOM') or line.startswith("HETATM"):

                    field = line[:6].strip()
                    atom_num = int(line[6:11])
                    atom_name = line[12:16].strip()

                    res_name = line[17:20].strip()
                    chain = line[21]
                    res_num = int(line[22:26])
                    insert_res = line[26:27]
                    xyz = np.array([float(line[30:38]),
                                    float(line[38:46]),
                                    float(line[46:54])])
                    if pqr_format:
                        alter_loc = ""
                        res_name = line[16:20].strip()
                        occ, beta = line[54:62].strip(), line[62:70].strip()
                        elem_symbol = ""
                    else:
                        alter_loc = line[16:17]
                        res_name = line[17:20].strip()
                        occ, beta = line[54:60].strip(), line[60:66].strip()
                        elem_symbol = line[76:78]

                    if occ == "":
                        occ = 0.0
                    else:
                        occ = float(occ)

                    if beta == "":
                        beta = 0.0
                    else:
                        beta = float(beta)

                    if res_num != old_res_num:
                        uniq_resid += 1
                        old_res_num = res_num

                    atom = {"field": field,
                            "num": atom_num,
                            "name": atom_name,
                            "alter_loc": alter_loc,
                            "res_name": res_name,
                            "chain": chain,
                            "res_num": res_num,
                            "uniq_resid": uniq_resid,
                            "insert_res": insert_res,
                            "xyz": xyz,
                            "occ": occ,
                            "beta": beta,
                            "elem": elem_symbol}

                    model_coor.atom_dict[atom_index] = atom
                    atom_index += 1

        logger.info('Read {} Model(s)'.format(len(self.coor_list)))
        logger.info("Succeed to read file {}, {} atoms found".format(
            os.path.relpath(pdb_in),
            self.coor_list[0].num))

    def write_pdb(self, pdb_out, check_file_out=True):
        """Write a pdb file.

        :param pdb_out: path of the pdb file to write
        :type pdb_out: str

        :param check_file_out: flag to check if output file already exists
        :type check_file_outt: bool (Default True)

        :Example:

        >>> TEST_OUT = str(getfixture('tmpdir'))
        >>> VIP_coor = Multi_Coor(os.path.join(TEST_PATH, '2rri.pdb'))\
        #doctest: +ELLIPSIS
        Read 20 Model(s)
        Succeed to read file ...2rri.pdb, 479 atoms found
        >>> VIP_coor.write_pdb(os.path.join(TEST_OUT, 'tmp.pdb'))\
        #doctest: +ELLIPSIS
        Succeed to save file ...tmp.pdb
        """

        if check_file_out and os_command.check_file_and_create_path(pdb_out):
            logger.info("PDB file {:>4} already exist, file not saved".format(
                pdb_out))
            return

        filout = open(pdb_out, 'w')
        if self.crystal_pack is not None:
            filout.write(self.crystal_pack)

        model = 1
        for model_coor in self.coor_list:
            filout.write("MODEL      {:^4}\n".format(model))
            model += 1

            for atom_num, atom in sorted(model_coor.atom_dict.items()):
                # print(pdb_dict[atom_num]["name"])

                # Atom name should start a column 14, with the type of atom ex:
                #   - with atom type 'C': ' CH3'
                # for 2 letters atom type, it should start at coulumn 13 ex:
                #   - with atom type 'FE': 'FE1'
                name = atom["name"]
                if len(name) <= 3 and name[0] in ['C', 'H', 'O',
                                                  'N', 'S', 'P']:
                    name = " " + name

                filout.write(
                    "{:6s}{:5d} {:4s}{:1s}{:3s} {:1s}{:4d}{:1s}"
                    "   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}\n".format(
                        atom["field"],
                        atom["num"],
                        name,
                        atom["alter_loc"],
                        atom["res_name"],
                        atom["chain"],
                        atom["res_num"],
                        atom["insert_res"],
                        atom["xyz"][0],
                        atom["xyz"][1],
                        atom["xyz"][2],
                        atom["occ"],
                        atom["beta"],
                        atom["elem"]))
            filout.write("ENDMDL\n".format())

        filout.write("TER\n")
        filout.close()

        logger.info("Succeed to save file %s" % os.path.relpath(pdb_out))
        return

    def compute_rmsd_to(self, atom_sel_2, selec_dict={'name': ['CA']}):
        """ Compute RMSD between two atom_dict
        Then return the RMSD value.

        :param atom_sel_1: atom dictionnary
        :type atom_sel_1: dict

        :param atom_sel_2: atom dictionnary
        :type atom_sel_2: dict

        :return: distance
        :rtype: float

        :Example:

        >>> TEST_OUT = str(getfixture('tmpdir'))
        >>> VIP_coor = Multi_Coor(os.path.join(TEST_PATH, '2rri.pdb'))\
        #doctest: +ELLIPSIS
        Read 20 Model(s)
        Succeed to read file ...2rri.pdb, 479 atoms found
        >>> aa_seq = VIP_coor.coor_list[0].get_aa_seq()['A']
        >>> print(aa_seq)
        HSDAVFTDNYTRLRKQMAVKKYLNSILNG
        >>> linear_pep = Coor()
        >>> pep_out = os.path.join(TEST_OUT, 'tmp_pep.pdb')
        >>> pdb_pep = linear_pep.make_peptide(aa_seq, pep_out)\
        #doctest: +ELLIPSIS
        -Make peptide: HSDAVFTDNYTRLRKQMAVKKYLNSILNG
        residue name:X
        residue name:H
        residue name:S
        residue name:D
        residue name:A
        residue name:V
        residue name:F
        residue name:T
        residue name:D
        residue name:N
        residue name:Y
        residue name:T
        residue name:R
        residue name:L
        residue name:R
        residue name:K
        residue name:Q
        residue name:M
        residue name:A
        residue name:V
        residue name:K
        residue name:K
        residue name:Y
        residue name:L
        residue name:N
        residue name:S
        residue name:I
        residue name:L
        residue name:N
        residue name:G
        Succeed to save file ...tmp_pep.pdb
        >>> linear_pep.read_pdb(pdb_pep) #doctest: +ELLIPSIS
        Succeed to read file ...tmp_pep.pdb ,  240 atoms found
        >>> rmsd_list = VIP_coor.compute_rmsd_to(linear_pep)
        >>> rmsd_str = ['{:.2f}'.format(i) for i in rmsd_list]
        >>> rmsd_str
        ['58.57', '58.40', '58.74', '58.35', '58.60', '58.53',\
 '58.49', '58.40', '58.45', '58.27', '58.52', '58.34', '58.57',\
 '58.33', '58.34', '58.63', '58.61', '58.40', '58.55', '58.32']
        """

        rmsd_list = []

        for atom_coor in self.coor_list:

            rmsd = atom_coor.compute_rmsd_to(atom_sel_2, selec_dict=selec_dict)

            rmsd_list.append(rmsd)

        return rmsd_list


if __name__ == "__main__":

    import doctest
    import shutil

    TEST_DIR = 'pdb_manip_test_out'
    TEST_OUT = os.path.join(TEST_DIR, 'pdb_manip_test')

    def getfixture(*args):
        return TEST_OUT

    print("-Test pdb_manip module:")
    print("pdb_manip:\t", doctest.testmod())

    # Erase all test files
    shutil.rmtree(TEST_DIR, ignore_errors=True)
