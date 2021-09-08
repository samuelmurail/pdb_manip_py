#!/usr/bin/env python3
# coding: utf-8

"""
Location of datafiles for gromacs_py tests

Use it as:

```
from pdb_manip_py.test.datafiles import *
```
"""

import os

# Autorship information
__author__ = "Samuel Murail, Pierre Tuffery"
__copyright__ = "Copyright 2020, RPBS"
__credits__ = ["Samuel Murail"]
__license__ = "GNU General Public License v2.0"
__maintainer__ = "Samuel Murail"
__email__ = "samuel.murail@u-paris.fr"
__status__ = "Production"

PYTEST_DIR = os.path.dirname(os.path.abspath(__file__))
TEST_FILE_PATH = os.path.join(PYTEST_DIR, "input/")

PDB_1D30 = os.path.join(TEST_FILE_PATH, '1D30.pdb')
PDB_4N1M = os.path.join(TEST_FILE_PATH, '4n1m.pdb')
