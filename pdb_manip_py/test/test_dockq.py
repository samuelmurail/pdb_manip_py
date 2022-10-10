# coding: utf-8

"""
Tests for DockQ functions
"""

import os
import pytest

from pdb_manip_py import pdb_manip
from pdb_manip_py import dockq

from .datafiles import PDB_1JD4, PDB_5MN6, PDB_1RXZ, PDB_1RXZ_MODEL

# Autorship information
__author__ = "Samuel Murail, Pierre Tuffery"
__copyright__ = "Copyright 2020, RPBS"
__credits__ = ["Samuel Murail"]
__license__ = "GNU General Public License v2.0"
__maintainer__ = "Samuel Murail"
__email__ = "samuel.murail@u-paris.fr"
__status__ = "Production"


def test_dockq_bad(tmp_path):
    """

    Raw DockQ results:
    ****************************************************************
    *                       DockQ                                  *
    *   Scoring function for protein-protein docking models        *
    *   Statistics on CAPRI data:                                  *
    *    0.00 <= DockQ <  0.23 - Incorrect                         *
    *    0.23 <= DockQ <  0.49 - Acceptable quality                *
    *    0.49 <= DockQ <  0.80 - Medium quality                    *
    *            DockQ >= 0.80 - High quality                      *
    *   Reference: Sankar Basu and Bjorn Wallner, DockQ: A quality *
    *   measure for protein-protein docking models, submitted      *
    *                                                              *
    *   For the record:                                            *
    *   Definition of contact <5A (Fnat)                           *
    *   Definition of interface <10A all heavy atoms (iRMS)        *
    *   For comments, please email: bjorn.wallner@.liu.se          *
    *                                                              *
    ****************************************************************
    Model  : pdb_manip_py/test/input/1jd4.pdb
    Native : pdb_manip_py/test/input/5m6n.pdb
    Number of equivalent residues in chain A 96 (receptor)
    Number of equivalent residues in chain B 96 (ligand)
    Fnat 0.000 0 correct of 33 native contacts
    Fnonnat 1.000 45 non-native of 45 model contacts
    iRMS 15.631
    LRMS 59.981
    DockQ 0.010 

    """

    model_coor = pdb_manip.Coor(PDB_1JD4)
    native_coor = pdb_manip.Coor(PDB_5MN6)

    dockq = model_coor.compute_dockQ(native_coor)

    assert pytest.approx(dockq['DockQ'], 0.5) == 0.010
    assert dockq['Fnat'] == 0.0
    assert dockq['Fnonnat'] == 1.0
    assert pytest.approx(dockq['LRMS'], 0.1) == 59.981
    assert pytest.approx(dockq['iRMS'], 0.5) == 15.631

def test_dockq_good(tmp_path):
    """

    Raw DockQ results:
    ****************************************************************
    *                       DockQ                                  *
    *   Scoring function for protein-protein docking models        *
    *   Statistics on CAPRI data:                                  *
    *    0.00 <= DockQ <  0.23 - Incorrect                         *
    *    0.23 <= DockQ <  0.49 - Acceptable quality                *
    *    0.49 <= DockQ <  0.80 - Medium quality                    *
    *            DockQ >= 0.80 - High quality                      *
    *   Reference: Sankar Basu and Bjorn Wallner, DockQ: A quality *
    *   measure for protein-protein docking models, submitted      *
    *                                                              *
    *   For the record:                                            *
    *   Definition of contact <5A (Fnat)                           *
    *   Definition of interface <10A all heavy atoms (iRMS)        *
    *   For comments, please email: bjorn.wallner@.liu.se          *
    *                                                              *
    ****************************************************************
    Model  : test.pdb
    Native : ../pdb_manip_py/test/input/1rxz.pdb
    Number of equivalent residues in chain A 245 (receptor)
    Number of equivalent residues in chain B 11 (ligand)
    Fnat 0.963 52 correct of 54 native contacts
    Fnonnat 0.088 5 non-native of 57 model contacts
    iRMS 0.618
    LRMS 1.050
    DockQ 0.934 


    """

    model_coor = pdb_manip.Coor(PDB_1RXZ_MODEL)
    native_coor = pdb_manip.Coor(PDB_1RXZ)

    dockq = model_coor.compute_dockQ(native_coor)

    assert pytest.approx(dockq['DockQ'], 0.001) == 0.934
    assert pytest.approx(dockq['Fnat'], 0.001) == 0.963
    assert pytest.approx(dockq['Fnonnat'], 0.01) == 0.088
    assert pytest.approx(dockq['LRMS'], 0.001) == 1.050
    assert pytest.approx(dockq['iRMS'], 0.001) == 0.618

def test_pdockq(tmp_path):
    """

    Raw DockQ results:
    ****************************************************************
    *                       DockQ                                  *
    *   Scoring function for protein-protein docking models        *
    *   Statistics on CAPRI data:                                  *
    *    0.00 <= DockQ <  0.23 - Incorrect                         *
    *    0.23 <= DockQ <  0.49 - Acceptable quality                *
    *    0.49 <= DockQ <  0.80 - Medium quality                    *
    *            DockQ >= 0.80 - High quality                      *
    *   Reference: Sankar Basu and Bjorn Wallner, DockQ: A quality *
    *   measure for protein-protein docking models, submitted      *
    *                                                              *
    *   For the record:                                            *
    *   Definition of contact <5A (Fnat)                           *
    *   Definition of interface <10A all heavy atoms (iRMS)        *
    *   For comments, please email: bjorn.wallner@.liu.se          *
    *                                                              *
    ****************************************************************
    Model  : test.pdb
    Native : ../pdb_manip_py/test/input/1rxz.pdb
    Number of equivalent residues in chain A 245 (receptor)
    Number of equivalent residues in chain B 11 (ligand)
    Fnat 0.963 52 correct of 54 native contacts
    Fnonnat 0.088 5 non-native of 57 model contacts
    iRMS 0.618
    LRMS 1.050
    DockQ 0.934 


    """

    model_coor = pdb_manip.Coor(PDB_1RXZ_MODEL)

    pdockq = model_coor.compute_pdockQ()

    assert pytest.approx(pdockq, 0.001) == 0.535
