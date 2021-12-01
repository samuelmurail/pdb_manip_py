#!/usr/bin/env python3
# coding: utf-8

"""
Tests for PDB2PQR functions
"""

import os


from pdb_manip_py import pdb_manip
from pdb_manip_py import pdb2pqr

from .datafiles import PDB_1D30, PDB_4N1M

# Autorship information
__author__ = "Samuel Murail, Pierre Tuffery"
__copyright__ = "Copyright 2020, RPBS"
__credits__ = ["Samuel Murail"]
__license__ = "GNU General Public License v2.0"
__maintainer__ = "Samuel Murail"
__email__ = "samuel.murail@u-paris.fr"
__status__ = "Production"


def test_pdb2pqr_ph4_7_12_charmm(tmp_path):
    """
    Test pdb2pqr at ph 4, 7 and 12 with the charmm ff.

    """

    ff = "CHARMM"

    pdb_manip.show_log()
    ##########################################
    # ##   Compute protonation  at PH 4.0  ###
    ##########################################
    out_file = os.path.join(tmp_path, '4N1M_charmm_ph4.pqr')
    pdb2pqr.compute_pdb2pqr(os.path.join(PDB_4N1M),
                            out_file, ff=ff, ph=4.0)
    prot_coor = pdb_manip.Coor(out_file)
    assert prot_coor.num == 2559

    HSD_index = prot_coor.get_index_selection(
        {'res_name': ['HSD'], 'name': ['CA']})
    HSE_index = prot_coor.get_index_selection(
        {'res_name': ['HSE'], 'name': ['CA']})
    HSP_index = prot_coor.get_index_selection(
        {'res_name': ['HSP'], 'name': ['CA']})

    assert len(HSD_index) == 0
    assert len(HSE_index) == 0
    assert len(HSP_index) == 5

    ##########################################
    # ##   Compute protonation  at PH 7.0  ###
    ##########################################
    out_file = os.path.join(tmp_path, '4N1M_charmm_ph7.pqr')
    pdb2pqr.compute_pdb2pqr(os.path.join(PDB_4N1M),
                            out_file, ff=ff, ph=7.0)
    prot_coor = pdb_manip.Coor(out_file)
    assert prot_coor.num == 2549

    HSD_index = prot_coor.get_index_selection(
        {'res_name': ['HSD'], 'name': ['CA']})
    HSE_index = prot_coor.get_index_selection(
        {'res_name': ['HSE'], 'name': ['CA']})
    HSP_index = prot_coor.get_index_selection(
        {'res_name': ['HSP'], 'name': ['CA']})

    assert len(HSD_index) == 4
    assert len(HSE_index) == 0
    assert len(HSP_index) == 1

    ##########################################
    # ##   Compute protonation  at PH 10.0 ###
    out_file = os.path.join(tmp_path, '4N1M_charmm_ph10.pqr')
    pdb2pqr.compute_pdb2pqr(os.path.join(PDB_4N1M),
                            out_file, ff=ff, ph=10.0)
    prot_coor = pdb_manip.Coor(out_file)
    assert prot_coor.num == 2548

    HSD_index = prot_coor.get_index_selection(
        {'res_name': ['HSD'], 'name': ['CA']})
    HSE_index = prot_coor.get_index_selection(
        {'res_name': ['HSE'], 'name': ['CA']})
    HSP_index = prot_coor.get_index_selection(
        {'res_name': ['HSP'], 'name': ['CA']})

    assert len(HSD_index) == 5
    assert len(HSE_index) == 0
    assert len(HSP_index) == 0


def test_pdb2pqr_ph4_7_12_amber(tmp_path):
    """
    Test pdb2pqr at ph 4, 7 and 12 with the amber ff.

    """

    ff = "AMBER"

    pdb_manip.show_log()
    ##########################################
    # ##   Compute protonation  at PH 4.0  ###
    ##########################################
    out_file = os.path.join(tmp_path, '4N1M_amber_ph4.pqr')
    pdb2pqr.compute_pdb2pqr(os.path.join(PDB_4N1M),
                            out_file, ff=ff, ph=4.0)
    prot_coor = pdb_manip.Coor(out_file)
    assert prot_coor.num == 2559

    HSD_index = prot_coor.get_index_selection(
        {'res_name': ['HID'], 'name': ['CA']})
    HSE_index = prot_coor.get_index_selection(
        {'res_name': ['HIE'], 'name': ['CA']})
    HSP_index = prot_coor.get_index_selection(
        {'res_name': ['HIP'], 'name': ['CA']})

    assert len(HSD_index) == 0
    assert len(HSE_index) == 0
    assert len(HSP_index) == 5

    ##########################################
    # ##   Compute protonation  at PH 7.0  ###
    ##########################################
    out_file = os.path.join(tmp_path, '4N1M_amber_ph7.pqr')
    pdb2pqr.compute_pdb2pqr(os.path.join(PDB_4N1M),
                            out_file, ff=ff, ph=7.0)
    prot_coor = pdb_manip.Coor(out_file)
    assert prot_coor.num == 2549

    HSD_index = prot_coor.get_index_selection(
        {'res_name': ['HID'], 'name': ['CA']})
    HSE_index = prot_coor.get_index_selection(
        {'res_name': ['HIE'], 'name': ['CA']})
    HSP_index = prot_coor.get_index_selection(
        {'res_name': ['HIP'], 'name': ['CA']})

    assert len(HSD_index) == 4
    assert len(HSE_index) == 0
    assert len(HSP_index) == 1

    ##########################################
    # ##   Compute protonation  at PH 10.0 ###
    out_file = os.path.join(tmp_path, '4N1M_amber_ph10.pqr')
    pdb2pqr.compute_pdb2pqr(os.path.join(PDB_4N1M),
                            out_file, ff=ff, ph=10.0)
    prot_coor = pdb_manip.Coor(out_file)
    assert prot_coor.num == 2547

    HSD_index = prot_coor.get_index_selection(
        {'res_name': ['HID'], 'name': ['CA']})
    HSE_index = prot_coor.get_index_selection(
        {'res_name': ['HIE'], 'name': ['CA']})
    HSP_index = prot_coor.get_index_selection(
        {'res_name': ['HIP'], 'name': ['CA']})

    assert len(HSD_index) == 5
    assert len(HSE_index) == 0
    assert len(HSP_index) == 0


def test_pdb2pqr_dna(tmp_path):
    """
    Check that pdb2pqr produce something with dna

    """

    ff = "CHARMM"

    pdb_manip.show_log()
    ##########################################
    # ##   Compute protonation  at PH 4.0  ###
    ##########################################
    out_file = os.path.join(tmp_path, '1D30_charmm_ph4.pqr')
    pdb2pqr.compute_pdb2pqr(os.path.join(PDB_1D30),
                            out_file, ff=ff, ph=4.0)
    prot_coor = pdb_manip.Coor(out_file)
    assert prot_coor.num == 486
