# coding: utf-8

"""
Tests for pdb_manip functions
"""

import os
import pytest

from pdb_manip_py import pdb_manip


from .datafiles import PDB_1RXZ, PDB_1RXZ_MODEL, PDB_1JD4, PDB_5MN6

# Autorship information
__author__ = "Samuel Murail, Pierre Tuffery"
__copyright__ = "Copyright 2020, RPBS"
__credits__ = ["Samuel Murail"]
__license__ = "GNU General Public License v2.0"
__maintainer__ = "Samuel Murail"
__email__ = "samuel.murail@u-paris.fr"
__status__ = "Production"


def test_RMSD_Score_good(tmp_path):
    """
    *****************************************************************************
    *                                 TM-SCORE                                  *
    * A scoring function to assess the similarity of protein structures         *
    * Based on statistics:                                                      *
    *       0.0 < TM-score < 0.17, random structural similarity                 *
    *       0.5 < TM-score < 1.00, in about the same fold                       *
    * Reference: Yang Zhang and Jeffrey Skolnick, Proteins 2004 57: 702-710     *
    * For comments, please email to: zhng@umich.edu                             *
    *****************************************************************************
    
    Structure1: A693559     Length=  256
    Structure2: B693559     Length=  256 (by which all scores are normalized)
    Number of residues in common=  256
    RMSD of  the common residues=    0.702
    
    TM-score    = 0.9864  (d0= 5.92)
    MaxSub-score= 0.9631  (d0= 3.50)
    GDT-TS-score= 0.9658 %(d<1)=0.8633 %(d<2)=1.0000 %(d<4)=1.0000 %(d<8)=1.0000
    GDT-HA-score= 0.8477 %(d<0.5)=0.5273 %(d<1)=0.8633 %(d<2)=1.0000 %(d<4)=1.0000
    
     -------- rotation matrix to rotate Chain-1 to Chain-2 ------
     i          t(i)         u(i,1)         u(i,2)         u(i,3)
     1     67.6669778271   0.6958845067   0.0690236783   0.7148289901
     2     32.1967378763   0.4136668048  -0.8521810010  -0.3204174093
     3      0.9865965376   0.5870472961   0.5186745352  -0.6215723600
    
    Superposition in the TM-score: Length(d<5.0)=256  RMSD=  0.70
    (":" denotes the residue pairs of distance < 5.0 Angstrom)

    MIDVIMTGELLKTVTRAIVALVSEARIHFLEKGLHSRAVDPANVAMVIVDIPKDSFEVYNIDEEKTIGVDMDRI
    ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    MIDVIMTGELLKTVTRAIVALVSEARIHFLEKGLHSRAVDPANVAMVIVDIPKDSFEVYNIDEEKTIGVDMDRI
    12345678901234567890123456789012345678901234567890123456789012345678901234

    FDISKSISTKDLVELIVEDESTLKVKFGSVEYKVALIDPSAIRKEPRIPELELPAKIVMDAGEFKKAIAAADKI
    ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    FDISKSISTKDLVELIVEDESTLKVKFGSVEYKVALIDPSAIRKEPRIPELELPAKIVMDAGEFKKAIAAADKI
    56789012345678901234567890123456789012345678901234567890123456789012345678

    SDQVIFRSDKEGFRIEAKGDVDSIVFHMTETELIEFNGGEARSMFSVDYLKEFCKVAGSGDLLTIHLGTNYPVR
    ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    SDQVIFRSDKEGFRIEAKGDVDSIVFHMTETELIEFNGGEARSMFSVDYLKEFCKVAGSGDLLTIHLGTNYPVR
    90123456789012345678901234567890123456789012345678901234567890123456789012

    LVFELVGGRAKVEYILAPRIESEKSTQATLERWF
    ::::::::::::::::::::::::::::::::::
    LVFELVGGRAKVEYILAPRIESEKSTQATLERWF
    3456789012345678901234567890123456
    """

    model_coor = pdb_manip.Coor(PDB_1RXZ_MODEL)
    native_coor = pdb_manip.Coor(PDB_1RXZ)

    rmsd, tmscore, index = model_coor.align_seq_coor_to(
        native_coor, chain_1=['B', 'C'], chain_2=['A', 'B'],
        tmscore_flag=True)

    assert pytest.approx(rmsd, 0.001) == 0.702
    assert pytest.approx(tmscore, 0.001) == 0.9864
    assert len(index[0]) == 256

    # Test TM score with small peptide
    rmsd, tmscore, index = model_coor.align_seq_coor_to(
        native_coor, chain_1=['C'], chain_2=['B'],
        tmscore_flag=True)

    assert pytest.approx(rmsd, 0.001) == 0.432
    assert pytest.approx(tmscore, 0.01) == 0.86
    assert len(index[0]) == 11

def test_RMSD_Score_bad(tmp_path):

    model_coor = pdb_manip.Coor(PDB_1JD4)
    native_coor = pdb_manip.Coor(PDB_5MN6)

    rmsd, tmscore, index = model_coor.align_seq_coor_to(
        native_coor, chain_1=['A', 'B'], chain_2=['A', 'B'],
        tmscore_flag=True)

    assert pytest.approx(rmsd, 0.001) == 18.157
    assert pytest.approx(tmscore, 0.001) == 0.1608
    assert len(index[0]) == 187
