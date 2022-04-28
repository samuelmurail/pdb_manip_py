#!/usr/bin/env python3
# coding: utf-8

"""
Tests for PDB2PQR functions
"""

import numpy as np

from pdb_manip_py import pdb_manip

from .datafiles import PDB_1JD4

# Autorship information
__author__ = "Samuel Murail, Pierre Tuffery"
__copyright__ = "Copyright 2020, RPBS"
__credits__ = ["Samuel Murail"]
__license__ = "GNU General Public License v2.0"
__maintainer__ = "Samuel Murail"
__email__ = "samuel.murail@u-paris.fr"
__status__ = "Production"


def test_plot_3D(tmp_path):
    """
    Test plot_pseudo_3D with chains.

    """

    prot_coor = pdb_manip.Coor(PDB_1JD4)

    plot_chain = prot_coor.plot_pseudo_3D('chain')

    plot_properties = plot_chain.properties()

    assert plot_properties['color'].shape == (191, 4)
    assert len(plot_properties['segments']) == 191
    np.testing.assert_allclose(
        plot_properties['color'][0],
        [0.66666667, 0., 0.10666667, 1.], rtol=1e-06)
    np.testing.assert_allclose(
        plot_properties['color'][-1],
        [1., 0.33333333, 0.83333333, 1.], rtol=1e-06)

    plot_beta = prot_coor.plot_pseudo_3D('beta')

    plot_properties = plot_beta.properties()
    assert plot_properties['color'].shape == (191, 4)
    assert len(plot_properties['segments']) == 191
    np.testing.assert_allclose(
        plot_properties['color'][0],
        [0., 0.666667, 0.576288, 1.], rtol=1e-06)
    np.testing.assert_allclose(
        plot_properties['color'][-1],
        [0.333333, 0.650895, 1., 1.], rtol=1e-05)
