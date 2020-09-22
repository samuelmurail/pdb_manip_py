#!/usr/bin/env python3
# coding: utf-8

#####################################
# #######     PDB2PQR      ##########
#####################################


import os
import logging

from os_command_py import os_command

# In case pdb2pqr is launched as main, relative import will failed
try:
    from . import pdb_manip
except ImportError:
    print("Relative import from . fails, use absolute import instead")
    import pdb_manip

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


PDB2PQR_MOD_DIRNAME = os.path.dirname(os.path.abspath(__file__))

# Check if Readthedoc is launched skip the program path searching
on_rtd = os.environ.get('READTHEDOCS') == 'True'
if on_rtd:
    print("pdb2pqr cannot be found")
    PDB2PQR_BIN = ""
else:
    # Add 'pdb2pqr_cli' in case it is installed with conda
    PDB2PQR_BIN = os_command.which('pdb2pqr.py', 'pdb2pqr_cli')
    # PDB2PQR_BIN = os_command.which('pdb2pqr_cli')

# Test folder path
PQR_LIB_DIR = os.path.dirname(os.path.abspath(__file__))
TEST_PATH = os.path.abspath(os.path.join(PQR_LIB_DIR, "test/input/"))


def compute_pdb2pqr(pdb_in, pdb_out, ff="CHARMM",
                    method='propka', ph=7.0, check_file_out=True):
    """
    Use pdb2pqr to define protonation state of each residue of a protein.

    :param pdb_in: path of input pdb file
    :type pdb_in: str

    :param pdb_out: path of output pdb file
    :type pdb_out: str

    :param ff: forcefield nomenclature for atom names
    :type ff: str, optional, default="CHARMM"

    :param method:  Method used to calculate ph values (propka or pdb2pka).
    :type method: str, optional, default="propka"

    :param ph: pH value to define AA residue
    :type ph: float, optional, default=7.0

    :param check_file_out: flag to check or not if file has already
        been created. If the file is present then the command break.
    :type check_file_out: bool, optional, default=True


    :Example:

    >>> pdb_manip.show_log()
    >>> TEST_OUT = str(getfixture('tmpdir'))
    >>> # Compute protonation with pdb2pqr:
    >>> compute_pdb2pqr(os.path.join(TEST_PATH,'4n1m.pdb'),
    ... os.path.join(TEST_OUT, '4n1m.pqr')) #doctest: +ELLIPSIS
    Succeed to read file ...4n1m.pdb ,  2530 atoms found
    Succeed to save file ...tmp_pdb2pqr.pdb
    pdb2pqr... --ff CHARMM --ffout CHARMM --chain --ph-calc-method=propka \
--with-ph=7.00 ...tmp_pdb2pqr.pdb ...4n1m.pqr
    0
    >>> prot_coor = pdb_manip.Coor()
    >>> prot_coor.read_pdb(os.path.join(TEST_OUT, '4n1m.pqr'), \
pqr_format = True) #doctest: +ELLIPSIS
    Succeed to read file ...4n1m.pqr ,  2549 atoms found
    >>> HSD_index = prot_coor.get_index_selection({'res_name' : ['HSD'],
    ... 'name':['CA']})
    >>> print(len(HSD_index))
    4
    >>> HSE_index = prot_coor.get_index_selection({'res_name' : ['HSE'],
    ... 'name':['CA']})
    >>> print(len(HSE_index))
    0
    >>> HSP_index = prot_coor.get_index_selection({'res_name' : ['HSP'],
    ... 'name':['CA']})
    >>> print(len(HSP_index))
    1

    .. note::
        Idealy I would need a pdb file with 3 different histidine protonation.
        I couldn't find one.

    """

    # print("Compute pdb2pqr on",pdb_in)

    # Check if output files exist and create directory:
    if check_file_out and os_command.check_file_and_create_path(pdb_out):
        logger.info("pdb2pqr not launched %s already exist" % pdb_out)
        return pdb_out

    out_folder = os_command.get_directory(pdb_out)
    # print("out_folder", out_folder)

    # WARING :
    # Many bugs are due to the REMARK field in pdb2pqr
    # The 2 following steps remove the REMARK field of the pdb

    tmp_coor = pdb_manip.Coor()
    tmp_coor.read_pdb(pdb_in)

    # Remove HETATM
    no_hetatm_pdb = tmp_coor.select_part_dict({'field': 'ATOM'})
    no_hetatm_pdb.write_pdb(os.path.join(out_folder + "/tmp_pdb2pqr.pdb"))

    cmd_pdb2pqr = os_command.Command(
        [PDB2PQR_BIN,
         "--ff", ff,
         "--ffout", ff,
         "--chain",
         "--ph-calc-method={}".format(method),
         "--with-ph={:.2f}".format(ph),
         os.path.join(out_folder + "/tmp_pdb2pqr.pdb"),
         pdb_out])

    cmd_pdb2pqr.display()
    out_data = cmd_pdb2pqr.run()
    os_command.delete_file(os.path.join(out_folder + "/tmp_pdb2pqr.pdb"))

    return out_data


if __name__ == "__main__":

    import doctest
    import shutil

    TEST_DIR = 'pdb_manip_test_out'
    TEST_OUT = os.path.join(TEST_DIR, 'pdb2pqr_test')

    def getfixture(*args):
        return TEST_OUT

    print("-Test pdb2pqr module:")
    print("pdb2pqr:\t", doctest.testmod())

    # Erase all test files
    shutil.rmtree(TEST_DIR, ignore_errors=True)
