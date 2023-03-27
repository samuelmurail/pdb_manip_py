#!/usr/bin/env python3
# coding: utf-8

#

#####################################
# #######      DOCKQ       ##########
#####################################


import os
import logging
import numpy as np

from os_command_py import os_command

# In case pdb2pqr is launched as main, relative import will failed
try:
    from . import pdb_manip
except ImportError:
    print("Relative import from . fails, use absolute import instead")
    import pdb_manip

# Autorship information
__author__ = "Samuel Murail, Pierre Tuffery"
__copyright__ = "Copyright 2020, RPBS"
__credits__ = ["Samuel Murail"]
__license__ = "GNU General Public License v2.0"
__version__ = "1.0.1"
__maintainer__ = "Samuel Murail"
__email__ = "samuel.murail@u-paris.fr"
__status__ = "Production"

# Logging
logger = pdb_manip.logger



def compute_dockQ(self, native_coor, rec_chain=None,
    lig_chain=None, native_rec_chain=None, native_lig_chain=None,
    back_atom=pdb_manip.BACK_ATOM):

    # Remove hydrogens and non protein atoms as well as altloc
    model_coor = self.select_part_dict(
        selec_dict={'name': pdb_manip.HEAVY_ATOM,
                    'res_name': pdb_manip.PROTEIN_RES,
                    'alter_loc': ['', 'A']})
    native_coor = native_coor.select_part_dict(
        selec_dict={'name': pdb_manip.HEAVY_ATOM,
                    'res_name': pdb_manip.PROTEIN_RES,
                    'alter_loc': ['', 'A']})

    # Remove incomplete backbone residue:
    self.remove_incomplete_backbone_residues(back_atom=back_atom)
    native_coor.remove_incomplete_backbone_residues(back_atom=back_atom)

    # Get shortest chain's sequence to identify peptide and receptor chains
    model_seq = model_coor.get_aa_seq()
    native_seq = native_coor.get_aa_seq()

    if lig_chain is None:
        lig_chain = [min(model_seq.items(), key=lambda x:len(x[1].replace('-', '')))[0]]
    logger.info(f'Model ligand chain : {" ".join(lig_chain)}')
    if rec_chain is None:
        rec_chain = [chain for chain in model_seq if chain not in lig_chain]
    logger.info(f'Model receptor chain : {" ".join(rec_chain)}')

    if native_lig_chain is None:
        native_lig_chain = [min(native_seq.items(), key=lambda x:len(x[1].replace('-', '')))[0]]
    logger.info(f'Native ligand chain : {" ".join(native_lig_chain)}')
    if native_rec_chain is None:
        native_rec_chain = [chain for chain in native_seq if chain not in native_lig_chain]
    logger.info(f'Native receptor chain : {" ".join(native_rec_chain)}')

    # Put lig chain at the end of the dict:
    model_coor = self.select_part_dict(
        selec_dict={'chain': rec_chain + lig_chain})
    native_coor = native_coor.select_part_dict(
        selec_dict={'chain': native_rec_chain + native_lig_chain})

    self.set_good_chain_order(lig_chain)
    native_coor.set_good_chain_order(native_lig_chain)

    # Align model on native structure using model back atoms:
    rmsd_prot, [align_rec_index, align_rec_native_index] = self.align_seq_coor_to(
        native_coor, chain_1=rec_chain, chain_2=native_rec_chain,
        back_names=back_atom)
    logger.info(f'Receptor RMSD: {rmsd_prot:.3f} A')

    ########################
    # Compute ligand RMSD: #
    ########################

    lrmsd, [align_lig_index, align_lig_native_index] = self.align_seq_coor_to(
        native_coor, chain_1=lig_chain, chain_2=native_lig_chain, align=False,
        back_names=back_atom)
    logger.info(f'Ligand   RMSD: {lrmsd:.3f} A')

    self.set_same_resid_in_common(native_coor,
        align_rec_index + align_lig_index,
        align_rec_native_index + align_lig_native_index)

    ## Delete non common atoms:
    model_del_index = self.get_index_selection({'res_num': [0]})
    self.del_atom_index(index_list=model_del_index)
    native_del_index = native_coor.get_index_selection({'res_num': [0]})
    native_coor.del_atom_index(index_list=native_del_index)
    logger.info(f'Delete atoms {len(model_del_index)} in Model and {len(native_del_index)} in Native')

    irmsd = self.interface_rmsd(native_coor, native_rec_chain, native_lig_chain)
    logger.info(f'Interface   RMSD: {irmsd:.3f} A')
    fnat, fnonnat = self.compute_native_contact(native_coor, rec_chain,
        lig_chain, native_rec_chain, native_lig_chain)
    logger.info(f'Fnat: {fnat:.3f}      Fnonnat: {fnonnat:.3f}')

    def scale_rms(rms, d):
        return(1. / (1 + (rms / d)**2))

    d1 = 8.5
    d2 = 1.5

    dockq = (fnat + scale_rms(lrmsd, d1) + scale_rms(irmsd, d2)) / 3

    logger.info(f'DockQ Score pdb_manip: {dockq:.3f} ')

    return(
        {'Fnat': fnat,
        'Fnonnat': fnonnat,
        'rRMS': rmsd_prot,
        'iRMS': irmsd,
        'LRMS': lrmsd,
        'DockQ': dockq})


def remove_incomplete_backbone_residues(self, back_atom=pdb_manip.BACK_ATOM):
    """ Remove residue with non all backbone atoms
    """

    no_alter_loc = self.select_part_dict(
            {'alter_loc': ['A', '']})

    first_atom_index = list(no_alter_loc.atom_dict.keys())[0]
    uniq_resid = no_alter_loc.atom_dict[first_atom_index]['uniq_resid']
    back_num = 0
    uniq_res_to_del = []
    index_del_list = []
    local_del_list = []

    
    for atom_num, atom in no_alter_loc.atom_dict.items():


        if uniq_resid != atom['uniq_resid']:
            
            if back_num != len(back_atom):
                uniq_res_to_del.append(uniq_resid)
                index_del_list += local_del_list
            uniq_resid = atom['uniq_resid']
            back_num = 0
            local_del_list = [atom_num]
        else:
            local_del_list.append(atom_num)

        if atom['name'] in back_atom:
            back_num += 1

    # print(f'Delete residues {uniq_res_to_del}')
    # print(f'Deletes atoms {index_del_list}')
    self.del_atom_index(index_del_list)

def set_good_chain_order(self, lig_chain):
    """ Change atom index in order of having ligand at the end.
    Update also uniq_resid as function of the new order
    """

    lig = self.select_part_dict(
        selec_dict={'chain': lig_chain}
    )
    lig_num = lig.num
    tot_num = self.num
    rec_num = tot_num - lig_num

    atom_uniq_res = -1

    old_res = list(self.atom_dict.values())[0]['uniq_resid']

    new_index = 0
    new_atom_dict = {}
    for atom_num, atom in sorted(self.atom_dict.items()):
        if atom['uniq_resid'] != old_res:
            atom_uniq_res += 1
            old_res = atom['uniq_resid']

        atom['uniq_resid'] = atom_uniq_res

        if atom['chain'] not in lig_chain:
            new_atom_dict[new_index] = atom
            new_index += 1


    for atom_num, atom in sorted(lig.atom_dict.items()): 
        if atom['uniq_resid'] != old_res:
            atom_uniq_res += 1
            old_res = atom['uniq_resid']

        atom['uniq_resid'] = atom_uniq_res
        new_atom_dict[new_index] = atom
        new_index += 1

    self.atom_dict = new_atom_dict
    return

def set_same_resid_in_common(self, native_coor,
    atom_index, native_atom_index):
    """ Set same resid for each corresponding residues:
    """

    assert len(atom_index) == len(native_atom_index), "Not same size of index"
    residue = np.unique(np.array(
        [self.atom_dict[index]['uniq_resid'] for index in atom_index]))

    native_residue = np.unique(np.array(
        [native_coor.atom_dict[index]['uniq_resid'] for index in native_atom_index]))

    if len(atom_index) != len(residue) * 4:
        print('issue with atom number of model')
    if len(native_atom_index) != len(native_residue) * 4:
        print('issue with atom number of native structure')

    assert len(residue) == len(native_residue), 'Issue with residue number in common'

    res_id = 1
    self.change_pdb_field(change_dict = {"res_num": 0})
    native_coor.change_pdb_field(change_dict = {"res_num": 0})

    rec_residue = []
    for i, j in zip(residue, native_residue):
        model_res_index = self.get_index_selection(
                        {'uniq_resid': [i]})
        native_res_index = native_coor.get_index_selection(
                        {'uniq_resid': [j]})

        self.change_index_pdb_field(
            index_list=model_res_index,
            change_dict={"res_num" : res_id})
        native_coor.change_index_pdb_field(
            index_list=native_res_index,
            change_dict={"res_num" : res_id})
        res_id += 1

    return

def interface_rmsd(self, native_coor, 
    native_rec_chain, native_lig_chain, cutoff=10.0):
    """
    """

    native_rec = native_coor.select_part_dict(
        selec_dict={'chain': native_rec_chain}
    )
    native_lig = native_coor.select_part_dict(
        selec_dict={'chain': native_lig_chain}
    )

    rec_native_index_interface = native_rec.get_index_dist_between(
        native_lig, cutoff_min=0, cutoff_max=cutoff)
    rec_interface_residue = native_rec.get_attribute_selection(
        attribute='res_num',
        index_list=rec_native_index_interface)

    lig_native_index_interface = native_lig.get_index_dist_between(
        native_rec, cutoff_min=0, cutoff_max=cutoff)
    lig_interface_residue = native_lig.get_attribute_selection(
        attribute='res_num',
        index_list=lig_native_index_interface)

    interface_residue = rec_interface_residue + lig_interface_residue

    native_inter_index = native_coor.get_index_selection(
        selec_dict={'res_num': interface_residue,
                    'name': pdb_manip.BACK_ATOM}
    )
    model_inter_index = self.get_index_selection(
        selec_dict={'res_num': interface_residue,
                    'name': pdb_manip.BACK_ATOM}
    )

    print(len(native_inter_index), len(model_inter_index))
    assert len(native_inter_index) == len(model_inter_index), \
        f'Issue with interface atom number {len(native_inter_index)} {len(model_inter_index)}'
    assert len(native_inter_index) !=0, \
        f'Issue with interface, no atom number {len(native_inter_index)} {len(model_inter_index)}'

    self.align_to(
        native_coor, index_list=[model_inter_index, native_inter_index],
        rot_kabsch=False)
    irmsd = self.compute_rmsd_to(
        native_coor, index_list=[model_inter_index, native_inter_index])
    return irmsd

def compute_native_contact(self, native_coor, rec_chain,
    lig_chain, native_rec_chain, native_lig_chain, cutoff=5.0):
    """
    """

    rec = self.select_part_dict(
        selec_dict={'chain': rec_chain}
    )
    lig = self.select_part_dict(
        selec_dict={'chain': lig_chain}
    )

    native_rec = native_coor.select_part_dict(
        selec_dict={'chain': native_rec_chain}
    )
    native_lig = native_coor.select_part_dict(
        selec_dict={'chain': native_lig_chain}
    )


    rec_index = native_rec.get_index_dist_between(
        native_lig, cutoff_min=0, cutoff_max=cutoff)
    rec_residue = native_rec.get_attribute_selection(
        attribute='res_num',
        index_list=rec_index)

    native_contact_list = []

    # Get native contac
    for residue in rec_residue:
        # Select only the residue
        native_rec_res_coor = native_rec.select_part_dict(
            selec_dict={'res_num': [residue]}
        )
        ref_lig_index = native_lig.get_index_dist_between(
            native_rec_res_coor, cutoff_min=0, cutoff_max=cutoff)

        inter_lig_res = native_lig.get_attribute_selection(
            attribute='res_num',
            index_list=ref_lig_index)

        native_contact_list += [[residue, lig_res] for lig_res in inter_lig_res]
    logger.info(f'Native contact number {len(native_contact_list)}')

    # Get Model contacts:
    rec_index = rec.get_index_dist_between(
        lig, cutoff_min=0, cutoff_max=cutoff)
    rec_residue = rec.get_attribute_selection(
        attribute='res_num',
        index_list=rec_index)

    model_contact_list = []

    for residue in rec_residue:
        # Select only the residue
        model_rec_res_coor = rec.select_part_dict(
            selec_dict={'res_num': [residue]}
        )
        ref_lig_index = lig.get_index_dist_between(
            model_rec_res_coor, cutoff_min=0, cutoff_max=cutoff)

        inter_lig_res = lig.get_attribute_selection(
            attribute='res_num',
            index_list=ref_lig_index)
        model_contact_list += [[residue, lig_res] for lig_res in inter_lig_res]

    logger.info(f'Model contact number {len(model_contact_list)}')

    nat_in_model = 0
    for native_contact in native_contact_list:
        if native_contact in model_contact_list:
            nat_in_model += 1

    if len(native_contact_list) != 0:
        fnat = nat_in_model/len(native_contact_list)
    else:
        fnat = 0.0

    nonnat_in_model = 0
    for model_contact in model_contact_list:
        if model_contact not in native_contact_list:
            nonnat_in_model += 1

    if len(model_contact_list) != 0:
        fnonnat = nonnat_in_model/len(model_contact_list)
    else:
        fnonnat = 0.0

    return(fnat, fnonnat)

def compute_pdockQ(self, rec_chain=None, lig_chain=None, cutoff=8.00):
    """ Compute pdockQ 
    inspired form:
    https://gitlab.com/ElofssonLab/FoldDock/-/blob/main/src/pdockq.py
    """

    model_seq = self.get_aa_seq()
    
    if lig_chain is None:
        lig_chain = [min(model_seq.items(), key=lambda x:len(x[1].replace('-', '')))[0]]
    logger.info(f'Model ligand chain : {" ".join(lig_chain)}')
    if rec_chain is None:
        rec_chain = [chain for chain in model_seq if chain not in lig_chain]
    logger.info(f'Model receptor chain : {" ".join(rec_chain)}')

    rec_index = self.get_index_selection(
        {'name': ['CB'],
         'chain': rec_chain})
    rec_index += self.get_index_selection(
        {'name': ['CA'],
         'res_name': ['GLY'],
         'chain': rec_chain})
    rec = self.select_from_index(rec_index)

    lig_index = self.get_index_selection(
        {'name': ['CB'],
         'chain': lig_chain})
    print(lig_index)

    lig_index += self.get_index_selection(
        {'name': ['CA'],
         'res_name': ['GLY'],
         'chain': lig_chain})
    lig = self.select_from_index(lig_index)

    rec_contact = rec.get_index_dist_between(
                lig, cutoff_max=cutoff,
                cutoff_min=0)
    lig_contact = lig.get_index_dist_between(
                rec, cutoff_max=cutoff,
                cutoff_min=0)

    rec_plddt = rec.get_attribute_selection(
            selec_dict={}, attribute='beta',
            index_list=rec_contact)
    lig_plddt = lig.get_attribute_selection(
            selec_dict={}, attribute='beta',
            index_list=lig_contact)

    contact_num = rec_contact.shape[0] + lig_contact.shape[0]
    if contact_num == 0:
        return 0

    # print(rec_plddt, lig_plddt)

    avg_plddt = np.average(rec_plddt + lig_plddt)

    x = avg_plddt*np.log10(contact_num)
    pdockq = 0.724 / (1 + np.exp(-0.052*(x-152.611)))+0.018

    return pdockq

pdb_manip.Coor.compute_dockQ = compute_dockQ
pdb_manip.Coor.set_same_resid_in_common = set_same_resid_in_common
pdb_manip.Coor.interface_rmsd = interface_rmsd
pdb_manip.Coor.compute_native_contact = compute_native_contact
pdb_manip.Coor.remove_incomplete_backbone_residues = remove_incomplete_backbone_residues
pdb_manip.Coor.set_good_chain_order = set_good_chain_order
pdb_manip.Coor.compute_pdockQ = compute_pdockQ
