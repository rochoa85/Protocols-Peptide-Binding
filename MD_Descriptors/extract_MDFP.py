#!/usr/bin/python

"""
Extract MD and peptide descriptors with MDFP and PepFun packages
NOTE: The protocol requires of auxiliary programas and Unix system commands - Tested on Ubuntu 16.04

From publication "Protocol for fragment-growing docking and MD-based scoring of peptide substrates"
Molecular Informatics, 2021
Authors: Rodrigo Ochoa, Angel Santiago, Melissa Alegría-Arcos, Lucy Jimenez

- For the MDFP functions, please refer to the publication:
Esposito, C., Wang, S., Lange, U. E., Oellien, F., & Riniker, S. (2020).
Combining machine learning and molecular dynamics to predict P-glycoprotein substrates.
Journal of Chemical Information and Modeling, 60(10), 4730-4749.
The classes_MDFP file is adapted from the repository: https://github.com/rinikerlab/mdfptools/

- For the PepFun functions, please refer to the publication:
Ochoa, R., & Cossio, P. (2021).
PepFun: Open Source Protocols for Peptide-Related Computational Analysis.
Molecules, 26(6), 1664.

Third-party tools required:

ParmEd: https://github.com/ParmEd/ParmEd
RDKit: https://github.com/rdkit/rdkit/releases
MDTraj: https://www.mdtraj.org/1.9.5/index.html
"""

########################################################################################
# Authorship
########################################################################################

__author__ = "Rodrigo Ochoa"
__license__ = "MIT"
__version__ = "1.0"
__email__ = "rodrigo.ochoa@udea.edu.co"

########################################################################################
# Modules
########################################################################################

import sys
import parmed as pmd
import pickle
# Classes from the MDFP package
from classes_MDFP import *
# Classes from the PepFun package
from pepfun import *

# Input files
traj_file = sys.argv[1]
energy_file = sys.argv[2]
coord_file = sys.argv[3]
sdf_file = sys.argv[4]
pdb_file= sys.argv[5]
system_file= sys.argv[6]
top_file = sys.argv[7]
pep_name = sys.argv[8]

# Initialize list of dictionaries. one for each compound
Fingerprint = []
Dict_Fingerprint = {} # initialize dictionary
Dict_Fingerprint.update({'cmpd_name': pep_name})

# Read in molecule
mol = Chem.MolFromMolFile(sdf_file)
smiles = Chem.MolToSmiles(mol)
Dict_Fingerprint.update({"smiles": smiles})

# Compute 2D counts and properties and add them to the dictionary
Dict_Fingerprint.update(ComposerGMX.Properties2D_From_Mol(mol))


# Load pdb and identify indices of solute atoms
pdb = md.load(pdb_file)
topology = pdb.topology
solute_atoms = ComposerGMX.solute_solvent_split(topology,solute_residue_name = 1)[0] # If more than two chains, change number to 2, 3, or more

# Read in the trajectory but only for the solute
solute_traj = md.load(traj_file, top=pdb_file, atom_indices = solute_atoms)

# Update indexes
indexes_map={}
counter=0
for a in solute_traj.top.atoms:
    if a.residue.index not in indexes_map:
        indexes_map[a.residue.index]=counter
        counter+=1

for a in solute_traj.top.atoms:
    if a.residue.index in indexes_map:
        a.residue.index=indexes_map[a.residue.index]

# Compute SASA and Rgyr and add the terms to the output dictionary
Dict_Fingerprint.update(ComposerGMX.extract_rgyr(solute_traj))
Dict_Fingerprint.update(ComposerGMX.extract_sasa(solute_traj))

# Get energy components
energy_df = ComposerGMX.read_xvg(energy_file)
print(energy_df.head())
Dict_Fingerprint.update(ComposerGMX.extract_protein_ligang_energy_terms_single(df_ene=energy_df))

# Calculate charges
prm_obj = pmd.gromacs.GromacsTopologyFile(top_file)

# extract partial charges of the solute atoms
solute_charges = [i.charge for idx,i in enumerate(prm_obj.atoms) if idx in solute_atoms ]

# Use the dipole_moment function of mdtraj to compute the dipole moment components
dip_moments = md.dipole_moments(solute_traj, solute_charges)
# Obtain the dipole magnitude
dipole_magnitude = np.sqrt(np.square(dip_moments).sum(axis=1))

# Calculate mean, standard deviation, and median of the dipole moment x-, y-, z- components
dip_components_av = list(dip_moments.mean(axis=0))
dip_components_std = list(dip_moments.std(axis=0))
dip_components_med = list(np.median(dip_moments, axis=0))
# Calculate mean, standard deviation, and median of the dipole magnitude
dip_mom_stats = ComposerGMX.get_stats(dipole_magnitude)

# Update the MDFP output dictionary and print out the dataframe
dict_dip_mom = {'av_mu_x': dip_components_av[0], 'av_mu_y': dip_components_av[1], 'av_mu_z': dip_components_av[2], 'std_mu_x': dip_components_std[0], 'std_mu_y': dip_components_std[1], 'std_mu_z': dip_components_std[2], 'med_mu_x': dip_components_med[0], 'med_mu_y': dip_components_med[1], 'med_mu_z': dip_components_med[2], 'av_mu': dip_mom_stats[0], 'std_mu': dip_mom_stats[1], 'med_mu': dip_mom_stats[2]}
Dict_Fingerprint.update(dict_dip_mom)

# Calculate hydrogen bonds
t = md.load(traj_file, top=system_file)
hbonds = md.baker_hubbard(t, freq=0.5, periodic=False)
label = lambda hbond : '%s -- %s' % (t.topology.atom(hbond[0]), t.topology.atom(hbond[2]))
counter_hbonds=0
for hbond in hbonds:
    if t.topology.atom(hbond[0]).residue.chain.index==2 or t.topology.atom(hbond[2]).residue.chain.index==2:
        if t.topology.atom(hbond[0]).residue.chain.index != t.topology.atom(hbond[2]).residue.chain.index:
            print(label(hbond))
            counter_hbonds+=1
print(counter_hbonds)

# Pepfun
pep=peptide_sequence(pep_name)
pep.compute_peptide_charges()
net_charge=pep.netCharge
pep.calculate_properties_from_sequence()
avg_hydro=pep.avg_hydro
isoelectric_point=pep.isoelectric_point
pep.solubility_rules()
sol_fail=pep.solubility_rules_failed
pep.synthesis_rules()
syn_fail=pep.synthesis_rules_failed

# Update the MDFP output dictionary and print out the dataframe
dict_pep_prop = {'c_hbonds': counter_hbonds, 'net_charge': net_charge, 'avg_hydro': avg_hydro, 'isoelectric_p': isoelectric_point, 'sol_fail': sol_fail, 'syn_fail': syn_fail}
Dict_Fingerprint.update(dict_pep_prop)

print(Dict_Fingerprint)

# Store data (serialize)
with open('{}.pkl'.format(pep_name), 'wb') as handle:
    pickle.dump(Dict_Fingerprint, handle, protocol=pickle.HIGHEST_PROTOCOL)

# Load data (deserialize)
with open('{}.pkl'.format(pep_name), 'rb') as handle:
    unserialized_data = pickle.load(handle)

print(Dict_Fingerprint == unserialized_data)
