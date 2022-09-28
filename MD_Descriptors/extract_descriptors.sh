#!/bin/bash

# Script to extract MD and peptide-based descriptors
# NOTE: The protocol requires of auxiliary programas and Unix system commands - Tested on Ubuntu 16.04

# From publication "Open protocols for docking and MD-based scoring of peptide substrates"
# Artifical Intelligence in the Life Sciences, 2022
# Authors: Rodrigo Ochoa, Angel Santiago, Melissa Alegría-Arcos
#
# Third-party tools required:
#
# Open Babel: https://sourceforge.net/projects/openbabel/ - Ubuntu package: openbabel
# Gromacs 5.1.4 (tested version): http://manual.gromacs.org/documentation/5.1.4/download.html

# Input variables
# NOTE: All the MD files should be named with the peptide sequence in order to be read by the script
peptide=$1
md_time=$2

# Creating the indexes required
gmx make_ndx -f ${peptide}.gro -o energy.ndx < indexing.txt

# Obtain the peptide, protein and full system in separate files
gmx trjconv -s ${peptide}.tpr -f ${peptide}.xtc -o peptide.pdb -n energy.ndx -b $md_time < choice_peptide.txt
gmx trjconv -s ${peptide}.tpr -f ${peptide}.xtc -o protein.pdb -n energy.ndx -b $md_time < choice_protein.txt
gmx trjconv -s ${peptide}.tpr -f ${peptide}.xtc -o system.pdb -n energy.ndx -b $md_time < choice_system.txt

# Convert peptide to sdf format
obabel -ipdb peptide.pdb -osdf > peptide.sdf

# Rerun trajectories to capture specific energy terms
gmx grompp -f md-RF.mdp -c ${peptide}.gro -r ${peptide}.gro -p topol.top -o reaction_field.tpr -n energy.ndx
gmx mdrun -s reaction_field.tpr -rerun ${peptide}.xtc -e ene_reaction_field.edr
gmx energy -f ene_reaction_field.edr -o ene_reaction_field.xvg < energy_terms_peptide.txt

# Call main python script to extract the descriptors
python extract_MDFP.py ${peptide}.xtc ene_reaction_field.xvg ${peptide}.gro peptide.sdf protein.pdb system.pdb topol.top $peptide
