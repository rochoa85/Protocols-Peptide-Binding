#!/usr/bin/python

"""
Fragment docking protocol to predict peptide bound conformations
NOTE: The protocol requires of auxiliary programas and Unix system commands - Tested on Ubuntu 16.04

From publication "Open protocols for docking and MD-based scoring of peptide substrates"
Artificial Intelligence in the Life Sciences, 2022
Authors: Rodrigo Ochoa, Angel Santiago, Melissa Alegría-Arcos

Third-party tools required:

BioPython: https://biopython.org/wiki/Download - Ubuntu package: python-rdkit
RDKit: https://github.com/rdkit/rdkit/releases - Ubuntu package: python-biopython
AutoDock Vina: http://vina.scripps.edu/download.html - Ubuntu package: autodock-vina
Open Babel: https://sourceforge.net/projects/openbabel/ - Ubuntu package: openbabel

Modeller: https://salilab.org/modeller/download_installation.html
MGL Tools: http://mgltools.scripps.edu/downloads
PDB2PQR: https://apbs-pdb2pqr.readthedocs.io/en/latest/downloads.html
"""

########################################################################################
# Authorship
########################################################################################

__author__ = "Rodrigo Ochoa"
__license__ = "MIT"
__version__ = "1.0"
__email__ = "rodrigo.ochoa@udea.edu.co"

########################################################################################
# Modules to import
########################################################################################

import numpy as np
import unittest
import subprocess
import itertools
import os
import multiprocessing
import argparse
import yaml

# RDKit
from rdkit import Chem
from rdkit.Chem import AllChem

# BioPython
from Bio.PDB import *

# Modeller
import modeller
from modeller.automodel import *

########################################################################################
# General Functions
########################################################################################

def aminoacidSMILES(amino):
    """
    Obtain the SMILES representation of a particular amino acid

    Arguments:
    amino -- One-letter code of the amino acid

    Return:
    smiles -- 1D Chemical representation of the amino acid
    """

    # Dictionary with the SMILES per amino acid
    aminoacids = {'G':{'SMILES': 'NCC(=O)O'},
                  'A':{'SMILES': 'N[C@@]([H])(C)C(=O)O'},
                  'R':{'SMILES': 'N[C@@]([H])(CCCNC(=N)N)C(=O)O'},
                  'N': {'SMILES': 'N[C@@]([H])(CC(=O)N)C(=O)O'},
                  'D': {'SMILES': 'N[C@@]([H])(CC(=O)O)C(=O)O'},
                  'C': {'SMILES': 'N[C@@]([H])(CS)C(=O)O'},
                  'E': {'SMILES': 'N[C@@]([H])(CCC(=O)O)C(=O)O'},
                  'Q': {'SMILES': 'N[C@@]([H])(CCC(=O)N)C(=O)O'},
                  'H': {'SMILES': 'N[C@@]([H])(CC1=CN=C-N1)C(=O)O'},
                  'I': {'SMILES': 'N[C@@]([H])(C(CC)C)C(=O)O'},
                  'L': {'SMILES': 'N[C@@]([H])(CC(C)C)C(=O)O'},
                  'K': {'SMILES': 'N[C@@]([H])(CCCCN)C(=O)O'},
                  'M': {'SMILES': 'N[C@@]([H])(CCSC)C(=O)O'},
                  'F': {'SMILES': 'N[C@@]([H])(Cc1ccccc1)C(=O)O'},
                  'P': {'SMILES': 'N1[C@@]([H])(CCC1)C(=O)O'},
                  'S': {'SMILES': 'N[C@@]([H])(CO)C(=O)O'},
                  'T': {'SMILES': 'N[C@@]([H])(C(O)C)C(=O)O'},
                  'W': {'SMILES': 'N[C@@]([H])(CC(=CN2)C1=C2C=CC=C1)C(=O)O'},
                  'Y': {'SMILES': 'N[C@@]([H])(Cc1ccc(O)cc1)C(=O)O'},
                  'V': {'SMILES': 'N[C@@]([H])(C(C)C)C(=O)O'}}

    # Store the SMILES in a variable
    smiles=aminoacids[amino]['SMILES']
    return smiles

################################################################

def modelling_fragments(pepT,pepM,step,receptor,number_chains,num_model):
    """
    Function to increasily model fragment of the desired peptide

    Arguments:
    pepT -- Template peptide sequence
    pepM -- Peptide sequence to model
    step -- Step of the modelling process during the growing of the sequence
    receptor -- Structure used to dock the peptide
    number_chains -- Number of chains in the protein structure (max 3)
    num_model -- Number of the model that will be used as input

    Return
    PDB structure of the peptide after growing flanking amino acids
    """

    # Copy the PDB files that are going to be used in the modelling
    os.system('cp output/step{}/model{}_step{}.pdb .'.format(step,num_model,step))
    if number_chains==1: os.system("sed -i 's/ A  / B  /g' model{}_step{}.pdb".format(num_model,step))
    if number_chains==2: os.system("sed -i 's/ A  / C  /g' model{}_step{}.pdb".format(num_model,step))
    os.system("cat {}.pdb model{}_step{}.pdb > complex1_step{}.pdb".format(receptor,num_model,step,step))

    # Start the Modeller environment
    code = 'complex1_step{}'.format(step)
    e = modeller.environ()
    m = modeller.model(e, file=code)
    aln = modeller.alignment(e)
    aln.append_model(m, align_codes=code)
    aln.write(file=code+'.seq')

    # Edit the information of the sequence to store the sequences in the Modeller format
    infoSeq=[x.strip() for x in open(code+'.seq')]
    header=[]
    sequenceLine=''
    for info in infoSeq:
        if ">" not in info and ":" not in info:
            if info:
                sequenceLine+=info
        else:
            if info: header.append(info)

    # Store the sequences in variables according to Modeller format
    last_line=sequenceLine.split("/")
    if len(last_line)==2:
        sequenceTemp=last_line[0]+"/"+pepT+"*"
        sequenceMod=last_line[0]+"/"+pepM+"*"
    if len(last_line)==3:
        sequenceTemp=last_line[0]+"/"+last_line[1]+"/"+pepT+"*"
        sequenceMod=last_line[0]+"/"+last_line[1]+"/"+pepM+"*"

    seqTempList = [sequenceTemp[i:i+75] for i in range(0, len(sequenceTemp), 75)]
    seqModList = [sequenceMod[i:i+75] for i in range(0, len(sequenceMod), 75)]

    # Create the alignment file
    alignmentFile=open("alignment.ali","w")
    for h in header: alignmentFile.write(h+"\n")
    for s in seqTempList: alignmentFile.write(s+"\n")
    alignmentFile.write("\n>P1;complex1_step{}_fill\nsequence:::::::::\n".format(step))
    for s in seqModList: alignmentFile.write(s+"\n")
    alignmentFile.close()

    # Directories for input atom files
    e.io.atom_files_directory = ['.', '../atom_files']
    a = automodel(e, alnfile='alignment.ali', knowns='complex1_step{}'.format(step), sequence='complex1_step{}_fill'.format(step))
    a.starting_model= 1
    a.ending_model  = 1
    a.make()

    # Read the structure obtained from Modeller
    parser = PDBParser()
    structure = parser.get_structure('REF',"complex1_step{}_fill.B99990001.pdb".format(step))
    model = structure[0]
    ch_selected=""
    if number_chains==1: ch_selected="B"
    if number_chains==2: ch_selected="C"
    out_structure=model[ch_selected]
    numbers_to_change=[]

    # Store the numbers of the residues that should be changed
    for num_aa,residue in enumerate(out_structure):
        numbers_to_change.append(int(residue.get_full_id()[3][1]))

    # Save a temportal PDB file
    io = PDBIO()
    io.set_structure(out_structure)
    io.save("post-modelled.pdb")

    # Replace the atoms and residues
    for newN,oldN in enumerate(numbers_to_change):
        diffN=len(str(oldN))-len(str(newN+1))
        if newN+1<10:
            if diffN==0: os.system("sed -i 's/ {}   {} / A   {} /g' post-modelled.pdb".format(ch_selected,oldN,newN+1))
            if diffN==1: os.system("sed -i 's/ {}  {} / A   {} /g' post-modelled.pdb".format(ch_selected,oldN,newN+1))
            if diffN==2: os.system("sed -i 's/ {} {} / A   {} /g' post-modelled.pdb".format(ch_selected,oldN,newN+1))
            #if diffN==3: os.system("sed -i 's/ {}{} / A {}    /g' post-modelled.pdb".format(ch_selected,oldN,newN+1))
        else:
            if diffN==0: os.system("sed -i 's/ {}  {} / A  {} /g' post-modelled.pdb".format(ch_selected,oldN,newN+1))
            if diffN==1: os.system("sed -i 's/ {} {} / A  {} /g' post-modelled.pdb".format(ch_selected,oldN,newN+1))
            #if diffN==2: os.system("sed -i 's/ {}{} / A  {} /g' post-modelled.pdb".format(ch_selected,oldN,newN+1))


    # Delete temporal files and save result
    os.system("rm complex1_step{}* model{}_step{}* alignment.ali".format(step,num_model,step))
    os.system("mv post-modelled.pdb pepM{}_step{}.pdb".format(num_model,step+1))

##########################################################################################

def generateConformer(sequence):
    """
    Function to generate a basic conformer based on the sequence

    Arguments:
    sequence -- Peptide fragment sequence

    Return:
    Peptide structure in PDB format
    """

    # Generate SMILES
    connect_smiles='O'
    for res in sequence:
         connect_smiles=connect_smiles[:-1]
         smiles=aminoacidSMILES(res)
         connect_smiles=connect_smiles+smiles
    final_smiles=connect_smiles

    # Generate molecule from smiles
    mol = Chem.MolFromSmiles(final_smiles)
    mol.SetProp("_Name",sequence)
    print("Generating the basic conformer for peptide {}".format(sequence))

    # Generate the conformer with the UFF force field
    AllChem.EmbedMolecule(mol)
    AllChem.UFFOptimizeMolecule(mol)
    writer = AllChem.SDWriter("pepM_{}.sdf".format(sequence))
    writer.write(mol)

    # Convert to PDB file
    os.system("obabel -isdf pepM_{peptide}.sdf -opdb > pepM_{peptide}.pdb".format(peptide=sequence))
    os.system("rm pepM_{}.sdf".format(sequence))

################################################################

def protonation(molecule,pH):
    """
    Function to protonate the molecules based on a particular pH

    Arguments:
    molecule -- Molecule that is going to be parameterized
    pH -- Definition of the pH to protonate

    Return:
    File after doing the protonation (.pqr format)
    """

    os.system("./scripts/pdb2pqr-linux/pdb2pqr --with-ph={} --ph-calc-method=propka --drop-water --apbs-input --ff=amber --verbose --chain --summary {mol}.pdb {mol}.pqr".format(pH,mol=molecule))
    #os.system("pdb2pqr --with-ph={} --ph-calc-method=propka --drop-water --apbs-input --ff=amber --verbose --chain --summary {mol}.pdb {mol}.pqr".format(pH,mol=molecule))
    os.system("rm {mol}.in {mol}-input.p {mol}.summary {mol}.propka".format(mol=molecule))

################################################################

def generate_box_initial(sequence,center_x,center_y,center_z,pdb,init_threshold,dist_threshold):
    """
    Function to generate the initial box for the docking

    Arguments:
    sequence -- Sequence of the peptide fragment that will be docking
    center_x -- center in x of the defined box
    center_y -- center in y of the defined box
    center_z -- center in z of the defined box
    pdb -- Structure file

    Return:
    Configuration file with the coordinates
    size_x -- Initial size of the cubic box in x
    size_y -- Initial size of the cubic box in y
    size_z -- Initial size of the cubic box in z
    """

    # Peptide lenght
    number_amino=len(sequence)

    # Read the structure
    parser = PDBParser()
    structure = parser.get_structure('PEP', pdb)
    model = structure[0]

    # Get the difference in coordinates
    #dist_threshold=2.5
    diff=model["A"][1]["CA"].coord-model["A"][number_amino]["CA"].coord
    diffValue=np.sqrt(np.sum(diff * diff))
    size_x=abs(diffValue*dist_threshold)
    size_y=abs(diffValue*dist_threshold)
    size_z=abs(diffValue*dist_threshold)

    #init_threshold=5
    if size_x<init_threshold: size_x==init_threshold
    if size_y<init_threshold: size_y==init_threshold
    if size_z<init_threshold: size_z==init_threshold

    # Exhaustiveness. This will depend on the number of cores
    exhaustiveness=int(multiprocessing.cpu_count())

    # Open the config file to write the information
    config=open("config.txt","w")
    config.write("center_x={}\n".format(center_x))
    config.write("center_y={}\n".format(center_y))
    config.write("center_z={}\n".format(center_z))
    config.write("size_x={}\n".format(size_x))
    config.write("size_y={}\n".format(size_y))
    config.write("size_z={}\n".format(size_z))
    config.write("exhaustiveness={}".format(exhaustiveness))
    config.close()

    # Return the cubic size
    return size_x,size_y,size_z

################################################################

def generate_box(sequence,center_x,center_y,center_z,pdb,ref_size_x,ref_size_y,ref_size_z,step,num_model,dist_threshold):
    """
    Function to generate the subsequent boxes for the docking

    Arguments:
    sequence -- Sequence of the peptide fragment that will be docking
    center_x -- center in x of the defined box
    center_y -- center in y of the defined box
    center_z -- center in z of the defined box
    pdb -- Structure file
    ref_size_x -- size of the previous box in x
    ref_size_y -- size of the previous box in y
    ref_size_z -- size of the previous box in z
    step -- Number of the docking stem
    num_model -- Number of the model obtained

    Return:
    Configuration file with the coordinates
    """

    # Peptide lenght
    number_amino=len(sequence)

    # Read the pdb structure
    parser = PDBParser()
    structure = parser.get_structure('PEP', pdb)
    model = structure[0]

    # Get the difference in coordinates
    diff_x=model["A"][1]["CA"].coord[0]-model["A"][number_amino]["CA"].coord[0]
    diff_y=model["A"][1]["CA"].coord[1]-model["A"][number_amino]["CA"].coord[1]
    diff_z=model["A"][1]["CA"].coord[2]-model["A"][number_amino]["CA"].coord[2]

    # Assign the sizes based on the growing direction of the peptide
    #dist_threshold=2.5
    if abs(diff_x*dist_threshold)>ref_size_x:
        size_x=abs(diff_x*dist_threshold)
    else:
        size_x=ref_size_x

    if abs(diff_y*dist_threshold)>ref_size_y:
        size_y=abs(diff_y*dist_threshold)
    else:
        size_y=ref_size_y

    if abs(diff_z*dist_threshold)>ref_size_z:
        size_z=abs(diff_z*dist_threshold)
    else:
        size_z=ref_size_z

    # Exhaustiveness. This will depend on the number of cores
    exhaustiveness=int(int(multiprocessing.cpu_count())/3)

    # Create the configuration file
    config=open("config{}.txt".format(num_model),"w")
    config.write("center_x={}\n".format(center_x))
    config.write("center_y={}\n".format(center_y))
    config.write("center_z={}\n".format(center_z))
    config.write("size_x={}\n".format(size_x))
    config.write("size_y={}\n".format(size_y))
    config.write("size_z={}\n".format(size_z))
    config.write("exhaustiveness={}".format(exhaustiveness))
    config.close()

    # Return the cubic size
    return size_x,size_y,size_z

################################################################

def docking(receptor,peptide,sequence,stepNumber,num_model):
    """
    Function to execute the docking

    Arguments:
    receptor -- Protein structure
    peptide -- Fragment that will be docked
    sequence -- Sequence of the peptide
    stepNumber -- Number of the step during the fragment growing
    num_model -- Number of the model
    """

    # List of amino acids
    aminoacids={"A":"ALA","D":"ASP","E":"GLU","F":"PHE","H":"HIS","I":"ILE","K":"LYS","L":"LEU","M":"MET","G":"GLY",
                "N":"ASN","P":"PRO","Q":"GLN","R":"ARG","S":"SER","T":"THR","V":"VAL","W":"TRP","Y":"TYR","C":"CYS"}

    # Check if this is the frist step
    if stepNumber==0:
        # Convert the molecules in PDBQT format
        os.system("./scripts/pythonsh scripts/prepare_receptor4.py -r {rec}.pqr -o {rec}.pdbqt -C -U waters".format(rec=receptor))
        os.system("./scripts/pythonsh scripts/prepare_ligand4.py -l {pep}.pqr  -C -U '' -B -o {pep}.pdbqt".format(pep=peptide))

        # Run Autodock Vina- NOTE: Install using: sudo apt-get install autodock-vina
        print("Docking step number {} ...".format(stepNumber))
        os.system("vina --receptor {}.pdbqt --ligand {}.pdbqt --log score.log --out out.pdbqt --config config.txt".format(receptor,peptide))

        # Split and get the first model
        os.system("csplit out.pdbqt /MODEL/ {{*}}; mv xx01 model1_step{}.pdbqt".format(stepNumber))
        os.system("mv xx02 model2_step{}.pdbqt".format(stepNumber))
        os.system("mv xx03 model3_step{}.pdbqt".format(stepNumber))
        os.system("rm xx* score.log out.pdbqt")

    else:
        # Select the amino acids that will be assigned as rigid
        rigidAA=[]
        for i in range(1,len(sequence)+1):
            if i!=1 and i!=2 and i!=len(sequence)-1 and i!=len(sequence):
                rigidAA.append(i)

        # Assign the rigid atoms
        rigidAtoms=[]
        for rigid in rigidAA:
            if rigid <=9:
                os.system("grep {} {}.pqr | grep 'A   {}' | awk '{{print $2}}' > atom{}.temp".format(aminoacids[sequence[rigid-1]],peptide,rigid,num_model))
            else:
                os.system("grep {} {}.pqr | grep 'A  {}' | awk '{{print $2}}' > atom{}.temp".format(aminoacids[sequence[rigid-1]],peptide,rigid,num_model))
            num=[x.strip() for x in open("atom{}.temp".format(num_model))]
            rigidAtoms=rigidAtoms+num
            os.system("rm atom{}.temp".format(num_model))

        # Get the atoms that will be inactivated
        inactive=[]
        os.system("./scripts/pythonsh scripts/prepare_ligand4.py -s -l {pep}.pqr -C -U '' -B -o {pep}.pdbqt".format(pep=peptide))
        os.system("grep REMARK {}.pdbqt | grep '  A  ' | awk '{{print $6\"\t\"$8}}' | sed 's/_/ /g' | awk '{{print $2\"_\"$4}}' > active{}.temp".format(peptide,num_model))
        bonds=[x.strip() for x in open("active{}.temp".format(num_model))]
        for b in bonds:
            data=b.split("_")
            if data[0] in rigidAtoms or data[1] in rigidAtoms:
                inactive.append("{}_{}".format(str(int(data[0])-1),str(int(data[1])-1)))
        os.system("rm active{}.temp".format(num_model))

        commandInactive="_".join(inactive)

        # Prepare the ligand based on the inactivated atoms
        os.system("./scripts/pythonsh scripts/prepare_ligand4.py -s -l {}.pqr -C -U '' -B -o {}.pdbqt -I {}".format(peptide,peptide,commandInactive))

        # Run Autodock Vina
        print("Docking step number {} ...".format(stepNumber))
        os.system("vina --receptor {}.pdbqt --ligand {}.pdbqt --log score.log --out out{}.pdbqt --config config{}.txt --num_modes 9".format(receptor,peptide,num_model,num_model))


################################################################

def new_coordinates(base,model_pdbqt):
    """
    Function to annotate the coordinates from the pqr file

    Arguments:
    base -- Structure protonated
    model_pdbqt -- Model requiring the information

    Return:
    PDB file with the new coordinates
    """

    # Read the structure
    print("Processing protein {} ...".format(model_pdbqt))
    parser = PDBParser()
    structure = parser.get_structure('PEP', "{}.pqr".format(base))
    model = structure[0]
    counter=1
    batoms=["C","N","O","H","S"]

    # Iterate over the residues
    for residue in model["A"]:
        resName=residue.get_resname()
        resNumber=residue.get_full_id()[3][1]

        for atom in residue:
            idAtom = atom.get_id()
            # Check the spaces of the PDB file based on the residue number
            if resNumber<=9:
                if idAtom in batoms:
                    bash="grep {} {}.pdbqt | grep 'A   {}' | grep '  {}  ' | awk '{{print $7}}'".format(resName,model_pdbqt,str(resNumber),idAtom)
                    new_x = subprocess.check_output(['bash','-c', bash])
                    bash="grep {} {}.pdbqt | grep 'A   {}' | grep '  {}  ' | awk '{{print $8}}'".format(resName,model_pdbqt,str(resNumber),idAtom)
                    new_y = subprocess.check_output(['bash','-c', bash])
                    bash="grep {} {}.pdbqt | grep 'A   {}' | grep '  {}  ' | awk '{{print $9}}'".format(resName,model_pdbqt,str(resNumber),idAtom)
                    new_z = subprocess.check_output(['bash','-c', bash])
                    atom.set_coord([float(new_x),float(new_y),float(new_z)])
                else:
                    if len(idAtom)==4: idAtom=idAtom[-1]+idAtom[:3]

                    bash="grep {} {}.pdbqt | grep 'A   {}' | grep ' {} ' | awk '{{print $7}}'".format(resName,model_pdbqt,str(resNumber),idAtom)
                    new_x = subprocess.check_output(['bash','-c', bash])
                    bash="grep {} {}.pdbqt | grep 'A   {}' | grep ' {} ' | awk '{{print $8}}'".format(resName,model_pdbqt,str(resNumber),idAtom)
                    new_y = subprocess.check_output(['bash','-c', bash])
                    bash="grep {} {}.pdbqt | grep 'A   {}' | grep ' {} ' | awk '{{print $9}}'".format(resName,model_pdbqt,str(resNumber),idAtom)
                    new_z = subprocess.check_output(['bash','-c', bash])
                    atom.set_coord([float(new_x),float(new_y),float(new_z)])
            else:
                if idAtom in batoms:
                    bash="grep {} {}.pdbqt | grep 'A  {}' | grep '  {}  ' | awk '{{print $7}}'".format(resName,model_pdbqt,str(resNumber),idAtom)
                    new_x = subprocess.check_output(['bash','-c', bash])
                    bash="grep {} {}.pdbqt | grep 'A  {}' | grep '  {}  ' | awk '{{print $8}}'".format(resName,model_pdbqt,str(resNumber),idAtom)
                    new_y = subprocess.check_output(['bash','-c', bash])
                    bash="grep {} {}.pdbqt | grep 'A  {}' | grep '  {}  ' | awk '{{print $9}}'".format(resName,model_pdbqt,str(resNumber),idAtom)
                    new_z = subprocess.check_output(['bash','-c', bash])
                    atom.set_coord([float(new_x),float(new_y),float(new_z)])
                else:
                    if len(idAtom)==4: idAtom=idAtom[-1]+idAtom[:3]

                    bash="grep {} {}.pdbqt | grep 'A  {}' | grep ' {} ' | awk '{{print $7}}'".format(resName,model_pdbqt,str(resNumber),idAtom)
                    new_x = subprocess.check_output(['bash','-c', bash])
                    bash="grep {} {}.pdbqt | grep 'A  {}' | grep ' {} ' | awk '{{print $8}}'".format(resName,model_pdbqt,str(resNumber),idAtom)
                    new_y = subprocess.check_output(['bash','-c', bash])
                    bash="grep {} {}.pdbqt | grep 'A  {}' | grep ' {} ' | awk '{{print $9}}'".format(resName,model_pdbqt,str(resNumber),idAtom)
                    new_z = subprocess.check_output(['bash','-c', bash])
                    atom.set_coord([float(new_x),float(new_y),float(new_z)])

            counter+=1

    # Save the final PDB structure
    io = PDBIO()
    io.set_structure(structure)
    io.save('%s.pdb' %model_pdbqt)

########################################################################################
########################################################################################
########################################################################################
# Main execution
########################################################################################
########################################################################################
########################################################################################
if __name__ == '__main__':

    # Script arguments
    parser = argparse.ArgumentParser(description='Fragment docking protocol to predict peptide bound conformations')
    parser.add_argument('-i', dest='config_file', type=argparse.FileType(mode='r'), required=True,
                        help='File containing all the necessary parameters to run the protocol')

    #####################################################################################
    # Assignment of parameters
    #####################################################################################
    args = parser.parse_args()
    if args.config_file:
        data = yaml.load(args.config_file)
        delattr(args, 'config_file')
        arg_dict = args.__dict__
        for key, value in data.items():
            if isinstance(value, list):
                for v in value:
                    arg_dict[key].extend(v)
            else:
                arg_dict[key] = value
    else:
        print("A config file is necessary to run the protocol. Exiting ...")
        exit()

    # Check the arguments
    if args.pep_seq:
        ref_peptide=args.pep_seq
    else:
        print("The parameter 'pep_seq' is required for the analysis. Exiting ...")
        exit()
    if args.pep_frag:
        peptide=args.pep_frag
    else:
        print("The parameter 'pep_frag' is required for the analysis. Exiting ...")
        exit()
    if args.target:
        receptor=args.target
    else:
        print("The parameter 'target' is required for the analysis. Exiting ...")
        exit()
    if args.num_chains:
        number_chains=int(args.num_chains)
    else:
        print("The parameter 'num_chains' is required for the analysis. Exiting ...")
        exit()
    if args.pep_ph:
        pH=args.pep_ph
    else:
        print("The parameter 'pep_ph' is required for the analysis. Exiting ...")
        exit()
    if args.center_x:
        center_x=args.center_x
    else:
        print("The parameter 'center_x' is required for the analysis. Exiting ...")
        exit()
    if args.center_y:
        center_y=args.center_y
    else:
        print("The parameter 'center_y' is required for the analysis. Exiting ...")
        exit()
    if args.center_z:
        center_z=args.center_z
    else:
        print("The parameter 'center_z' is required for the analysis. Exiting ...")
        exit()
    try:
        if args.dist_box:
            dist_threshold=args.dist_box
    except:
        print("The distance box threshold has not been provided. In that case the software will use a default value of 2.5 for the physical distance per dimension")
        dist_threshold=2.5
    try:
        if args.init_box:
            init_threshold=args.init_box
    except:
        print("The initial box threshold has not been provided. In that case the software will use a default value of 5 for the initial physical distance per dimension")
        init_threshold=5

    ####################################################################################
    # Start functions
    ####################################################################################

    # Predict peptide conformer
    generateConformer(peptide)

    # Protonation
    protonation(receptor,pH)
    protonation("pepM_{}".format(peptide),pH)

    # Docking
    size_x,size_y,size_z=generate_box_initial(peptide,center_x,center_y,center_z,"pepM_{}.pdb".format(peptide),init_threshold,dist_threshold)
    docking(receptor,"pepM_{}".format(peptide),peptide,0,1)

    # Count the number of models
    bash="ls | grep '_step0.pdbqt' | grep 'model' | wc -l"
    num_models = subprocess.check_output(['bash','-c', bash])
    if num_models==0:
        print("Docking error ... exiting the protocol")
        exit()

    new_coordinates("pepM_{}".format(peptide),"model1_step0")
    new_coordinates("pepM_{}".format(peptide),"model2_step0")
    new_coordinates("pepM_{}".format(peptide),"model3_step0")

    # Move to a folder with step0
    os.system("mkdir output/step0; mv pepM_{}.* *_step0.pdbqt *_step0.pdb config.txt output/step0".format(peptide))

    # list of future fragments for growing the peptide
    list_fragments=[]
    ref_limit=ref_peptide.index(peptide)
    count=0
    limit=0
    while limit==0:
        count+=1
        low_limit=ref_limit-count
        if low_limit<0: low_limit=0
        up_limit=ref_limit+len(peptide)+count
        if up_limit>len(ref_peptide): up_limit=len(ref_peptide)
        list_fragments.append(ref_peptide[low_limit:up_limit])
        if low_limit==0 and up_limit==len(ref_peptide):limit=1

    # Print the list of fragments
    print(list_fragments)

    # Iterate over the fragments
    final_step=0
    pep_reference=peptide
    for i,frag in enumerate(list_fragments):
        # Get the index of the fragment in the sequence
        index_sub=frag.index(pep_reference)
        pepTemplate=""
        control_size=0
        for pos,aa in enumerate(frag):
            if pos<index_sub:
                pepTemplate=pepTemplate+"-"
            else:
                control_size+=1
                if control_size > len(pep_reference):
                    pepTemplate=pepTemplate+"-"
                else:
                    pepTemplate=pepTemplate+aa
        pep_reference=frag

        for num in range(1,4):
            # Model the fragment
            modelling_fragments(pepTemplate,frag,i,receptor,number_chains,num)
            protonation("pepM{}_step{}".format(num,i+1),pH)

            # Generate the box
            if i==0:
                if num==1:
                    size_x1,size_y1,size_z1=generate_box(frag,center_x,center_y,center_z,"pepM{}_step{}.pdb".format(num,i+1),size_x,size_y,size_z,i+1,num,dist_threshold)
                if num==2:
                    size_x2,size_y2,size_z2=generate_box(frag,center_x,center_y,center_z,"pepM{}_step{}.pdb".format(num,i+1),size_x,size_y,size_z,i+1,num,dist_threshold)
                if num==3:
                    size_x3,size_y3,size_z3=generate_box(frag,center_x,center_y,center_z,"pepM{}_step{}.pdb".format(num,i+1),size_x,size_y,size_z,i+1,num,dist_threshold)
            else:
                if num==1:
                    size_x1,size_y1,size_z1=generate_box(frag,center_x,center_y,center_z,"pepM{}_step{}.pdb".format(num,i+1),size_x1,size_y1,size_z1,i+1,num,dist_threshold)
                if num==2:
                    size_x2,size_y2,size_z2=generate_box(frag,center_x,center_y,center_z,"pepM{}_step{}.pdb".format(num,i+1),size_x2,size_y2,size_z2,i+1,num,dist_threshold)
                if num==3:
                    size_x3,size_y3,size_z3=generate_box(frag,center_x,center_y,center_z,"pepM{}_step{}.pdb".format(num,i+1),size_x3,size_y3,size_z3,i+1,num,dist_threshold)

        # Apply here multiprocessing
        pool=multiprocessing.Pool()
        for num in range(1,4):
            #num_models=docking(receptor,"pepM_step{}".format(i+1),frag,i+1)
            pool.apply_async(docking, args=(receptor,"pepM{}_step{}".format(num,i+1),frag,i+1,num,))
        pool.close()
        pool.join()

        # Check the pose with the RMSD closer to the preivous fragment from the top 10
        for num in range(1,4):
            os.system("csplit out{}.pdbqt /MODEL/ {{*}} -f xx{}".format(num,num))
            bash="ls xx{}* | wc -l | awk '{{print $1}}'".format(num)
            num_results = subprocess.check_output(['bash','-c', bash]).strip().decode("utf-8")
            rmsd_values=[]
            for n in range(1,int(num_results)):
                os.system("mv xx{}0{} xx{}0{}.pdbqt".format(num,n,num,n))
                new_coordinates("pepM1_step{}".format(i+1),"xx{}0{}".format(num,n))
                os.system("gmx rms -s output/step{}/model{}_step{}.pdb -f xx{}0{}.pdb -fit none -o output.xvg < rmsd_inp.txt".format(i,num,i,num,n))
                bash="tail -n 1 output.xvg | awk '{print $2}'"
                rmsd = subprocess.check_output(['bash','-c', bash]).strip().decode("utf-8")
                os.system("rm output.xvg")
                rmsd_values.append(float(rmsd))
            print(rmsd_values)
            ind_min=rmsd_values.index(min(rmsd_values))+1
            os.system("mv xx{}0{}.pdbqt model{}_step{}.pdbqt".format(num,ind_min,num,i+1))
            os.system("rm xx* score.log out{}.pdbqt".format(num))

            # Count number models
            bash="ls | grep '_step{}.pdbqt' | grep 'model{}' | wc -l".format(i+1,num)
            num_models = subprocess.check_output(['bash','-c', bash]).strip().decode("utf-8")
            if num_models==0:
               print("Docking error ... exiting the protocol")
               exit()

        # Update the coordinates and save the results
        new_coordinates("pepM1_step{}".format(i+1),"model1_step{}".format(i+1))
        new_coordinates("pepM2_step{}".format(i+1),"model2_step{}".format(i+1))
        new_coordinates("pepM3_step{}".format(i+1),"model3_step{}".format(i+1))

        os.system("mkdir output/step{num_step}; mv pepM1_step{num_step}.* pepM2_step{num_step}.* pepM3_step{num_step}.* *_step{num_step}.pdb *_step{num_step}.pdbqt *_step{num_step}.pqr config*.txt output/step{num_step}".format(num_step=str(i+1)))
        final_step=i+1

    # Save the final complexes
    for num in range(1,4):
        os.system('cp output/step{}/model{}_step{}.pdb .'.format(final_step,num,final_step))
        if number_chains==1: os.system("sed -i 's/ A  / B  /g' model{}_step{}.pdb".format(num,final_step))
        if number_chains==2: os.system("sed -i 's/ A  / C  /g' model{}_step{}.pdb".format(num,final_step))
        os.system("cat {}.pdb model{}_step{}.pdb > output/final_complex{}_{}.pdb".format(receptor,num,final_step,num,ref_peptide))
        os.system("rm model{}_step{}.pdb".format(num,final_step))
