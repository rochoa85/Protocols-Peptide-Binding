# Protocols for fragment-growing docking and MD-based scoring of peptide substrates

### General Information and Third-Party Tools

- From publication "Protocol for fragment-growing docking and MD-based scoring of peptide substrates"
- Molecular Informatics, 2021
- Authors: Rodrigo Ochoa, Angel Santiago, Melissa Alegría-Arcos, Lucy Jiménez

Here we present a set of protocols to dock peptides using a fragment-growing docking protocol for the *de novo* prediction of peptide conformations, a routine to capture descriptors from protein-peptide MD trajectories, and a script to predict observables such as average scoring values. As an application, a granzyme B protease was docked to a library of known peptide substrates and random sequences, and each complex was subjected to short MD simulations. Then a set of descriptors was calculated to generate a regression model able to predict with enough accuracy binding observables such as average scores from AutoDock Vina. The code to run the proposed protocols is available in this repository with some examples of execution.

**Third-party tools required:**

These were tested under an Ubuntu 20.04 operating system. **They can be installed using Conda to generate a virtual environment with all the requirements.**

- BioPython: https://biopython.org/wiki/Download - Ubuntu package: python-biopython
- RDKit: https://github.com/rdkit/rdkit/releases - Ubuntu package: python-rdkit
- AutoDock Vina: http://vina.scripps.edu/download.html - Ubuntu package: autodock-vina
- Open Babel: https://sourceforge.net/projects/openbabel/ - Ubuntu package: openbabel
- Modeller: https://salilab.org/modeller/download\_installation.html
- MGL Tools: http://mgltools.scripps.edu/downloads
- PDB2PQR: https://apbs-pdb2pqr.readthedocs.io/en/latest/downloads.html
- ParmEd: https://github.com/ParmEd/ParmEd
- MDTraj: https://www.mdtraj.org/1.9.5/index.html
- Scikit-Learn: https://scikit-learn.org/stable/

The project is split into three main protocols that are part from the publication. Some generalities and examples are provided in each section.

## 1. Fragment-docking protocol

This is a fragment-growing docking protocol for the *de novo* prediction of peptide conformations available in the folder `Fragment_Docking`. The protocol contain a subfolder called `scripts` with necessary software, and input files to run a basic example.

The script syntax is as follows::

`usage: fragment_docking.py [-h] -i CONFIG_FILE`

where the arguments are:

```
optional arguments:
  -h, --help      show this help message and exit
  -i CONFIG_FILE  File containing all the necessary parameters to run the
                  protocol
 ```
To run the script, we require a target PDB structure file, and an output folder where the docking results step by step will be stored. Based on the file `receptor.pdb`, and input file is provided to run the protocol:

```
pep_seq: TKSPYQEF
pep_frag: PYQ
target: receptor
num_chains: 1
pep_ph: 7.0
center_x: 18.139
center_y: 32.759
center_z: 73.242
dist_box: 2.5
init_box: 5
```

Where **pep_seq** is the peptide sequence that will be docked, **pep_frag** is the middle fragment the peptide will grow, **target** is the name of the PDB file with the protein target, **num_chains** is the number of chain of the protein, **pep_ph** is the desired pH to protonate the molecules, **center_[x,y,z]** are the box coordinates to start the docking, **dist_box** is the distance from the extreme atoms to generate the box and **init_box** is the minimum distance per dimension to generate the box. With that, the script can be called as `python fragment_docking.py -i input.txt`

The method will start the modelling and docking of the initial fragment, which grows by adding one amino acid at each flanking after a new step, until the peptide obtain the final desired size. The code can be modified to update box sizes, as well as verify if the conformation of the growing ligand is according to previous findings of the biological system.

## 2. MD-based and peptide descriptors

PENDING
`bash extract_descriptors.sh LGPDESKQ 10000`

## 3. Machine learning analysis

PENDING

## Support

In case the protocol is useful for other research projects and require some advice, please contact us to the email: rodrigo.ochoa@udea.edu.co
