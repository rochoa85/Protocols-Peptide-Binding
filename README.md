# Protocol for fragment-growing docking and MD-based scoring of peptide substrates

### General Information and Third-Party Tools

Fragment docking protocol to predict peptide bound conformations
NOTE: The protocol requires of auxiliary programas and Unix system commands - Tested on Ubuntu 16.04

From publication "Protocol for fragment-growing docking and MD-based scoring of peptide substrates"
Molecular Informatics, 2021
Authors: Rodrigo Ochoa, Angel Santiago, Melissa Alegría-Arcos, Lucy Jimenez

**Third-party tools required:**

BioPython: https://biopython.org/wiki/Download - Ubuntu package: python-biopython
RDKit: https://github.com/rdkit/rdkit/releases - Ubuntu package: python-rdkit
AutoDock Vina: http://vina.scripps.edu/download.html - Ubuntu package: autodock-vina
Open Babel: https://sourceforge.net/projects/openbabel/ - Ubuntu package: openbabel

Modeller: https://salilab.org/modeller/download\_installation.html
MGL Tools: http://mgltools.scripps.edu/downloads
PDB2PQR: https://apbs-pdb2pqr.readthedocs.io/en/latest/downloads.html - Right now there is a local version

ParmEd: https://github.com/ParmEd/ParmEd
MDTraj: https://www.mdtraj.org/1.9.5/index.html
Scikit-Learn: https://scikit-learn.org/stable/

NOTE: These packages can be installed with CONDA

### How to access the fragment-docking script

An example of the script syntax is as follows::

`usage: fragment_docking.py [-h] -i CONFIG_FILE`

where the arguments are:

```
optional arguments:
  -h, --help      show this help message and exit
  -i CONFIG_FILE  File containing all the necessary parameters to run the
                  protocol

 ```

In addition of the script, the folder require of the target PDB structure file, a folder with a set of necessary scripts, and an output folder where the docking results step by step will be stored.

```
[target].pdb output scripts
```
An example to run the protocol script using the structure provided in the docker folder `/home/docking` is here:

`
python fragment_docking.py -i input.txt"
`

The method will start the modelling of the initial script, and the docking of the fragment after each growing step, until the peptide obtain the final desired size. The code can be modified to modifiy box sizes, as well as verify if the conformation of the growing ligand is according to previous findings of the biological system. An example of the output docked result, and the configuration each docking step is available in the output\_1NTV folder as follows:

`final_complex1_NFDNPVYRKT.pdb final_complex2_NFDNPVYRKT.pdb final_complex3_NFDNPVYRKT.pdb step0  step1  step2	step3  step4`


XXXXXXXXXXXXXX

10 nanoseconds
`bash extract_descriptors.sh LGPDESKQ 10000`
