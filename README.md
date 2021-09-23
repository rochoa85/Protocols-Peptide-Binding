# Protocols for fragment-growing docking and MD-based scoring of peptide substrates

### General Information and Third-Party Tools

- From publication "Protocols for fragment-growing docking and MD-based scoring of peptide substrates"
- Molecular Informatics, 2021
- Authors: Rodrigo Ochoa, Angel Santiago, Melissa Alegría-Arcos, Lucy Jiménez

Here we present a set of protocols to dock peptides using a fragment-growing docking protocol for the *de novo* prediction of peptide conformations, a routine to capture descriptors from protein-peptide MD trajectories, and a script to predict observables such as average scoring values. As an application, a granzyme B protease was docked to a library of known peptide substrates and random sequences, and each complex was subjected to short MD simulations. Then a set of descriptors was calculated to generate a regression model able to predict with enough accuracy binding observables such as average scores from AutoDock Vina. The code to run the proposed protocols is available in this repository with some examples of execution.

**Third-party tools required:**

These were tested under an Ubuntu 20.04 operating system. **They can be installed using Conda to generate a virtual environment with all the requirements.**

- BioPython: https://biopython.org/wiki/Download - Ubuntu package: python-biopython
- RDKit: https://github.com/rdkit/rdkit/releases - Ubuntu package: python-rdkit
- AutoDock Vina: http://vina.scripps.edu/download.html - Ubuntu package: autodock-vina
- Open Babel: https://sourceforge.net/projects/openbabel/ - Ubuntu package: openbabel
- Modeller: https://salilab.org/modeller/download_installation.html
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

Where **pep_seq** is the peptide sequence that will be docked, **pep_frag** is the middle fragment the peptide will grow, **target** is the name of the PDB file with the protein target, **num_chains** is the number of chain of the protein, **pep_ph** is the desired pH to protonate the molecules, **center_[x,y,z]** are the box coordinates to start the docking, **dist_box** is the distance from the extreme atoms to generate the box and **init_box** is the minimum distance per dimension to generate the box.

With that input file, the script can be called as `python fragment_docking.py -i input.txt`

The method will start the modelling and docking of the initial fragment, which grows by adding one amino acid at each flanking after a new step, until the peptide obtain the final desired size. The code can be modified to update box sizes, as well as verify if the conformation of the growing ligand is according to previous findings of the biological system.

## 2. MD-based and peptide descriptors

Based on MD simulations of protein-peptide complexes with the Gromacs package, it is possible to extract a set of descriptors derived from the trajectories, as well as features from the peptide based on physico-chemical properties.

The MDFP tools libraries are implemented to extract a set of descriptors from the MD trajectories(https://github.com/rinikerlab/mdfptools/). Each descriptor is split into 3 positions on the vector, which include the average, median and standard deviation value of the calculated property among the MD frames. To achieve this, the calculated trajectory is rerun with Gromacs to add new energy terms in the outputs per frame. After that, multiple descriptors are captured. These include the Coulomb and Lennard-Jones energy contributions between the peptide, the receptor and the water molecules. Other properties are the SASA and radius of gyration, the charges calculated with the ParmED module, the dipole moments and evolution of hydrogen bonds with the MDtraj module, and bioinformatics properties of the peptide using the PepFun package. A total of 70 descriptors per complex are calculated.

The command to run the script is: `bash extract_descriptors.sh LGPDESKQ 10000`

In this scenario, it is required to have in the same folder the MD files with extensions **xtc, tpr and gro** obtained from the protein-peptide simulations, naming each file with the peptide sequence. Examples of these input files are located in the folder `example_MD_files` to reproduce the pipeline.

After running the scripts per protein-peptide complex, a vector is stored as a pickled object, which can be read later by machine learning models (see next section).

## 3. Machine learning analysis

Previous to run the ML model, it is required to obtain a response variable to train the model. In our example, an average docking score obtained from the MD frames was calculated per protein-peptide complex. The score used is the one from AutoDock Vina (*Vina*), but it can be energies predicted from more exhaustive simulations such as MM/PBSA and MM/GBSA calculations.

With the descriptors and defined response, a script is provided to train and run a regression model using two methods: a gradient booster regressor and a linear regression model. For both cases, the vectors are available in the `dict_objects` folder through **pickled** files obtained from the previous protocol. The script also requires a `total_peptides.txt` file with all the sequences, and a `scores_vina_total.txt` with the response variable, in this case the average *Vina* scores calculated from the MD simulations. **The paths and file names can be modified directly in the script.**

The command to run the script is: `python prediction_scores.py`

Two plots are generated, one with the gradient booster deviance of both the training and test sets, and a second plot with the feature and permutation importance. In addition, the R2 and MSE (Mean Square Errors) metrics are provided to assess the quality of the models.

**NOTE: All the three protocols can be combined to generate a pipeline of docking, sampling and scoring of peptide substrates for a protein target of interest.**

## Support

In case the protocol is useful for other research projects and require some advice, please contact us to the email: rodrigo.ochoa@udea.edu.co
