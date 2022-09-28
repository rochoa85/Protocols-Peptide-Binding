#!/usr/bin/python

"""
Prepare and run machine learning model with output energy variable
NOTE: The protocol requires of auxiliary programas and Unix system commands - Tested on Ubuntu 16.04

From publication "Open protocols for docking and MD-based scoring of peptide substrates"
Artificial Intelligence in the Life Sciences, 2022
Authors: Rodrigo Ochoa, Angel Santiago, Melissa Alegría-Arcos

Third-party tools required:

Scikit-Learn: https://scikit-learn.org/stable/
RDKit: https://github.com/rdkit/rdkit/releases
"""

import pickle
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from sklearn import datasets, ensemble
from sklearn.linear_model import LinearRegression
from sklearn.inspection import permutation_importance
from sklearn.metrics import r2_score
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import train_test_split
from scipy.stats import pearsonr
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit import SimDivFilters

#from classes_ML_training import *

# Classes
class TrainTestSplit:

    def __init__(self):
        pass

    @classmethod
    def max_chem_diversity(self, df_dataset, smiles_column = "smiles", test_set_size = None):
        """It splits the dataset into a training and a test set such that the test set contains a maximumally diverse set of compounds.
        To use this function, the input dataset has to contain a smiles_column.
        Help function for diversity_balanced and diversity_stratified functions.

        The chemical similarity is quantified using the ECFP4 Tanimoto coefficient.
        To obtain the test set, a maximally diverse subset of compounds is selected using the MaxMin algorithm in the RDKit.
        The remaining compounds are included in the training set.
        """
        random_seed = None
        if smiles_column not in list(df_dataset):
            print("Error: the column {} is not contained in the input dataframe. A column containing SMILES is required to compute Morgan fingerprints.".format(smiles_column))
            return

        if 'dataset' in df_dataset:
            df_dataset = df_dataset.drop(columns = ['dataset'])

        df_dataset.reset_index(inplace = True, drop = True)
        indexes = list(df_dataset.index)

        mols =  [Chem.MolFromSmiles(m) for m in df_dataset[smiles_column]]
        morgan_fps  = [AllChem.GetMorganFingerprint(x,2) for x in mols]
        nfps = len(morgan_fps)

        if test_set_size != None:
            n_to_pick = int(test_set_size)
        else:
            n_to_pick = int(df_dataset.shape[0]/5)

        if n_to_pick > nfps:
            print("Error: the specified test_set_size is larger than the total number of instances in the dataset. Only 20% of the dataset instances are included in the test set.")
            n_to_pick = int(df_dataset.shape[0]/5)

        if random_seed == None:
            random_seed = int(np.random.randint(n_to_pick, size=1))

        def distij(i,j,morgan_fps=morgan_fps):
            return (1-DataStructs.DiceSimilarity(morgan_fps[i],morgan_fps[j]))

        picker = SimDivFilters.MaxMinPicker()
        pickIndices = picker.LazyPick(distij, nfps, n_to_pick, seed=random_seed)
        test_indexes = list(pickIndices)
        test_indexes = list(np.take(indexes, test_indexes))
        train_indexes = list(set(df_dataset.index) ^ set(test_indexes))
        df_training_set = df_dataset.iloc[train_indexes]
        df_test_set = df_dataset.iloc[test_indexes]
        df_training_set.reset_index(inplace = True, drop = True)
        df_test_set.reset_index(inplace = True, drop = True)

        return(df_training_set, df_test_set)


# Input files
list_peptides=[x.strip() for x in open("total_peptides.txt")]
score_peptides=[x.strip() for x in open("scores_vina_total.txt")]

# Read the scores
score_dictionary={}
for s in score_peptides:
    info=s.split()
    score_dictionary[info[0]]={"score1":info[1]}

total_peptides=[]

# Read the object with all the descriptors
descriptors=[]
for i,pep in enumerate(list_peptides):
    # Load data (deserialize)
    with open('dict_objects/{}.pkl'.format(pep), 'rb') as handle:
        unserialized_data = pickle.load(handle)
        if i==0:
            for key in unserialized_data:
                if key != "cmpd_name" and key != "smiles":# and key != "prot_lig_lj_av_PL" and key != "prot_lig_lj_med_PL":
                    descriptors.append(key)

        unserialized_data["score"]=float(score_dictionary[unserialized_data["cmpd_name"]]["score1"])
        total_peptides.append(unserialized_data)

print(len(total_peptides))
df_peptides=pd.DataFrame(total_peptides)
print(df_peptides)

# Split in train and test using 75% and 25% respectively
df_peptides_train, df_peptides_test = TrainTestSplit.max_chem_diversity(df_peptides, smiles_column='smiles', test_set_size = df_peptides.shape[0]/4)
print(descriptors)

# List with peptide descriptors
peptide_descriptors=['MW','HA_count','RB_count','N_count','O_count','F_count','P_count','S_count','Cl_count','Br_count','I_count',
                     'HBD_count','HBA_count','2d_shape','2d_psa','net_charge','avg_hydro','isoelectric_p','sol_fail','syn_fail']

interaction_descriptors=[]
for d in descriptors:
    if d not in peptide_descriptors:
        interaction_descriptors.append(d)

mdfp_train = np.array(df_peptides_train[descriptors])
mdfp_test = np.array(df_peptides_test[descriptors])

# This is to use only MD descriptors instead of both options
#mdfp_train = np.array(df_peptides_train[interaction_descriptors])
#mdfp_test = np.array(df_peptides_test[interaction_descriptors])

response_train = np.array(df_peptides_train["score"])
response_test = np.array(df_peptides_test["score"])

# Gradiend Boosting regressor
params = {'n_estimators': 500,
          'max_depth': 4,
          'min_samples_split': 5,
          'learning_rate': 0.01,
          'loss': 'ls'}

reg = ensemble.GradientBoostingRegressor(**params)
reg.fit(mdfp_train, response_train)

print(reg.predict(mdfp_test))
mse = mean_squared_error(response_test, reg.predict(mdfp_test))
print("The mean squared error (MSE) on test set: {:.4f}".format(mse))
pr = pearsonr(response_test, reg.predict(mdfp_test))[0]
print("Pearson correlation is: {:.4f}".format(pr))
r2 = r2_score(response_test, reg.predict(mdfp_test))
print("R2 score is: {:.4f}".format(r2))

# Print model deviance for training and test sets
test_score = np.zeros((params['n_estimators'],), dtype=np.float64)
for i, y_pred in enumerate(reg.staged_predict(mdfp_test)):
    test_score[i] = reg.loss_(response_test, y_pred)

fig = plt.figure(figsize=(6, 6))
plt.subplot(1, 1, 1)
plt.title('Deviance')
plt.plot(np.arange(params['n_estimators']) + 1, reg.train_score_, 'b-',
         label='Training Set Deviance')
plt.plot(np.arange(params['n_estimators']) + 1, test_score, 'r-',
         label='Test Set Deviance')
plt.legend(loc='upper right')
plt.xlabel('Boosting Iterations')
plt.ylabel('Deviance')
fig.tight_layout()
plt.show()

# Plot feature_importance
feature_importance = reg.feature_importances_
sorted_idx = np.argsort(feature_importance)
pos = np.arange(sorted_idx.shape[0]) + .5
fig = plt.figure(figsize=(12, 12))
plt.subplot(1, 2, 1)
plt.barh(pos, feature_importance[sorted_idx], align='center')
plt.yticks(pos, np.array(descriptors)[sorted_idx],fontsize=6)
plt.title('Feature Importance (MDI)')

# Plot permutation_importance
result = permutation_importance(reg, mdfp_test, response_test, n_repeats=10,
                                random_state=42, n_jobs=2)
sorted_idx = result.importances_mean.argsort()
plt.subplot(1, 2, 2)
plt.yticks(fontsize=6)
plt.boxplot(result.importances[sorted_idx].T,
            vert=False, labels=np.array(descriptors)[sorted_idx])
plt.title("Permutation Importance (test set)")
fig.tight_layout()
plt.show()

# Test Linear Regression model
print("\nResults from the Linear Regression Model")
model = LinearRegression()
model.fit(mdfp_train, response_train)
predictions = model.predict(mdfp_test)
print("Pearson Lineal: ",pearsonr(response_test, predictions)[0])
print("MSE: ",mean_squared_error(response_test, predictions))
print("R2: ",r2_score(response_test, predictions))
