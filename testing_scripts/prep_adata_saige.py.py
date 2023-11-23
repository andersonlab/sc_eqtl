#!/usr/bin/env python

__author__ = 'Bradley Harris'
__date__ = '2023-11-23'
__version__ = '0.0.1'
####### Python script to prepare input files for SAIGE

# Load libraries
##################################################
#### Bradley August 2023
# Atlassing the rectum scRNAseq
# Using all of the expression data, not just limited to that which is genotyped.
##################################################

# Change dir
import os
cwd = os.getcwd()
print(cwd)
os.chdir("/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/bradley_analysis/results")
cwd = os.getcwd()
print(cwd)

# Load packages
import sys
sys.path.append('/software/team152/bh18/pip')
sys.path.append('/usr/local/')
print("System path")
print(sys.path)
import numpy as np
print("Numpy file")
print(np.__file__)
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib as mp
from matplotlib import pyplot as plt
from matplotlib.pyplot import rc_context
import kneed as kd
import scvi
import csv
import datetime
import seaborn as sns
import matplotlib.pyplot as plt
import math
import scipy.stats as st
import re
from scipy import stats
from sklearn.mixture import GaussianMixture
from scipy.stats import norm
from sympy import symbols, Eq, solve
from sklearn.preprocessing import StandardScaler
from scipy.optimize import fsolve
from scipy.optimize import brentq
import torch
from scib_metrics.benchmark import Benchmarker
print("Loaded libraries")

# Inherit options
phenotype__file = sys.argv[1] 
aggregate_on = sys.argv[2] 
genotype_pc__file = sys.argv[3] 
genotype_id = sys.argv[4] 
sample_id = sys.argv[5] 
general_file_dir = sys.argv[6] 
nperc = float(sys.argv[7]) 

# Testing in this instance, so define these
#phenotype__file = "/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/freeze_003/ti-cd_healthy-fr003_004/anderson_ti_freeze003_004-eqtl_processed.h5ad"
#aggregate_on = "category__machine"
#general_file_dir = "/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/bradley_analysis/results/TI/SAIGE_runfiles"
#genotype_pc__file = f"{general_file_dir}/14pcs.tsv"
#genotype_id = "Corrected_genotyping_ID"
#sample_id = "sanger_sample_id"
#nperc=10

# Load in the adata object
adata = ad.read_h5ad(phenotype__file)

# Load the genotype PCs
geno_pcs = pd.read_csv(genotype_pc__file, sep = "\t")
# Temp/ remove the prefix from the index
geno_pcs.index = geno_pcs.index.str.replace('label__machine-Tuft_cell-dMean_', '')
geno_pcs.index = geno_pcs.index.str.split('_').str[0]

# For each category of the aggregation column, we want to:
# 1. Filter for these cells
# 2. Identify the genes with expression in > n% of cells
# 3. Append this to the covariates (genotype PCs, Age, Sex encoded)
# 4. Save this in a new directory, including the aggregation category
cats = np.unique(adata.obs[aggregate_on])
for c in cats:
    print(f"~~~~~~~~~~~~Working on: {c}~~~~~~~~~~~~~~~~")
    print("Filtering anndata")
    temp = adata[adata.obs[aggregate_on] == c]
    # Filter for intersection with genotypes
    temp = temp[temp.obs[genotype_id].isin(geno_pcs.index)]
    # Identify genes expressed in > n% of cells
    print("Filtering lowly expressed genes")
    counts=temp.X
    count_per_column = np.array((counts > 1).sum(axis=0))[0]
    min_cells = math.ceil(temp.shape[0]*(nperc/100))
    keep_genes = np.where(count_per_column > counts.shape[0] * nperc/100)[0]
    # Subset counts for these genes
    counts = counts[:,keep_genes]
    # Bind this with the age, sex, sample and gt info after extracting gene names
    keep_genes_names = temp.var.index[keep_genes].values
    print("Adding age and sex covariates to genotype PCs")
    to_add = temp.obs[['sex', 'age']]
    # Recode sex
    to_add['sex'] = to_add['sex'].replace({'M': 1, 'F': 0})
    # Dummy variable for age
    to_add['Age_Tertile'] = pd.qcut(to_add['age'], q=[0, 1/3, 2/3, 1], labels=['Age_0', 'Age_1', 'Age_2'])
    to_add['Age_1'] = (to_add['Age_Tertile'] == 'Age_1').astype(int)
    to_add['Age_2'] = (to_add['Age_Tertile'] == 'Age_2').astype(int)
    to_add = to_add[["sex", "Age_1", "Age_2"]]
    # Bind this onto the counts
    counts = pd.DataFrame.sparse.from_spmatrix(counts)
    counts.index = temp.obs.index
    counts.columns = keep_genes_names
    counts = counts.merge(to_add, left_index=True, right_index=True)
    # Add the donor ID (genotyping ID so that we match the genotypes)
    counts[genotype_id] = temp.obs[genotype_id]
    # Save this as a dataframe in a directory specific to this resolution - make this if not already
    print("Saving")
    savedir=f"{general_file_dir}/{c}"
    if os.path.exists(savedir) == False:
        os.mkdir(savedir)

    # Save   
    counts.to_csv(f"{savedir}/saige_filt_expr_input.txt", sep = "\t", index=False)
    
    
    
    
    
