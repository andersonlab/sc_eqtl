#!/usr/bin/env python

__author__ = 'Bradley Harris'
__date__ = '2023-11-27'
__version__ = '0.0.1'

# Load packages
import anndata as ad
import scanpy as sc
import argparse
import numpy as np
import os

# Inherit options
parser = argparse.ArgumentParser(
        description="""
            Prepping files for the SAIGEQTL analysis
            """
    )

parser.add_argument(
        '-p', '--phenotype__file',
        action='store',
        dest='phenotype__file',
        required=True,
        help=''
    )

parser.add_argument(
        '-a', '--aggregate_on',
        action='store',
        dest='aggregate_on',
        required=True,
        help=''
    )

parser.add_argument(
        '-o', '--general_file_dir',
        action='store',
        dest='general_file_dir',
        required=True,
        help=''
    )

parser.add_argument(
        '-id', '--genotype_id',
        action='store',
        dest='genotype_id',
        required=True,
        help=''
    )

# Load the phenotype file
# phenotype__file="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/freeze_003/ti-cd_healthy-fr003_004/anderson_ti_freeze003_004-eqtl_processed.h5ad"
# general_file_dir="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/bradley_analysis/results/TI/SAIGE_nf_test"
# aggregate_on="category__machine"
adata=ad.read_h5ad(phenotype__file)

# Extract the category levels
cats = np.unique(adata.obs[aggregate_on])

# Save this in the general file dir, making this directory uf not already present
general_file_dir_cat=f"{general_file_dir}/{aggregate_on}"
if os.path.exists(general_file_dir_cat) == False:
    os.makedirs(general_file_dir_cat, exist_ok=True)
    
with open(f"{general_file_dir_cat}/levels.txt", 'w') as file:
        for item in cats:
            file.write(f"{item}\n")
            
# Save a list of samples
# Insert any code here if samples have missing covariates etc
with open(f"{general_file_dir}/samples.txt", 'w') as file:
        for item in np.unique(adata.obs[genotype_id]):
            file.write(f"{item}\n")