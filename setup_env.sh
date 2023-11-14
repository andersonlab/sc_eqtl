# Bradley 14/11/23
# Setting up the conda/mamba environment to run SAIGETL
# This method is an application of SAIGE (https://github.com/weizhouUMICH/SAIGE), to run eQTL analyses on scRNAseq data
# It runs using code from: https://github.com/weizhou0/qtl/blob/7f54b239db26a37851208d2ca4961b63f8fae2d9/extdata/step1_fitNULLGLMM_qtl.R#L4
# And https://github.com/weizhou0/qtl/blob/7f54b239db26a37851208d2ca4961b63f8fae2d9/extdata/step2_tests_qtl.R#L4
# As described in https://github.com/annacuomo/SAIGE_QTL_analyses/blob/main/saige_qtl_runners/qtl_cmd.sh

# Clone the relevent environment (SAIGEQTL) [Change depending on where you keep scripts etc]
cd /lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/bradley_analysis/scripts/scRNAseq/SAIGE/qtl/conda_env
git clone git@github.com:weizhou0/qtl.git

# Install from environment yaml (following qtl/conda_env/createCondaEnvSAIGE_steps.txt)
mamba env create -f qtl/conda_env/environment-RSAIGE.yml
