##### Bradley November 2023
# Testing the performance of SAIGEQTL on healthy TI samples
# Gathering the files neccesary to perform testing

# Specify the general_file_dir
general_file_dir="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/bradley_analysis/results/TI/SAIGE_runfiles"
geno_dir=${general_file_dir}/genotypes

# Specify the input vcf, phenotype_file, genotype pc file (processed h5ad object). NOTE: I am just using a fixed number of PCs here, will ultimately want to iterate across this
file__vcf="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/tobi_qtl_analysis/data/genotypes/dna_genotypes/2023_Jun/tobi_impute/CCF_OTAR-plates_1_2_3-imputed-all_chr.vcf.gz"
phenotype__file="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/freeze_003/ti-cd_healthy-fr003_004/anderson_ti_freeze003_004-eqtl_processed.h5ad"
genotype_pc__file="${general_file_dir}/14pcs.tsv"
genotype_id="Corrected_genotyping_ID"
sample_id="sanger_sample_id"
aggregate_on="category__machine"
nperc_expression="1"
n_geno_pcs="5"
condition_col="disease_status"
condition="healthy"
covariates="age,sex"
expression_pca=True
scale_covariates=True

# Make a plink file of each vcf chromosome [NOTE: This requires numeric chromosomes, so have updated '--output-chr' flag to specify this and therefore differs from the plink file generates in the eqtl nextflow pipeline]
# bcf conda env
pgen_or_bed="dosage=DS --make-bed"
plink2_filters='--allow-extra-chr 0 --chr 1-22 XY --output-chr MT --snps-only --rm-dup exclude-all'
plink2 --vcf ${file__vcf} ${pgen_or_bed} ${plink2_filters} --hwe 0.0000001 --out ${geno_dir}/plink_genotypes

# Run script to subset anndata object based on this aggregation column, identify genes to test, make the neccessary input files for SAIGE
# scvi_gpu3 environment
python prep_adata_saige.py \
    --phenotype__file $phenotype__file \
    --aggregate_on $aggregate_on \
    --genotype_pc__file $genotype_pc__file \
    --genotype_id $genotype_id \
    --sample_id $sample_id \
    --general_file_dir $general_file_dir \
    --nperc_expression $nperc_expression \
    --condition_col $condition_col \
    --condition $condition \
    --covariates $covariates \
    --expression_pca $expression_pca