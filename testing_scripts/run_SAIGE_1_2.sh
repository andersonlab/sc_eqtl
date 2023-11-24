#!/usr/bin/env bash
# Perform the SAIGEQTL analysis of single cell expression from TI (test)

# Load modules and docker
module load ISG/singularity/3.9.0
saige_eqtl=/software/team152/bh18/singularity/singularity/saige.simg

# Define options for this test (will ultimately be inherited) and general options
category="Secretory"

phenotype__file="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/freeze_003/ti-cd_healthy-fr003_004/anderson_ti_freeze003_004-eqtl_processed.h5ad"
aggregate_on="category__machine"
general_file_dir="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/bradley_analysis/results/TI/SAIGE_runfiles"
genotype_pc__file=f"{general_file_dir}/14pcs.tsv"
genotype_id="Corrected_genotyping_ID"
sample_id="sanger_sample_id"
nperc=20
n_geno_pcs="5"
condition_col="disease_status"
condition="healthy"
covariates="age,sex"
expression_pca=True
n_geno_pcs=5
general_file_dir="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/bradley_analysis/results/TI/SAIGE_runfiles"
catdir=${general_file_dir}/${category}/${condition_col}/${condition}

# Construct covariance list by adding genotype PC strings to covs
for ((i=1; i<=n_geno_pcs; i++)); do
    covariates+=",PC$i"
done

# And expression PCs
if [[ "$expression_pca" == "True" ]]; then
    # Define file
    knee_file=${catdir}/knee.txt
    knee=$(<"$knee_file")
    # Generate the 'PC' string up to the numeric value
    covariates_cell=$(printf "xPC%d," $(seq "$knee") | sed 's/,$//')
fi

# Perform SAIGE with default options on a number of genes from the gene list file [This can be parallelised across all genes of a given condition] - GENOME WIDE
gene=ENSG00000160075
singularity exec -B /lustre -B /software $saige_eqtl step1_fitNULLGLMM_qtl.R &> ${catdir}/${gene}_step1.log  \
        --useSparseGRMtoFitNULL=FALSE  \
        --useGRMtoFitNULL=FALSE \
        --phenoFile=${catdir}/saige_filt_expr_input.txt	\
        --phenoCol=${gene}       \
        --covarColList="$covariates,$covariates_cell"    \
        --sampleCovarColList=${covariates}      \
        --sampleIDColinphenoFile=${genotype_id} \
        --traitType=count \
        --outputPrefix=${catdir}/${gene}_npc${n_geno_pcs} \
        --skipVarianceRatioEstimation=FALSE  \
        --isRemoveZerosinPheno=FALSE \
        --isCovariateOffset=FALSE  \
        --isCovariateTransform=TRUE  \
        --skipModelFitting=FALSE  \
        --tol=0.00001   \
        --plinkFile=${general_file_dir}/genotypes/plink_genotypes      \
        --IsOverwriteVarianceRatioFile=TRUE

# Perform the analysis genome-wide
step1prefix=${catdir}/${gene}_npc${n_geno_pcs}
step2prefix=${catdir}/${gene}__npc${n_geno_pcs}_gw.txt
singularity exec -B /lustre -B /software $saige_eqtl Rscript /usr/local/bin/step2_tests_qtl.R &> ${catdir}/${gene}_npc${n_geno_pcs}_log.txt \
        --bedFile=${general_file_dir}/genotypes/plink_genotypes.bed      \
        --bimFile=${general_file_dir}/genotypes/plink_genotypes.bim      \
        --famFile=${general_file_dir}/genotypes/plink_genotypes.fam      \
        --SAIGEOutputFile=${step2prefix}     \
        --minMAF=0.05 \
        --minMAC=20 \
        --LOCO=FALSE    \
        --GMMATmodelFile=${step1prefix}.rda     \
        --SPAcutoff=2 \
        --varianceRatioFile=${step1prefix}.varianceRatio.txt    \
        --markers_per_chunk=10000

# Run the ACAT test on the genome wide results
singularity exec -B /lustre -B /software $saige_eqtl Rscript /usr/local/bin/step3_gene_pvalue_qtl.R \
        --assocFile=${step2prefix}        \
        --geneName=$gene       \
        --genePval_outputFile=${catdir}/${gene}_gw_ACAT.txt




######## ~~~~~~~~~ SCRAP ~~~~~~~~~~~~ #############

# Run the analysis on this set (cis only) [This could be ran on chunks of genes ]
regionFile=${catdir}/ENSG00000187608_cis_region.txt
echo -e "1\t913497\t1114540" > ${regionFile}
step1prefix=${catdir}/ENSG00000187608
step2prefix=${catdir}/ENSG00000187608_cis.txt
singularity exec -B /lustre -B /software $saige_eqtl Rscript /usr/local/bin/step2_tests_qtl.R &> ${catdir}/ENSG00000187608_step2_cis_log.txt \
        --bedFile=${general_file_dir}/genotypes/plink_genotypes.bed      \
        --bimFile=${general_file_dir}/genotypes/plink_genotypes.bim      \
        --famFile=${general_file_dir}/genotypes/plink_genotypes.fam      \
        --SAIGEOutputFile=${step2prefix}     \
        --chrom=1       \
        --minMAF=0 \
        --minMAC=20 \
        --LOCO=FALSE    \
        --GMMATmodelFile=${step1prefix}.rda     \
        --SPAcutoff=2 \
        --varianceRatioFile=${step1prefix}.varianceRatio.txt    \
        --rangestoIncludeFile=${regionFile}     \
        --markers_per_chunk=10000

