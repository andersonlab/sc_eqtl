# Bradley 14/11/23
# Setting up the conda/mamba environment to run SAIGETL
# This method is an application of SAIGE (https://github.com/weizhouUMICH/SAIGE), to run eQTL analyses on scRNAseq data
# It runs using code from: https://github.com/weizhou0/qtl/blob/7f54b239db26a37851208d2ca4961b63f8fae2d9/extdata/step1_fitNULLGLMM_qtl.R#L4
# And https://github.com/weizhou0/qtl/blob/7f54b239db26a37851208d2ca4961b63f8fae2d9/extdata/step2_tests_qtl.R#L4
# As described in https://github.com/annacuomo/SAIGE_QTL_analyses/blob/main/saige_qtl_runners/qtl_cmd.sh

# Try with singularity
# From (https://weizhou0.github.io/SAIGE-QTL-doc/docs/Installation_docker.html)
module load ISG/singularity/3.9.0
# Make symbolic link from HOMEDIR to where want to store cache
ln -s /software/team152/bh18/singularity/cache ~/.singularity/cache
cd /software/team152/bh18/singularity
singularity build saige.simg docker://wzhou88/saigeqtl:0.1.0
saige_eqtl=/software/team152/bh18/singularity/singularity/saige.simg
singularity exec $saige step1_fitNULLGLMM_qtl.R --help # Check this prints a help message
singularity exec $saige step2_tests_qtl.R --help
singularity exec $saige step3_gene_pvalue_qtl.R --help
singularity exec $saige makeGroupFile.R  --help

#####################################################
########## Try the example commands #################
#####################################################
# ~~~~~~~~~~~~~~~~~~~~ Step 1 ~~~~~~~~~~~~~~~~~~~~~~~
# Need to adjust the .fam file 
cd /lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/bradley_analysis/scripts/scRNAseq/SAIGE/qtl/extdata/
for f in ./input/genotype_100markers_2chr*
do
  echo $f
  newfile="./input/bh_$(basename "$f")"
  echo $newfile
  cp "$f" "$newfile"
done
awk -v OFS='\t' '{ $2 = substr($2, 2); print }' ./input/bh_genotype_100markers_2chr.fam > tmp && mv tmp ./input/bh_genotype_100markers_2chr.fam
singularity exec -B /lustre -B /software $saige_eqtl step1_fitNULLGLMM_qtl.R &> step1.log  \
        --useSparseGRMtoFitNULL=FALSE  \
        --useGRMtoFitNULL=FALSE \
        --phenoFile=./input/seed_1_nfam_5_nindep_0_ncell_100_lambda_50_withCov_Poisson.txt	\
        --phenoCol=y       \
        --covarColList=X1,X2,PC1,PC2    \
        --sampleCovarColList=PC1,PC2      \
        --sampleIDColinphenoFile=IND_ID \
        --traitType=count \
        --outputPrefix=./output/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_1 \
        --skipVarianceRatioEstimation=FALSE  \
        --isRemoveZerosinPheno=FALSE \
        --isCovariateOffset=FALSE  \
        --isCovariateTransform=TRUE  \
        --skipModelFitting=FALSE  \
        --tol=0.00001   \
        --plinkFile=./input/bh_genotype_100markers_2chr      \
        --IsOverwriteVarianceRatioFile=TRUE

# ~~~~~~~~~~~~~~~~~~~~ Step 2 ~~~~~~~~~~~~~~~~~~~~~~~
# Make the file that encodes the cis region
regionFile=./input/gene_1_cis_region.txt
echo -e "2\t300001\t610001" > ${regionFile}

##step1 output prefix and step 2 output prefix
step1prefix=./output/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_1
step2prefix=./output/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_1_cis

singularity exec -B /lustre -B /software $saige_eqtl Rscript /usr/local/bin/step2_tests_qtl.R &> step2.log \
        --bedFile=./input/bh_genotype_100markers_2chr.bed      \
        --bimFile=./input/bh_genotype_100markers_2chr.bim      \
        --famFile=./input/bh_genotype_100markers_2chr.fam      \
        --SAIGEOutputFile=${step2prefix}     \
        --chrom=2       \
        --minMAF=0 \
        --minMAC=20 \
        --LOCO=FALSE    \
        --GMMATmodelFile=${step1prefix}.rda     \
        --SPAcutoff=2 \
        --varianceRatioFile=${step1prefix}.varianceRatio.txt    \
        --rangestoIncludeFile=${regionFile}     \
        --markers_per_chunk=10000

# ~~~~~~~~~~~~~~~~~~~~ Step 2 - Genome wide ~~~~~~~~~~~~~~~~~~~~~~~
# This requires the use of the downloaded results from: https://weizhou0.github.io/SAIGE-QTL-doc/docs/genomewide-eQTL.html into 'SAIGEQTL_step1_example_output'
# Create the mapping file
rm ./input/step1_output_formultigenes.txt
for i in {1..100}
  do
     step1prefix=./SAIGEQTL_step1_example_output/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_${i}
     echo -e "gene_${i} ${step1prefix}.rda ${step1prefix}.varianceRatio.txt" >> ./input/step1_output_formultigenes.txt
  done

#run step 2 for all 100 genes 
step2prefix=./output/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_chr2
singularity exec -B /lustre -B /software $saige_eqtl Rscript /usr/local/bin/step2_tests_qtl.R       \
    --bedFile=./input/bh_genotype_100markers_2chr.bed      \
    --bimFile=./input/bh_genotype_100markers_2chr.bim      \
    --famFile=./input/bh_genotype_100markers_2chr.fam      \
    --SAIGEOutputFile=${step2prefix}     \
    --GMMATmodel_varianceRatio_multiTraits_File=./input/step1_output_formultigenes.txt      \
    --chrom=2       \
    --minMAF=0.05 \
    --LOCO=FALSE    \
    --SPAcutoff=2 \
    --markers_per_chunk=10000

# ~~~~~~~~~~~~~~~~~~~~ Step 3 ~~~~~~~~~~~~~~~~~~~~~~~
# Obtain gene-level p-values using the ACAT test
# https://weizhou0.github.io/SAIGE-QTL-doc/docs/gene_step3.html
# Run this on a single association (this just seems to be subsetting for the smallest p-value per gene)
singularity exec -B /lustre -B /software $saige_eqtl Rscript /usr/local/bin/step3_gene_pvalue_qtl.R \
        --assocFile=./output/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_1_cis        \
        --geneName=gene_1       \
        --genePval_outputFile=./output/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_1_cis_genePval