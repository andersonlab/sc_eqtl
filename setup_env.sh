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
saige=/software/team152/bh18/singularity/singularity/saige.simg
singularity exec $saige step1_fitNULLGLMM_qtl.R --help # Check this prints a help message
singularity exec $saige step2_tests_qtl.R --help
singularity exec $saige step3_gene_pvalue_qtl.R --help
singularity exec $saige makeGroupFile.R  --help
# Try an example command
# NOTE: Help commands work but not any execution. This is because it is reliant on a SAIGE function which is not installed.
# install the singularity image from the original saige
mv saige.simg saige_eqtl.simg
saige_eqtl=/software/team152/bh18/singularity/singularity/saige_eqtl.simg
singularity build saige.simg docker://wzhou88/saige:0.45
saige=/software/team152/bh18/singularity/singularity/saige.simg

# Copy the saige package from the old container to the eqtl one
# Export the package from the source container
singularity shell -B /lustre -B /software $saige
cd /usr/local/lib/R/library/
tar -cvf /software/team152/bh18/temp/SAIGE.tar.gz SAIGE
singularity exec -B /software $saige_eqtl cp /software/team152/bh18/temp/SAIGE.tar.gz /software/team152/bh18/R/x86_64-pc-linux-gnu/3.6
singularity shell -B /lustre -B /software $saige_eqtl
cd /software/team152/bh18/R/x86_64-pc-linux-gnu/3.6 
tar -xvf SAIGE.tar.gz
# Also copy the whole lib of the saige container available to the eqtl one
singularity shell -B /lustre -B /software $saige
cp /usr/local/lib/R/lib/* /software/team152/bh18/temp/
singularity exec -B /software $saige_eqtl cp /software/team152/bh18/temp/libRlapack.so /software/team152/bh18/R/x86_64-pc-linux-gnu/3.6
# Make the saige eqtl containr writable
saige_write=/software/team152/bh18/singularity/saige_eqtl_sandbox.simg
singularity build -B /lustre -B /software --sandbox $saige_write $saige_eqtl
# Add the libRlapack file to the LD_LIBRARY_PATH of the writable container
singularity shell -B /lustre -B /software $saige_write
export LD_LIBRARY_PATH=/software/team152/bh18/temp/:$LD_LIBRARY_PATH


# Second try (all R packages)
# SAIGE package is present within the '/usr/local/lib/R/library/' directory of $saige
singularity exec -B /lustre -B /software $saige R -e 'lib <- .libPaths(); file.copy(lib, "/software/team152/bh18/temp/", recursive = TRUE)'
singularity exec  -B /lustre -B /software --bind /software/team152/bh18/temp/library:/usr/local/lib/R/site-library $saige_eqtl R

###### INvestigate why the libRlapack.so  is not available ######## <Check ChatGPT>

# NOTE: To work within singularity on the farm with acess to lustre, run: singularity run -B /lustre -B /software $saige
# This first step requires the original SAIGE, not SAIGEeQTL
export $PATH=/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/bradley_analysis/scripts/scRNAseq/SAIGE/qtl/extdata/:$PATH
cd /lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/bradley_analysis/scripts/scRNAseq/SAIGE/qtl/extdata/
singularity exec -B /lustre -B /software $saige_eqtl step1_fitNULLGLMM_qtl.R \
        --useSparseGRMtoFitNULL=FALSE  \
        --useGRMtoFitNULL=FALSE \
        --phenoFile=./input/seed_1_nfam_5_nindep_0_ncell_100_lambda_50_withCov_Poisson.txt	\
        --phenoCol=gene_1       \
        --covarColList=X1,X2,pf1,pf2    \
        --sampleCovarColList=X1,X2      \
        --sampleIDColinphenoFile=IND_ID \
        --traitType=count \
        --outputPrefix=./output/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_1 \
        --skipVarianceRatioEstimation=FALSE  \
        --isRemoveZerosinPheno=FALSE \
        --isCovariateOffset=FALSE  \
        --isCovariateTransform=TRUE  \
        --skipModelFitting=FALSE  \
        --tol=0.00001   \
        --plinkFile=./input/genotype_100markers_marker_plink_temp       \
        --IsOverwriteVarianceRatioFile=TRUE



#### STEP 2
# Downloaded example output from step1 from https://weizhou0.github.io/SAIGE-QTL-doc/docs/genomewide-eQTL.html into:
#  /lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/bradley_analysis/scripts/scRNAseq/SAIGE/qtl/extdata
cd /lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/bradley_analysis/scripts/scRNAseq/SAIGE/qtl/extdata
for i in {1..100}
  do
     step1prefix=../SAIGEQTL_step1_example_output/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_gene_${i}
     echo -e "gene_${i} ${step1prefix}.rda ${step1prefix}.varianceRatio.txt" >> ./input/step1_output_formultigenes.txt
  done

step2prefix=./output/nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_chr2
singularity exec -B /lustre -B /software $saige_eqtl Rscript step2_tests_qtl.R       \
      --bedFile=./input/genotype_100markers_2chr.bed      \
      --bimFile=./input/genotype_100markers_2chr.bim      \
      --famFile=./input/genotype_100markers_2chr.fam      \
      --SAIGEOutputFile=${step2prefix}     \
      --GMMATmodel_varianceRatio_multiTraits_File=./input/step1_output_formultigenes.txt      \
      --chrom=2       \
      --minMAF=0.05 \
      --LOCO=FALSE    \
      --SPAcutoff=2 \
      --markers_per_chunk=10000
