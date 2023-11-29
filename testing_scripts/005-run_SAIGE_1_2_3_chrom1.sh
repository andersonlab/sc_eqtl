#!/usr/bin/env bash
# Perform the SAIGEQTL analysis of single cell expression from TI (test)
# bsub -o logs/saige_test-%J-output.log -e logs/saige_test-%J-error.log -q long -G team152 -n 4 -M 20000 -a "memlimit=True" -R "select[mem>20000] rusage[mem=20000] span[hosts=1]" -J "saige_test" < testing_scripts/run_SAIGE_1_2.sh 


# Load modules and docker
module load ISG/singularity/3.9.0
saige_eqtl=/software/team152/bh18/singularity/singularity/saige.simg

# Define options for this test (will ultimately be inherited) and general options
category="Secretory"
phenotype__file="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/freeze_003/ti-cd_healthy-fr003_004/anderson_ti_freeze003_004-eqtl_processed.h5ad"
aggregate_on="category__machine"
general_file_dir="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/bradley_analysis/results/TI/SAIGE_runfiles"
genotype_pc__file=${general_file_dir}/genotypes/plink_genotypes.eigenvec
genotype_id="Corrected_genotyping_ID"
sample_id="sanger_sample_id"
nperc=1
condition_col=""
condition=""
covariates="age_imputed,sex,Keras:predicted_celltype_probability"
covariates_cell="Keras:predicted_celltype_probability"
expression_pca="true"
annotation__file="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/tobi_qtl_analysis/repos/nf-hgi_eqtl/eqtl/assets/gene_counts_Ensembl_105_phenotype_metadata.annotation_file.txt"
cis_only=true
cis_window=1000000

echo "Prepping the directory variables"
# Construct the category directory path
if [ -n "$condition_col" ]; then
        catdir=${general_file_dir}/${category}/${condition_col}/${condition}
else
        catdir=${general_file_dir}/${category}
fi

# If submitting as an array job for chromosome 1 genes only:
gene_list=${catdir}/chr1_genes.txt
gene=$(head $gene_list -n ${LSB_JOBINDEX} | tail -n 1)
# bsub -o logs/saige_array_test-%J-%I-output.log -e logs/saige_array_test-%J-%I-error.log -q normal -G team152 -n 4 -M 10000 -a "memlimit=True" -R "select[mem>10000] rusage[mem=10000] span[hosts=1]" -J "saige_array_test[1-5]%200" < testing_scripts/005-run_SAIGE_1_2_3_chrom1.sh 

# Perform SAIGE with default options on a number of genes from the gene list file
for n_geno_pcs in 5 10 15 20; do
        echo "Working for genotype PC number"
        echo $n_geno_pcs
        echo "Prepping the covariates"
        
        for ((i=1; i<=n_geno_pcs; i++)); do
            covariates+=",PC$i"
        done

        # And expression PCs to the cell-level covariates
        if [[ "$expression_pca" == "True" ]]; then
        # Define file
        knee_file=${catdir}/knee.txt
        knee=$(<"$knee_file")
        # Generate the 'PC' string up to the numeric value
        covariates_cell="${covariates_cell},$(printf "xPC%d," $(seq "$knee") | sed 's/,$//')"
        fi

        # Fix covariate issue (replacing ':' with '_') in both the covariates and covatiates_cell
        covariates="${covariates//:/_}"
        covariates_cell="${covariates_cell//:/_}"
        # Combine
        covariates_sample_cell=$(echo "$covariates,$covariates_cell" | tr ',' '\n' | sort -u | tr '\n' ',' | sed 's/,$//')
        # Specify sample-level covariates
        covariates_sample=$(echo "$covariates" | tr ',' '\n' | sort | comm -23 - <(echo "$covariates_cell" | tr ',' '\n' | sort) | tr '\n' ',' | sed 's/,$//')


        echo "Estimating the variance"
        singularity exec -B /lustre -B /software $saige_eqtl step1_fitNULLGLMM_qtl.R \
                --useSparseGRMtoFitNULL=FALSE  \
                --useGRMtoFitNULL=FALSE \
                --phenoFile=${catdir}/saige_filt_expr_input.txt	\
                --phenoCol=${gene}       \
                --covarColList=${covariates_sample_cell}    \
                --sampleCovarColList=${covariates_sample}      \
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

        # Perform the analysis cis-only or genome-wide
        echo "Testing eQTLs"
        step1prefix=${catdir}/${gene}_npc${n_geno_pcs}
        if [ "$cis_only" = true ]; then
                step2prefix=${catdir}/${gene}__npc${n_geno_pcs}_cis
                # Find the coordinates/chromosome of the given gene
                gene_chr=$(awk -v search="$gene" '$1 == search {print $5}' "$annotation__file")
                gene_start=$(awk -v search="$gene" '$1 == search {print $2}' "$annotation__file")
                end_region=$(awk -v value="$gene_start" -v add="$cis_window" 'BEGIN {print value + add}')
                start_region=$(awk -v value="$gene_start" -v add="$cis_window" 'BEGIN {new_value = value - add; if (new_value < 0) new_value = 0; print new_value}')
                # Save this into a temporary file to perform the cis-only analysis
                echo -e "${gene_chr}\t${start_region}\t${end_region}" > ${step2prefix}_region_file.txt

                # Now perform the cis-only analysisq2
                echo "Performing the cis-only analysis"
                singularity exec -B /lustre -B /software $saige_eqtl Rscript /usr/local/bin/step2_tests_qtl.R \
                        --bedFile=${general_file_dir}/genotypes/plink_genotypes_chr${gene_chr}.bed      \
                        --bimFile=${general_file_dir}/genotypes/plink_genotypes_chr${gene_chr}.bim      \
                        --famFile=${general_file_dir}/genotypes/plink_genotypes_chr${gene_chr}.fam      \
                        --SAIGEOutputFile=${step2prefix}.txt     \
                        --chrom=${gene_chr}       \
                        --minMAF=0.05 \
                        --minMAC=20 \
                        --LOCO=FALSE    \
                        --GMMATmodelFile=${step1prefix}.rda     \
                        --SPAcutoff=2 \
                        --varianceRatioFile=${step1prefix}.varianceRatio.txt    \
                        --rangestoIncludeFile=${step2prefix}_region_file.txt     \
                        --markers_per_chunk=10000
                
                # Also perform the ACAT test for this gene
                singularity exec -B /lustre -B /software $saige_eqtl Rscript /usr/local/bin/step3_gene_pvalue_qtl.R \
                        --assocFile=${step2prefix}.txt      \
                        --geneName=$gene       \
                        --genePval_outputFile=${step2prefix}_ACAT.txt

                # Add gene name to the output
                awk -v new_val=${gene} 'BEGIN {OFS="\t"} {print $0, new_val}' "${step2prefix}.txt" > tmp && mv tmp ${step2prefix}.txt 

                # Remove the intermediate files (step 1 only)
                rm ${step1prefix}*
                echo "Finished analysis, removed intermediate files"
        else
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
        fi
done
