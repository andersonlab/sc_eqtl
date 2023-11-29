##### Bradley November 2023
# Testing the performance of SAIGEQTL on healthy TI samples
# Gathering the files neccesary to perform testing

# Specify the general_file_dir
general_file_dir="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/bradley_analysis/results/TI/SAIGE_runfiles"
geno_dir=${general_file_dir}/genotypes

# Specify the input vcf, phenotype_file, genotype pc file (processed h5ad object). NOTE: I am just using a fixed number of PCs here, will ultimately want to iterate across this
file__vcf="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/tobi_qtl_analysis/data/genotypes/dna_genotypes/2023_Jun/tobi_impute/CCF_OTAR-plates_1_2_3-imputed-all_chr.vcf.gz"
phenotype__file="/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/freeze_003/ti-cd_healthy-fr003_004/anderson_ti_freeze003_004-eqtl_processed.h5ad"
genotype_id="Corrected_genotyping_ID"
sample_id="sanger_sample_id"
aggregate_on="label__machine"
nperc_expression=1
condition_col=""
condition=""
covariates="age,sex"
expression_pca="true"
scale_covariates=True

testoutfile=${geno_dir}/plink_genotypes_chr9.fam
if [ ! -f "$testoutfile" ]; then
    # Make a plink file of each vcf chromosome [NOTE: This requires numeric chromosomes, so have updated '--output-chr' flag to specify this and therefore differs from the plink file generates in the eqtl nextflow pipeline]
    # bcf conda env
    pgen_or_bed="dosage=DS --make-bed"
    plink2_filters='--allow-extra-chr 0 --chr 1-22 XY --output-chr MT --snps-only --rm-dup exclude-all'
    plink2 --vcf ${file__vcf} ${pgen_or_bed} ${plink2_filters} --hwe 0.0000001 --out ${geno_dir}/plink_genotypes

    # Compute 20 genotype PCs and save into the geno dir (ultimately, these will be inherited)
    # Also only include samples in the anndata file used in analysis (a bit hacky, but won't be a problem once in the pipeline) - keep_samples made from the filtered anndata object of both conditions
    plink2 --bfile ${geno_dir}/plink_genotypes --keep ${general_file_dir}/keep_samples.txt --pca 20 --out ${geno_dir}/plink_genotypes

    # Also divide the plink file by chromosome for ease when running cis-only
    for chr_num in {1..22} X Y
    do
        # Use PLINK to extract data for each chromosome
        plink2 --bfile ${geno_dir}/plink_genotypes --chr "$chr_num" --make-bed --out ${geno_dir}/plink_genotypes_chr${chr_num}
    done
fi


