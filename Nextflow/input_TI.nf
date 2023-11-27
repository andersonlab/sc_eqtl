{
  "n_geno_pcs_values": [5, 10, 15, 20], // These will be iterated over and the one with the maximum number of cis-eGenes detected will be kept for downstream analysis
  "general_file_dir": "/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/bradley_analysis/results/TI/SAIGE_runfiles_nextflow", // work/results directory
  "file__vcf": "/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/tobi_qtl_analysis/data/genotypes/dna_genotypes/2023_Jun/tobi_impute/CCF_OTAR-plates_1_2_3-imputed-all_chr.vcf.gz",
  "phenotype__file": "/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/freeze_003/ti-cd_healthy-fr003_004/anderson_ti_freeze003_004-eqtl_processed.h5ad",
  "genotype_id": "Corrected_genotyping_ID", // Column in anndata object with the genotyping ID (must match vcf)
  "sample_id": "sanger_sample_id", // Non-genotype ID
  "aggregate_on": "category__machine", // Anndata object column on which cells are to be grouped by for the analyss
  "nperc_expression": "1", // The % of cells which must express a gene for it to be included in the analysis at the given resolution
  "condition_col": "", // Column of anndata object to filter on. If left black, all cells within a given resolution will be utilised
  "condition": "", // Value of condition_col to be used for the analysis
  "covariates": "age_imputed,sex", // Covariates used in the eQTL testing, these MUST be complete for all samples and match exactly one column in anndata object
  "expression_pca": true, // Include expression PCs in the eQTL analysis? If so, the number of PCs up to the 'knee' is used
  "scale_covariates": False, // Scale covariates for the analysis? two-level, categorical traits (e.g sex) will be binarised
  "saige_eqtl": "/software/team152/bh18/singularity/singularity/saige.simg" // Location of docker container for saige-eQTL (see: https://weizhou0.github.io/SAIGE-QTL-doc/docs/Installation_docker.html)
}