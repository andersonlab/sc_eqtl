params.inputFile = file("input_TI.nf")

// Read input parameters from the JSON file
inputParams = file(params.inputFile).text | parseJson

// Define the process to execute the get_samps_resolution_adata.py script
process PREP_ADATA_CATEGORY {

  input:
    val phenotype__file
    val aggregate_on
    val general_file_dir
    val genotype_id

  output:
    path("${inputParams.general_file_dir}/{aggregate_on}/levels.txt"), emit: level_file
    path("${inputParams.general_file_dir}/samples.txt"), emit: sample_file

  script:
    """
        # Run the Python script with input parameters
        python get_samps_resolution_adata.py \
            --phenotype__file ${phenotype__file} \
            --aggregate_on ${aggregate_on} \
            --general_file_dir ${general_file_dir} \
            --genotype_id ${genotype_id}
    """
}

// Define process to convert the genotypes and compute the PCA matrix using these samples. NOTE: These are different parameters from the origin nf-core/eqtl
process PLINK_2_BED_PCA {

    input:
        val file__vcf
        val general_file_dir

    output:
        path("${general_file_dir}/genotypes/plink_genotypes.eigenvec"), emit: pca_file

  script:
    """
    # Prep dir
    if [ ! -d "${general_file_dir}/genotypes" ]; then
        mkdir -p "${general_file_dir}/genotypes"
        echo "Directory ${general_file_dir}/genotypes created."
    fi

    # Execute
    plink2 --vcf ${file__vcf} dosage=DS --make-bed --allow-extra-chr 0 --chr 1-22 XY --output-chr MT --snps-only --rm-dup exclude-all --hwe 0.0000001 --out ${general_file_dir}/genotypes/plink_genotypes

    # Then compute the PCA using the included samples only
    plink2 --bfile ${general_file_dir}/genotypes/plink_genotypes --keep ${general_file_dir}/keep_samples.txt --pca 20 --out ${general_file_dir}/genotypes/plink_genotypes

    # Also divide the plink file by chromosome for ease when running cis-only
    for chr_num in {1..22} X Y
    do
        # Use PLINK to extract data for each chromosome
        plink2 --bfile ${general_file_dir}/genotypes/plink_genotypes --chr "$chr_num" --make-bed --out ${general_file_dir}/genotypes/plink_genotypes_chr${chr_num}
    done

    """
}

process PROCESS_ADATA_PER_LEVEL {
    
    input:
        path levelsFile from path("${inputParams.general_file_dir}/{aggregate_on}/levels.txt")
        val phenotype__file
        val aggregate_on
        val general_file_dir
        val genotype_id
        val genotype_pc__file
        val sample_id
        val nperc
        val n_geno_pcs
        val covariates
        val expression_pca
        val scale_covariates

    output:
        path("${general_file_dir}/{aggregate_on}/saige_filt_expr_input.txt") emit: saige_filt_expr_input,
        path("${general_file_dir}/{aggregate_on}/test_genes.txt") emit: test_genes,
        path("${general_file_dir}/{aggregate_on}/knee.txt") emit: knee

    // Specify the maximum number of concurrent tasks in the array job
    maxForks 10

    script:
    """
        python prep_adata_saige.py \
            --phenotype__file ${phenotype__file} \
            --aggregate_on ${aggregate_on} \
            --genotype_pc__file ${genotype_pc__file} \
            --genotype_id ${genotype_id} \
            --sample_id ${sample_id} \
            --general_file_dir ${general_file_dir} \
            --nperc_expression ${nperc_expression} \
            --condition_col ${condition_col} \
            --condition ${condition} \
            --covariates ${covariates} \
            --expression_pca ${expression_pca} \
            
}