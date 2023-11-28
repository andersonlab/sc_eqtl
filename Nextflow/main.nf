// Define the process to execute the get_samps_resolution_adata.py script
process PREPADATACATEGORY {
    label 'process_high_memory'

    input:
        file phenotype__file
        val aggregate_on
        path general_file_dir
        val genotype_id

    output:
        path("$inputParams.general_file_dir/$aggregate_on/levels.txt"), emit: level_file
        path("$inputParams.general_file_dir/samples.txt"), emit: sample_file

    script:
    """
        # Run the Python script with input parameters
        python bin/get_samps_resolution_adata.py \
            --phenotype__file $phenotype__file \
            --aggregate_on $aggregate_on \
            --general_file_dir $general_file_dir \
            --genotype_id $genotype_id
    """
}

// Define process to convert the genotypes and compute the PCA matrix using these samples. NOTE: These are different parameters from the origin nf-core/eqtl
process PLINK2BEDPCA {
    label 'process_medium'

    input:
        file__vcf
        general_file_dir

    output:
        path("$general_file_dir/genotypes/plink_genotypes.eigenvec"), emit: pca_file

    script:
    """
        # Prep dir
        if [ ! -d "$general_file_dir/genotypes" ]; then
            mkdir -p "$general_file_dir/genotypes"
            echo "Directory $general_file_dir/genotypes created."
        fi

        # Execute
        plink2 --vcf $file__vcf dosage=DS --make-bed --allow-extra-chr 0 --chr 1-22 XY --output-chr MT --snps-only --rm-dup exclude-all --hwe 0.0000001 --out $general_file_dir/genotypes/plink_genotypes

        # Then compute the PCA using the included samples only
        plink2 --bfile $general_file_dir/genotypes/plink_genotypes --keep $general_file_dir/keep_samples.txt --pca 20 --out $general_file_dir/genotypes/plink_genotypes

        # Also divide the plink file by chromosome for ease when running cis-only
        for chr_num in {1..22} X Y
        do
            # Use PLINK to extract data for each chromosome
            plink2 --bfile $general_file_dir/genotypes/plink_genotypes --chr "$chr_num" --make-bed --out $general_file_dir/genotypes/plink_genotypes_chr$chr_num
        done
    """
}


/*
process PROCESSADATAPERLEVEL {
    label 'process_high_memory'
    
    // Specify the maximum number of concurrent tasks in the array job
    maxForks 10

    input:
        levelOption from PREPADATACATEGORY.out
        phenotype__file
        aggregate_on
        genotype_pc_file
        genotype_id
        sample_id
        general_file_dir
        nperc
        condition_col
        condition
        covariates
        expression_pca
        scale_covariates

    output:
        path("$general_file_dir/$aggregate_on/$levelOption/saige_filt_expr_input.txt") emit: saige_filt_expr_input,
        path("$general_file_dir/$aggregate_on/$levelOption/test_genes.txt") emit: test_genes,
        path("$general_file_dir/$aggregate_on/$levelOption/knee.txt") emit: knee

    script:
    """
        python bin/prep_adata_saige.py 
            --level $levelOption
            --phenotype__file $phenotype__file 
            --aggregate_on $aggregate_on 
            --genotype_pc_file $genotype_pc_file 
            --genotype_id $genotype_id 
            --sample_id $sample_id 
            --general_file_dir $general_file_dir 
            --nperc $nperc 
            --condition_col $condition_col 
            --condition $condition 
            --covariates $covariates 
            --expression_pca $expression_pca 
            --scale_covariates $scale_covariates
    """
}


// Define the SAIGE run module
process RUNSAIGE {
    label 'process_low'

    // Specify the number of forks (10k)
    maxForks 1000

    // Input files for the array job
    input:
        levelOption from PREPADATACATEGORY.out
        phenotype__file
        aggregate_on
        general_file_dir
        n_geno_pcs
        covariates
        covariates_cell
        genotype_id
        sample_id
        annotation__file
        condition_col
        condition
        cis_only
        cis_window
        knee
        gene from path("$general_file_dir/$aggregate_on/$levelOption/test_genes.txt")

    // Output files
    output:
        path("${general_file_dir}/$aggregate_on/$levelOption/chr*_nPC_$n_geno_pcs.txt")

    // Define the Bash script to run for each array job
    script:
    """
        # Execute with the bash executable in an array (one job per gene within level)
        bash bin/RUNSAIGE_1_2.sh 
            -c levelOption 
            -p $phenotype__file 
            -a $aggregate_on 
            -d $general_file_dir 
            -w $n_geno_pcs 
            -e $covariates 
            -k $covariates_cell 
            -i $genotype_id 
            -s $sample_id 
            -m $annotation__file 
            -o $condition_col 
            -t $condition 
            -x $cis_only 
            -y $cis_window 
            -k $knee 
            -g gene
    """
}
*/

// Define workflow
workflow SAIGE_workflow{

    // Input parameters
    phenotype__file
    aggregate_on
    general_file_dir
    genotype_id
    file__vcf
    n_geno_pcs
    covariates
    covariates_cell
    expression_pca
    annotation__file
    cis_only
    cis_window


    main:

        // Connect processes to create the pipeline
        PREPADATACATEGORY.out --> PLINK2BEDPCA
        //PLINK2BEDPCA.out.sample_file --> PROCESSADATAPERLEVEL
        //PROCESSADATAPERLEVEL.out --> RUNSAIGE

}