/*
 * Skrypt przygotowany przez Julia Denkewicz (119644) na potrzeby projektu zaliczeniowego z Bioinformatyki roślin
 */

/*
 * Pipeline input parameters
 */

params.reads = "$baseDir/sample_data/reads/*_{R1,R2}.fq"
params.transcriptome_fa = "$baseDir/sample_data/refs/22.fa"
params.transcriptome_gtf = "$baseDir/sample_data/refs/22.gtf"
params.out_dir = "$baseDir/results"

log.info """"\
            (hisat2) RNA seq PIPELINE [HBR AND UTR SAMPLES]
            ===============================================================================
            transcriptome : ${params.transcriptome_fa}
            reads         : ${params.reads}
            results       : ${params.out_dir}
        """
        .stripIndent()


/*
 * PART I: Index
 * Create a FASTA genome index using hisat2 (2 CPU threads were used in the calculations)
 */
process index {
    cpus 2

    input:
    path transcriptome from params.transcriptome_fa

    output:
    path 'hisat_ref_genom' into indexes_map

    script:
    """
    mkdir hisat_ref_genom
    hisat2-build -p $task.cpus $transcriptome hisat_ref_genom/transcriptome.baseName
    """
}

Channel.fromFilePairs(params.reads, checkIfExists: true)
    .set{read_pairs}

/* PART II: Mapping
 * Mapping of reads to a reference genome (getting .bam files using samtools)
 */
process mapping {
    tag "$aligned_out"

    input:
    path index from indexes_map
    tuple aligned_out, path(reads) from read_pairs

    output:
    path ("*.bam") into map

    script:
    """
    INDEX=`find -L ./ -name "*.1.ht2" | sed 's/\\.1.ht2\$//'`
    hisat2 \\
            -x \$INDEX \\
            -1 ${reads[0]} \\
            -2 ${reads[1]} \\
            --summary-file ${aligned_out}.hisat2.summary.log \\
            --no-mixed \\
            | samtools view -bS > ${aligned_out}.bam
    """
}

/* PART III: Quantification
 * Using the feature counts program and the package from the Bioconductor (Rsubread) featureCounts library to count the readings
 */
process quantification {
    container = 'quay.io/biocontainers/subread:2.0.1–-hed695b0_0'

    input:
    path (bam_file) from map.collect()
    path annotation from params.transcriptome_gtf

    output:
    path ("*featureCounts.txt")
    path ("*featureCounts.txt.summary")
    path ("simple.counts.txt") into all_matrixCount

    script:
    """
    featureCounts \\
        -p \\
        -a $annotation \\
        -o output.featureCounts.txt \\
        *.bam
    cat output.featureCounts.txt | cut -f 1,7,8,9,10,11,12 > simple.counts.txt
    """
}


/* PART IV: DESeq2
 * Demonstration of diffrent expression using DESeq2 for data normalization (R script)
 */
process DESeq {
    container = 'quay.io/biocontainers/bioconductor-deseq2:1.38.0--r42hc247a5b_0'
    publishDir params.out_dir, mode: 'copy'
    
    input:
    path count_matrix from all_matrixCount

    output:
    path ("*.csv")

    script:
    """
    Rscript "$baseDir/script_deseq2.R" ${count_matrix}
    """
}
