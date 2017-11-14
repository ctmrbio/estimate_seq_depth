#!/usr/bin/env nextflow
// vim: syntax=groovy expandtab
/****************************************
 * CTMR metagenome simulation and 
 * sequencing depth estimation pipeline 
 * ------------------------------------
 * Copyright (c) Authors 2017 
 * Authors:
 *  Fredrik Boulund <fredrik.boulund@ki.se>
 *  Luisa Hugerth <luisa.warchavchik.hugerth@ki.se>
 ****************************************/

Channel
    .fromFilePairs(params.input_reads)
    .ifEmpty{ exit 1, "Found no reference metagenome to subsample, did you specify --input_reads? I got: '${params.input_reads}'"}
    .set {
        input_reads
    }


process simulate_metagenome {
    tag "$sample_type"
    publishDir "${params.outdir}/simulated_metagenomes", mode: 'move'

    input:
    set pair_id, file(reads) from input_reads

    output:
    file "${sample_type}_10000000reads_{1,2,3}rep_{1,2}.fastq.gz"
    file "${sample_type}_1000000reads_{1,2,3}rep_{1,2}.fastq.gz"
    file "${sample_type}_100000reads_{1,2,3}rep_{1,2}.fastq.gz"
    file "${sample_type}_10000reads_{1,2,3}rep_{1,2}.fastq.gz"
    file "${sample_type}_1000reads_{1,2,3}rep_{1,2}.fastq.gz"

    script:
    sample_type = reads[0].baseName.split('_')[0]
    """
    for level in 10000000 1000000 100000 10000 1000; do
        for replicate in 1 2 3; do
            reformat.sh \
                in1=${reads[0]} \
                in2=${reads[1]} \
                out=${sample_type}_\${level}reads_\${replicate}rep_#.fastq.gz \
                samplereadstarget=\${level} 
        done
    done
    """
} 


