#!/usr/bin/env nextflow
// vim: syntax=groovy expandtab
/****************************************
 * Annotate reference contigs using TIGRFAM
 * and Pfam models.
 * Copyright (c) Authors 2017 
 * Authors:
 *  Fredrik Boulund <fredrik.boulund@ki.se>
 *  Luisa Hugerth <luisa.warchavchik.hugerth@ki.se>
 ****************************************/

Channel
    .fromFilePairs(params.input_reads)
    .ifEmpty{ exit 1, "Found no input contigs, did you specify --input_reads? I got: '${params.input_reads}'"}
    .into {input_reads_tigrfam;
           input_reads_pfam}

params.tigrfams_lib = '/home/ctmr/db/TIGRFAMs/latest/LIB/TIGRFAMs_15.0_HMM.LIB'


/****************************************
 *                TIGRFAMs
 ****************************************/
process hmmsearch_tigrfam {
    tag {pair_id}
    publishDir "${params.outdir}/hmmsearch_tigrfam", mode: 'move'

    input:
    set pair_id, file(reads) from input_reads_tigrfam

    output:
    file "${pair_id}.tbl.txt"
    file "${pair_id}.hmmsearch.stdout"
    file "${pair_id}.tigrfam_counts.tsv"

    script:
    """
    translate6frames.sh \
        in=${reads[0]} \
        in2=${reads[1]} \
        out=${pair_id}_6translated.fa \
    && \
    hmmsearch \
        --cpu ${task.cpus} \
        --seed 1337 \
        --tblout ${pair_id}.tbl.txt \
        ${params.tigrfams_lib} \
        ${pair_id}_6translated.fa \
        > ${pair_id}.hmmsearch.stdout \
    && \
    count_tigrfam_annotations.py \
        --tbl ${pair_id}.tbl.txt \
        --cutoffs ${params.tigrfam_cutoffs} \
        --output ${pair_id}.tigrfam_counts.tsv
    """
} 
