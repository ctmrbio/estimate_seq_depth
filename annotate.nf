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
    set pair_id, "${pair_id}.domtbl.txt" into input_tigrfam_annotations
    file "${pair_id}.tbl.txt"
    file "${pair_id}.hmmsearch.stdout"

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
        --domtblout ${pair_id}.domtbl.txt \
        ${params.tigrfams_lib} \
        ${pair_id}_6translated.fa \
        > ${pair_id}.hmmsearch.stdout
    """
} 

process parse_annotations {
    tag {pair_id}
    publishDir "${params.outdir}/tigrfam_annotations", mode: 'move'

    input:
    set pair_id, file(domain_table) from input_tigrfam_annotations

    output:
    file "${pair_id}.tigrfam_counts.tsv"

    script:
    """
    count_tigrfam_annotations.py \
        --tbl ${domain_table} \
        --cutoffs ${params.tigrfam_cutoffs} \
        --output ${pair_id}.tigrfam_counts.tsv
    """
}



