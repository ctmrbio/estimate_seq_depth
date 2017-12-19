#!/usr/bin/env nextflow
// vim: syntax=groovy expandtab
/****************************************
 * Annotate reference contigs using TIGRFAM
 * Copyright (c) Authors 2017 
 * Authors:
 *  Fredrik Boulund <fredrik.boulund@ki.se>
 *  Luisa Hugerth <luisa.warchavchik.hugerth@ki.se>
 ****************************************/

Channel
    .fromPath(params.input_contigs)
    .ifEmpty{ exit 1, "Found no input contigs, did you specify --input_contigs? I got: '${params.input_contigs}'"}
    .set {input_contigs_mgm}


/****************************************
 *                TIGRFAMs
 ****************************************/
process orf_prediction {
    tag {file_id}
    publishDir "${params.outdir}/metagenemark", mode: 'copy'

    input:
    file(contigs) from input_contigs_mgm

    output:
    set file_id, file("${file_id}.predicted_proteins.fasta") into input_proteins_tigrfam
    file "${file_id}.predicted_nucleotides.fasta"
    file "${file_id}.gmhmmp.lst"

    script:
    file_id = contigs.baseName
    """
    gmhmmp \
        -r \
        -m ${params.mgm_model} \
        -o ${file_id}.gmhmmp.lst \
        -A ${file_id}.predicted_proteins.fasta \
        -D ${file_id}.predicted_nucleotides.fasta \
        ${contigs} 
    """
}


process hmmsearch_tigrfam {
    tag {file_id}
    publishDir "${params.outdir}/hmmsearch_tigrfam", mode: 'copy'

    input:
    set file_id, file(proteins) from input_proteins_tigrfam

    output:
    file "${file_id}.tbl.txt"
    file "${file_id}.hmmsearch.stdout"
    file "${file_id}.tigrfam_counts.tsv"

    script:
    """
    hmmsearch \
        --cpu ${task.cpus} \
        --seed 1337 \
        --tblout ${file_id}.tbl.txt \
        ${params.tigrfams_lib} \
        ${proteins} \
        > ${file_id}.hmmsearch.stdout \
    && \
    count_tigrfam_annotations.py \
        --tbl ${file_id}.tbl.txt \
        --cutoffs ${params.tigrfam_cutoffs} \
        --output ${file_id}.tigrfam_counts.tsv
    """
} 
