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
    .fromFilePairs(params.input_reads)
    .ifEmpty{ exit 1, "Found no input read pairs, did you specify --input_reads? I got: '${params.input_reads}'"}
    .into {input_reads_megahit;
           input_reads_map}


process assemble {
    tag {pair_id}
    publishDir "${params.outdir}/assembled_contigs", mode: 'copy'

    input:
    set pair_id, file(reads) from input_reads_megahit

    output:
    set pair_id, file("${pair_id}.contigs.fasta") into input_contigs_mgm

    script:
    """
    megahit \
        -1 ${reads[0]} \
        -2 ${reads[1]} \
        -o megahit_out
    mv megahit_out/final.contigs.fa ${pair_id}.contigs.fasta
    """
}


process orf_prediction {
    tag {file_id}
    publishDir "${params.outdir}/metagenemark_contigs", mode: 'copy'

    input:
    set file_id, file(contigs) from input_contigs_mgm

    output:
    set file_id, file("${file_id}.predicted_proteins.fasta") into input_proteins_tigrfam
    set file_id, file("${file_id}.predicted_nucleotides.fasta") into orf_nucleotides
    file "${file_id}.gmhmmp.lst"

    script:
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


// Match up predicted orfs with the correct set of reads
// by grouping the two channels on their pair_id (position 0 in the tuple)
orf_nucleotides
    .combine(input_reads_map, by: 0) 
    .set {input_orfs_reads}


process map_reads_to_orfs {
    tag {pair_id}
    publishDir "${params.outdir}/metagenemark_contigs", mode: 'copy'

    input:
    set pair_id, file(orfs), file(reads) from input_orfs_reads

    output:
    file "${pair_id}.mapped_reads.sam.gz" 
    file "${pair_id}.covstats.txt" 
    file "${pair_id}.rpkm.txt" 

    script:
    """
    bbmap.sh \
        in1=${reads[0]} \
        in2=${reads[1]} \
        ref=${orfs} \
        nodisk=t \
        out=${pair_id}.mapped_reads.sam.gz \
        covstats=${pair_id}.covstats.txt \
        rpkm=${pair_id}.rpkm.txt \
    """
}


process hmmsearch_tigrfam {
    tag {file_id}
    publishDir "${params.outdir}/hmmsearch_tigrfam_contigs", mode: 'copy'

    input:
    set file_id, file(proteins) from input_proteins_tigrfam

    output:
    file "${file_id}.tbl.txt"
    file "${file_id}.hmmsearch.stdout"

    script:
    """
    hmmsearch \
        --cpu ${task.cpus} \
        --seed 1337 \
        --tblout ${file_id}.tbl.txt \
        ${params.tigrfams_lib} \
        ${proteins} \
        > ${file_id}.hmmsearch.stdout \
    """
} 


