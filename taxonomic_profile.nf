#!/usr/bin/env nextflow
// vim: syntax=groovy expandtab
/****************************************
 * CTMR taxonomic composition estimation
 * using multiple methods 
 * Copyright (c) Authors 2017 
 * Authors:
 *  Fredrik Boulund <fredrik.boulund@ki.se>
 *  Luisa Hugerth <luisa.warchavchik.hugerth@ki.se>
 ****************************************/

// Create file objects
kaiju_db = file(params.kaiju_db)
kaiju_nodes = file(params.kaiju_nodes)
kaiju_names = file(params.kaiju_names)
metaphlan_pickle = file(params.metaphlan_pickle)

// Channels with paired input reads
Channel
    .fromFilePairs(params.input_reads)
    .ifEmpty{ exit 1, "Found no input reads, did you specify --input_reads? I got: '${params.input_reads}'"}
    .into {input_reads_kaiju;
           input_reads_metaphlan2;
           input_centrifuge}


/****************************************
 *    TAXONOMIC COMPOSITION ESTIMATION
 ****************************************/
process kaiju {
    tag {pair_id}
    publishDir "${params.outdir}/kaiju", mode: 'move'

    input:
    set pair_id, file(reads) from input_reads_kaiju
    file database from kaiju_db
    file nodes from kaiju_nodes
    file names from kaiju_names

    output:
    file "${pair_id}.kaiju"
    file "${pair_id}.krona"
    file "${pair_id}.krona.html"
    file "${pair_id}.summary.{family,genus,species}"

    """
    rename.sh \
        in1=${reads[0]} \
        in2=${reads[1]} \
        out=${pair_id}_#.fq \
        prefix="" 

    kaiju \
        -z ${task.cpus} \
        -t $nodes \
        -f $database \
        -i ${pair_id}_1.fq \
        -j ${pair_id}_2.fq \
        -o ${pair_id}.kaiju

    kaiju2krona \
        -t $nodes \
        -n $names \
        -i ${pair_id}.kaiju \
        -o ${pair_id}.krona

    ktImportText \
        -o ${pair_id}.krona.html \
        ${pair_id}.krona

    for rank in family genus species; do
        kaijuReport \
            -t $nodes \
            -n $names \
            -i ${pair_id}.kaiju \
            -r \$rank \
            -o ${pair_id}.summary.\$rank
    done
    """
} 


process metaphlan2 {
    tag {pair_id}
    publishDir "${params.outdir}/metaphlan2", mode: 'copy'

    input:
    set pair_id, file(reads) from input_reads_metaphlan2
    file mpa_pickle from metaphlan_pickle

    output:
    file "${pair_id}.bowtie2.bz2"
    file "${pair_id}.metaphlan2.txt" 
    set pair_id, file("${pair_id}.krona") into input_metaphlan_krona

    """
    source activate metaphlan2
    metaphlan2.py \
        --nproc ${task.cpus} \
        --bowtie2out ${pair_id}.bowtie2.bz2 \
        --input_type fastq \
        --mpa_pkl $mpa_pickle \
        --bowtie2db ${params.metaphlan_bowtie2_db} \
        --sample_id ${pair_id} \
        ${reads[0]},${reads[1]} \
        ${pair_id}.metaphlan2.txt 
    metaphlan2krona.py \
        -p ${pair_id}.metaphlan2.txt \
        -k ${pair_id}.krona
    """
}


process make_metaphlan_krona {
    tag {pair_id}
    publishDir "${params.outdir}/metaphlan2", mode: 'move'

    input:
    set pair_id, file(metaphlan2_krona) from input_metaphlan_krona

    output:
    file "${pair_id}.krona.html" 

    """
    ktImportText \
        -o ${pair_id}.krona.html \
        $metaphlan2_krona
    """
}


process centrifuge {
    tag {pair_id}
    publishDir "${params.outdir}/centrifuge", mode: 'move'

    input:
    set pair_id, file(reads) from input_centrifuge

    output:
    file "${pair_id}.centrifuge.tab"
    file "${pair_id}.centrifuge_report.tsv"

    """
    centrifuge \
        --time \
        --threads ${task.cpus} \
        -x ${params.centrifuge_index_prefix} \
        -1 ${reads[0]} \
        -2 ${reads[1]} \
        -S ${pair_id}.centrifuge.tab \
        --report-file ${pair_id}.centrifuge_report.tsv 
    """
}
