#!/bin/bash

source activate py2.7
#1) Select genomes
head -n50000000 /home/ctmr/db/refseq/refseq_bacteria.fasta | grep '^>' | awk '{print $1}' > 1026_genome_headers.txt
perl /home/ctmr/src/amplicons/extract_sample.pl --sample=1026_genome_headers.txt --database=/home/ctmr/db/refseq/refseq_bacteria.fasta > full_target_genomes.fasta

#2) Split them up
python /home/ctmr/src/CONCOCT/scripts/cut_up_fasta.py -c 120 -o 0 full_target_genomes.fasta > target_genome_reads120bp.fa

#3) Remove potential human DNA.
##OBS! We want to remove everything that could *pontetially* be human, so rather have false positives than false negatives
##The standard procedure would be: bbmap.sh minid=0.95 maxindel=3 bwr=0.16 bw=12 quickmatch fast minhits=2 path=/home/ctmr/db/hg19/hg19_main_mask_ribo_animal_allplant_allfungus.fa   untrim -Xmx23g in=target_genome_reads120bp.fa outu=clean_reads.fq outm=human_reads.fq
bbmap.sh minid=0.8 maxindel=3 bwr=0.16 bw=12 quickmatch fast minhits=2 path=/home/ctmr/db/hg19/hg19.gz untrim -Xmx23g in=target_genome_reads120bp.fa outu=clean_reads.fq outm=human_reads.fq

##OBS! Only 43 short fragments (adapters?) were removed with either approach

#4) 
