# estimate_seq_depth
A set of small workflows that:

1. simulates metagenomes for different sample types
2. runs taxonomic profiling 
3. runs basic functional profiling

It's intended use is to help estimate the required sequencing depth for
different type of metagenomic microbiome samples.

# Run
Run workflows using Nextflow, e.g.:

    nextflow run path/to/estimate_seq_depth/taxonomic_profile.nf -profile ctmr_nas --input_reads 'path/to/reads/*_{1,2}.fastq.gz'

This will start the taxonomic profiling workflow using the `ctmr_nas` profile, as defined in `conf/ctmr_nas.config`. 

## Dependencies
The workflows uses a lot of external scripts:

* BBMap (`bbmap.sh`, `randomreads.sh`, `reformat.sh`)
* Python packages, e.g. `BioPython`, `NumPy`, ...
* FastQC
* Kaiju
* MetaPhlAn2 
* Centrifuge

