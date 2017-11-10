#!/bin/bash 
# Script to generate simulated metagenomes for different sample types
# Fredrik Boulund and Luisa Hugerth 2017

# "Unofficial Bash strict mode"
set -euo pipefail
IFS=$'\n\t'

# Define input arguments
if [ $# -eq 0 ]; then
	echo "usage: script.sh FASTA N"
	echo ""
	echo "  FASTA   path to a single FASTA file with reference"
	echo "          genome sequences. The file basename is split on"
	echo '          the underscore character (_) to produce the basename'
	echo "          for the output filenames, e.g. 'biopsies_m1_f500.fasta'"
	echo "          produces output files 'biopsies_50M_{1,2}.fastq.gz'."
	echo '  N       a unique BBMap database build number (int) for'
	echo "          the BBMap index of genomes.fasta."
	echo ""
	echo "The script is designed to be easy to run with GNU parallel, e.g.:"
	echo "parallel './run_randomreads.sh {} {#}' ::: ../type1/type1_m1_f500.fasta ../type2/type2_m1_f500.fasta"
	exit 1
fi
path_to_genomes_fasta=$1
db_build_no=$2

# Set mutation rates
snprate=0.03
insrate=0.01
delrate=0.01
subrate=0.01
length=125
reads=50000000

# Figure out a good output name 
sample_type_basename=$(basename $path_to_genomes_fasta)
clean_type_name=${sample_type_basename%%_*}

# Generate reads
randomreads.sh \
	build=$db_build_no \
	ref=$path_to_genomes_fasta \
	out=${clean_type_name}_50M_#.fastq.gz \
	snprate=$snprate \
	insrate=$insrate \
	delrate=$delrate \
	subrate=$subrate \
	length=$length \
	paired=true \
	reads=$reads 

fastqc ${clean_type_name}_50M_1.fastq.gz 
fastqc ${clean_type_name}_50M_2.fastq.gz 
