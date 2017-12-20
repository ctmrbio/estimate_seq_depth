#!/usr/bin/bash 
# Compute N50 length of assembly
# Author: Luisa Hugerth, Fredrik Boulund
# Date: 2017-12-20

set -eou pipefail

if [ $# -eq 0 ]
then
	echo "Compute N50 length of assembly produced by MegaHIT"
	echo "usage: n50.sh contigs.fasta"
	exit 1
fi

echo -n "$1	"
lim=$(grep ">" $1 | cut -f4 -d= | sort -nr | datamash sum 1)
grep ">" $1  | cut -f4 -d= | sort -n | awk -v mean=$lim '{SUM+=$0; if(SUM>=mean){print $0; exit}}'
