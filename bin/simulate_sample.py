#!/usr/bin/env python

"""simulate-sample.py: extracts genomes from a database to simulate a natural sample.
"""

from __future__ import division
from collections import defaultdict
import argparse
import csv
import math
from Bio import SeqIO
from random import choices
from numpy import random

__author__ = "Luisa W Hugerth"
__email__ = "luisa.hugerth@scilifelab.se"


		
def parse_taxonomy(composition, minimum, factor):
	taxcounts = dict()
	with open (composition) as csvfile:
		reader = csv.reader(csvfile, delimiter="\t")
		for row in reader:
			tax=row[0]
			raw = int(row[1])
			count = int(raw/factor)
			if count < minimum:
				count = minimum
			taxcounts[tax] = count
			#print(tax, count)
	return taxcounts


def translate_taxonomy(taxcounts, taxonomy):
	genera = defaultdict(int)
	used = set()
	with open(taxonomy) as csvfile:
		reader = csv.reader(csvfile, delimiter="\t")
		for row in reader:
			fulltaxon = row[1]
			eachtaxa = fulltaxon.split(";")
			if(len(eachtaxa)>=6):
				genus = eachtaxa[5]
				if(genus != ''):
					if genus.split("_")[0] == "Candidatus":
						genus = " ".join(genus.split("_")[0:1])
					else:
						genus = genus.split("_")[0]
					for i in range(2, 7):
						name = ";".join(eachtaxa[0:i])
						name = name + ";"
						if name in taxcounts and name not in used:
							used.add(name)
							genera[genus] += taxcounts[name]
							#print(name, str(taxcounts[name]), genus, str(genera[genus]))
	return genera
				
def print_genomes(genera, genomes, noplasmid):
	ids = defaultdict(set)
	tax = dict()
	seqs = dict()
	with open(genomes) as handle:
		for record in SeqIO.parse(handle, "fasta"):
			desc = record.description.split(" ")
			genus = desc[1]
			species = " ".join(desc[1:3])
			seqid = record.id
			seq = record.seq
			seqs[seqid] = seq
			if genus in genera:
				if ('plasmid' not in desc) or (not noplasmid):
					ids[genus].add(seqid)
					tax[seqid] = species
	for genus, counts in genera.items():
		if (ids[genus]):
			#keys = choices(list(list(ids[genus]) for _ in list(range(counts))))
			#flat_keys = [item for sublist in keys for item in sublist]
			keys = random.choice(list(ids[genus]), counts)
			#print (counts, keys)
			for key in keys:
				species = tax[key]
				print("".join([">", " ".join([key, species])]))
				print(seqs[key])


def main(composition, taxonomy, genomes, minimum, factor, noplasmid):
	if (not factor):
		factor = 1
	if (not minimum):
		minimum = 0
	taxcounts = parse_taxonomy(composition, minimum, factor)	#Ok! Table makes sense
	genera = translate_taxonomy(taxcounts, taxonomy)		#Ok! Table makes sense
	print_genomes(genera, genomes, noplasmid)


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Requires a composition table, a taxonomy DB and a genome DB')
	parser.add_argument('-c', '--composition', help='TSV table with taxon\tcounts')
	parser.add_argument('-tax', '--taxonomy', help='Taxonomy database used in the generation of the composition file')
	parser.add_argument('-g', '--genomes', help="Fasta DB to extract the community from")
	parser.add_argument('-m', '--minimum', type=int, nargs='?', const=0, help="Minimal number of genomes per taxon, default 0")
	parser.add_argument('-f', '--factor', type=int, nargs='?', const=1, help="Integer to divide taxon counts by, limiting the number of genomes, default 1")
	parser.add_argument('--noplasmid', type=bool, nargs='?', const=False, help="Whether plasmids in the genome DB may be included as genomes, default False")
	args = parser.parse_args()
	
	main(args.composition, args.taxonomy, args.genomes, args.minimum, args.factor, args.noplasmid)

