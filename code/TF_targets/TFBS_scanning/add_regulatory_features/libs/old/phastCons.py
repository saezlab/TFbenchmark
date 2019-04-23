#!/usr/bin/env python2
# coding: utf-8


import re
from os import listdir
from os.path import isfile, join
from pycellbase.cbclient import CellBaseClient
cbc = CellBaseClient()
cbcregion = cbc.get_genomic_region_client()


def GetConservation (region):
	query = str(region[0]) + ":" + str(region[1]) + "-" + str(region[2])
	conservation = cbcregion.get_conservation(query, assembly="GRCh38")
	for result in conservation[0]['result']:
		if 'phastCons' in result['source']:
			phastCons = result['values']
		if 'phylop' in result['source']:
			phylop = result['values']
	positions = [ str(i) for i in range(region[1], region[2], 1) ]
	phastCons_coverage = float(len(phastCons))/float(len(positions))
	phylop_coverage = float(len(phylop))/float(len(positions))
	return phastCons,phylop,phastCons_coverage,phylop_coverage



def GetPhastCons_old (region, PhastConscores):
	positions = [str(region[0]) + '_' + str(i) for i in range(region[1], region[2], 1) ]
	scores = [float(PhastConscores[position]) for position in positions if position in PhastConscores ]
	coverage = float(len(scores))/float(len(positions))
	return scores,coverage


def loadPhastConscores_new2():
	phastCons_file = '/Users/luzgaral/data/phastCons/hg38.phastCons20way.wigFix'
	PhastConscores = {}
	with open(phastCons_file) as f:
		for line in f:
			if 'fixedStep' in line:
				position = int(line.strip().split()[2].replace('start=', ''))
				ch = line.strip().split()[1].replace('chrom=chr', '')
			else:
				position += 1
				key = ch + '_' + str(position)
				val = line.strip()
				PhastConscores[key] = val
	return PhastConscores



def loadPhastConscores_new():
	phastCons_file = '/Users/luzgaral/data/phastCons/hg38.phastCons20way.wigFix'
	PhastConscores = {}
	with open(phastCons_file) as f:
		for line in f:
			if 'fixedStep' in line:
				position = int(line.strip().split()[2].replace('start=', ''))
				ch = line.strip().split()[1].replace('chrom=chr', '')
			else:
				position += 1
				key = ch + '_' + str(position)
				val = line.strip()
				PhastConscores[key] = val
	return PhastConscores


def loadPhastConscores():
	print 'Annotating promoters only'
	promoter_positions = set()
	promoters_file = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/data/TF_target_sources/TFBS_scanning/promoter_sequences/Hsapiens.UCSC.hg38.knownGene.1000up.250down.hgncsymbol.fasta'
	with open(promoters_file) as f:
		for line in f:
			if '>' in line:
				s = line.strip().split(';')[2]
				ch = s.split(':')[0].replace('chr', '')
				start = int(s.split(':')[1].split('-')[0])
				end = int(s.split(':')[1].split('-')[1])
				for i in range(start, end, 1):
					region = ch + '_' + str(i)
					promoter_positions.add(region)
	phastCons_path = '/Users/luzgaral/data/phastCons/'
	phastCons_files = [f for f in listdir(phastCons_path) if isfile(join(phastCons_path, f))]
	PhastConscores = collections.defaultdict(float)
	for fi in phastCons_files:
		if 'chr' in fi:
			print '\t-' + fi
			with open(phastCons_path+fi) as f:
				for line in f:
					if 'fixedStep' in line:
						position = int(line.strip().split()[2].replace('start=', ''))
						ch = line.strip().split()[1].replace('chrom=chr', '')
					else:
						position += 1
						key = ch + '_' + str(position)
						if key in promoter_positions:
							val = line.strip()
							PhastConscores[key] = val
	return PhastConscores
