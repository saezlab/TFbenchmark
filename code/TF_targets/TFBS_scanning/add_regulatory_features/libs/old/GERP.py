#!/usr/bin/env python2
# coding: utf-8


import re



def GetGERPs (region, GERPscores):
	positions = [str(region[0]) + '_' + str(i) for i in range(region[1], region[2], 1) ]
	scores = [GERPscores[position] for position in positions if position in GERPscores ]
	coverage = float(len(scores))/float(len(positions))
	return scores,coverage



def loadGERPscores():
	GERPscores = {}
	GERPfile = '/Users/luzgaral/software/annovar/humandb/hg38_gerp.txt'
	with open(GERPfile) as f:
		next(f)
		for line in f:
			(key, val) = line.strip().split()
			GERPscores[key] = val
	return GERPscores




