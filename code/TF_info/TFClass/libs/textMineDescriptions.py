#!/usr/bin/env python2
# coding: utf-8


import re




def ParseActiveSitesDescriptions(accession2name,accession2descriptions, outfileP):
	outf = open(outfileP, 'w')
	for prot in accession2descriptions:
		print ' - ' + prot
		for activeSite in accession2descriptions[ prot ]:
			newlines = [prot] +  [accession2name[prot]] + [activeSite]
			outf.write('\t'.join(newlines) + '\n')
	outf.closed

