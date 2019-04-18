#!/usr/bin/env python2
# coding: utf-8


import itertools
import re
import urllib2, urllib
from collections import defaultdict
from bs4 import BeautifulSoup




def GetActiveSitesDescriptions(uniprots=None, infile=None):
	
	if uniprots:
		uniprots = uniprots.split(',')
		print uniprots
	elif infile:
		uniprots = loadFileasList(infile)
		print ' - Loading proteins from:' + infile


	accession2name = {}
	accession2descriptions = defaultdict(list)

	URL = "http://www.uniprot.org/uniprot/{}.xml"
	for uniprot in uniprots:	
		uniprot = uniprot.upper()
		url = URL.format( uniprot )
		print ' - Reading ' +  url	
		req = urllib2.Request(url)
		data = urllib2.urlopen(req).read()
		soup = BeautifulSoup(data)
		for entry in soup.find_all('entry'):
			if entry.gene is not None:
				# Extract descriptions incuding PTMs
				features = [f for f in entry.find_all('feature') if f['type'] in ['active site'] ]
				if len(features) > 0:
					print '   Number of features: ' + str(len(features))
					for ft in features:
						activeSite = parseComment(ft)
						# Save accession
						accession =  str(entry.accession.get_text())
						# Save activeSite
						accession2descriptions[ accession ].append(activeSite)
						# Save gene names
						names = [ str(g.get_text()) for g in entry.gene.find_all("name")]
						accession2name[ accession ] = names[0]
	return (accession2name,accession2descriptions)



def parseComment(feature):
	activeSite = []
	if feature.find_all('location') is not None:
		if feature.find('location').find('position') is not None:
			activeSite.append(feature.find('location').find('position')['position'])
	else:
		activeSite.append(0)
	if feature.has_attr('description'):
		activeSite.append(feature['description'])
	else:
		activeSite.append('Unknown')
	return '\t'.join(activeSite)


def loadFileasList(file):
	l = [line.strip() for line in open(file)]
	return l