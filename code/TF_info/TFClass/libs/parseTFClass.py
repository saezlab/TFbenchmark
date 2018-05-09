#!/usr/bin/env python2
# coding: utf-8


import itertools
import re
import urllib2, urllib
from collections import defaultdict
from bs4 import BeautifulSoup




def ParseTFClassHTMLurl(outfile, URL):

	out_file = open(outfile, 'w')
	header = [ 'uniprot', 'name', 'superclass_id', 'superclass_name', 'class_id', 'class_name', 'family_id', 'family_name', 'subfamily_id', 'subfamily_name' ]
	out_file.write('\t'.join(header) + '\n')


	print 'Loading Uniprot 2 Gene Symbol annotation'
	uniprot2name = LoadUpniproXref()

	url = URL.format( URL )
	print 'Quering ' +  url	
	req = urllib2.Request(url)
	data = urllib2.urlopen(req).read()
	soup = BeautifulSoup(data, "html.parser")
	for tr in soup.find_all('tr'):
		TDs = tr.find_all('td')
		if tr.has_attr('class'):
			if 'superclass_tr' in tr['class']:
				for td in TDs:
					if td.has_attr('class'):
						if 'superclass_no_td' in td['class']:
							superclass_id = td.get_text().strip().encode('ascii', 'ignore')
						if 'superclass_descr_td' in td['class']:
							superclass_name = td.get_text().replace('Superclass: ', '').strip().replace(' ', '_').encode('ascii', 'ignore')
			if 'class_tr' in tr['class']:
				for td in TDs:
					if td.has_attr('class'):
						if 'class_no_td' in td['class']:
							class_id = td.get_text().strip().encode('ascii', 'ignore')
					if td.has_attr('colspan'):
						if 'Class' in td.get_text():
							class_name = td.get_text().replace('Class: ', '').strip().replace(' ', '_').encode('ascii', 'ignore')
			if 'family_tr' in tr['class']:
				for td in TDs:
					if td.has_attr('class'):
						if 'family_no_td' in td['class']:
							family_id = td.get_text().strip().encode('ascii', 'ignore')
					if td.has_attr('colspan'):
						if 'Family' in td.get_text():
							family_name = td.get_text().replace('Family: ', '').strip().replace(' ', '_').encode('ascii', 'ignore')
			if 'subfamily_tr' in tr['class']:
				for td in TDs:
					if td.has_attr('class'):
						if 'subfamily_no_td' in td['class']:
							subfamily_id = td.get_text().strip().encode('ascii', 'ignore')
					if td.has_attr('colspan'):
						if 'Subfamily' in td.get_text():
							subfamily_name = td.get_text().replace('Subfamily: ', '').strip().replace(' ', '_').encode('ascii', 'ignore')
			
			if 'genus_tr' in tr['class']:
				for td in TDs:
					if td.has_attr('class'):
						if 'uniprot' in td['class']:
							x = td.find_all('a')[0]['href'].replace('http://www.uniprot.org/uniprot/', '')
							uniprot = x.encode('utf-8').strip()
							name = 'unknown'
							if uniprot in uniprot2name:
								name = uniprot2name[uniprot]

							line = [ uniprot, name, superclass_id, superclass_name, class_id, class_name, family_id, family_name, subfamily_id, subfamily_name ]
							out_file.write('\t'.join(line) + '\n')
	out_file.close()





def LoadUpniproXref():
	uniprot2name = {}
	file = '/Volumes/GoogleDrive/My Drive/databases/uniprot/hsa_allids_reviewed2'
	with open(file, 'r') as f:
		for line in f:
			s = line.strip().split('\t')
			uniprot2name[ s[0] ] = s[4]
	f.close()
	return uniprot2name