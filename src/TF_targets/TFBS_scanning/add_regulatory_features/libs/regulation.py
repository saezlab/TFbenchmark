#!/usr/bin/env python2
# coding: utf-8


from pycellbase.cbclient import CellBaseClient
cbc = CellBaseClient()
cbcregion = cbc.get_genomic_region_client()



def GetRegulatory(region):
	query = str(region[0]) + ":" + str(region[1]) + "-" + str(region[2])
	regulatory = cbcregion.get_regulatory(query, assembly="GRCh38")
	ccat_histone = list()
	Regulatory_Build = list()
	SWEmbl_R0025 = list()
	SWEmbl_R015 = list()
	SWEmbl_R0005_IDR = list()
	for result in regulatory[0]['result']:
		overlap = get_overlap(region, result)
		if float(overlap) > 0.25:
			if 'ccat_histone' in result['source']:
				ccat_histone.append(result['featureType'])
			if 'Regulatory_Build' in result['source']:
				Regulatory_Build.append(result['featureType'])
			if 'SWEmbl_R0025' in result['source']:
				SWEmbl_R0025.append(result['featureType'])
			if 'SWEmbl_R015' in result['source']:
				SWEmbl_R015.append(result['featureType'])
			if 'SWEmbl_R0005_IDR' in result['source']:
				SWEmbl_R0005_IDR.append(str(result['featureType']))

	out = [ ','.join(list(set(ccat_histone))) , 
	','.join(list(set(Regulatory_Build))),
	','.join(list(set(SWEmbl_R0025))), 
	','.join(list(set(SWEmbl_R015))),  
	','.join(list(set(SWEmbl_R0005_IDR))) ] 

	return out


def get_overlap(region, cbc_entry):
	total = region[2] - region[1] + 1
	if cbc_entry['start'] > region[2]:
		overlap = 0
	elif cbc_entry['end'] < region[1]:
		overlap = 0
	elif cbc_entry['start'] <= region[1] and cbc_entry['end'] >= region[2]:
		overlap = total
	elif cbc_entry['start'] > region[1] and cbc_entry['end'] >= region[2]:
		overlap = region[2] - cbc_entry['start'] + 1
	elif cbc_entry['start'] <= region[1] and cbc_entry['end'] < region[2]:
		overlap = cbc_entry['end'] - region[1] + 1
	elif cbc_entry['start'] > region[1] and cbc_entry['end'] < region[2]:
		overlap = cbc_entry['end'] - cbc_entry['start'] + 1
	else:
		print '-----\nproblem found\n'
		print cbc_entry
		print region
	return str("{:10.2f}".format(overlap/total)).strip()
