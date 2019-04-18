#!/usr/bin/env python2
# coding: utf-8


from pycellbase.cbclient import CellBaseClient
cbc = CellBaseClient()
cbcregion = cbc.get_genomic_region_client()


def GetConservation (region):
	query = str(region[0]) + ":" + str(region[1]) + "-" + str(region[2])
	conservation = cbcregion.get_conservation(query, assembly="GRCh38")
	phastCons = ()
	phylop = ()
	for result in conservation[0]['result']:
		if 'phastCons' in result['source']:
			phastCons = result['values']
		if 'phylop' in result['source']:
			phylop = result['values']
	phastCons = [x for x in phastCons if x is not None]
	phylop = [x for x in phylop if x is not None]

	positions = region[2] - region[1] + 1
	phastCons_coverage = float(len(phastCons))/positions
	phylop_coverage = float(len(phylop))/positions
	return phastCons,phylop,phastCons_coverage,phylop_coverage
