#!/usr/bin/env python2
# coding: utf-8



import sys
import numpy
from libs.ensembl_regulation import MatchENSR,loadRegulatoryFeatures
from libs.conservation import GetConservation
from libs.regulation import GetRegulatory
from collections import defaultdict
from promoters_xref import loadPromotersXRef


def AnnotateTFBS (infile, outfile):
	print 'Loading ENSMBL regulatory regions'
	(RegulatoryFeatures_coord,RegulatoryFeatures_name) = loadRegulatoryFeatures()
	PromotersXRef = loadPromotersXRef()
	conservation_bucket = defaultdict(list)
	regulatory_bucket = defaultdict(list)
	TFs_bucket = defaultdict(list)
	print 'Annotating TFBS from fimo'
	outf = open(outfile, 'w')
	with open(infile) as f:
		header = f.readline().strip()
		outf.write(header + '\tgene\tchr\tgenome_start\tgenome_end\tPhylop_scores\tPhylop_coverage\tPhastCons_scores\tPhastCons_coverage\tENSR\tENSR_type\tENSR_coverage\tccat_histone\tEnsembl_Regulatory_Build\tSWEmbl_R0025\tSWEmbl_R015\tSWEmbl_R0005_IDR\n')
		for line in f:
			if '_alt' in line:
				continue
			content = line.strip()
			# annotate only top 5000 hits
			TF = content.split('\t')[0]
			if TF in TFs_bucket:
				TFs_bucket[TF] += 1
			else:
				TFs_bucket[TF] = 1
			if TFs_bucket[TF] > 5000:
				continue
			print 'REGION: ' + content
			# obtain informatifon from the target region
			transcript = content.split('\t')[2].split(';')[0]
			gene = PromotersXRef[transcript].split(';')[4]
			genomic_coordinates = getGenomicCoordinates_fimo_v4_12(PromotersXRef[transcript], content)
			ID = '_'.join(str(x) for x in genomic_coordinates)
			# retrieve regulatory and conservation information
			if ID not in conservation_bucket:
				#print '- Retrieving data from Cellbase ...'
				#print '\tConservation'
				conservation_bucket[ID] = GetConservation(genomic_coordinates)
				#print '\tRegulation'
				regulatory_bucket[ID] = GetRegulatory(genomic_coordinates)

			#print '... done!\n'
			(phastCons,phylop,phastCons_coverage,phylop_coverage) = conservation_bucket[ID]
			(features,featuretypes,REGcoverage) = MatchENSR (genomic_coordinates, RegulatoryFeatures_coord,RegulatoryFeatures_name)
			regulatory = regulatory_bucket[ID]
			# print in outfile
			out = content.split('\t')
			out.append(gene)
			out.append(str(genomic_coordinates[0]))
			out.append(str(genomic_coordinates[1]))
			out.append(str(genomic_coordinates[2]))
			if len(phylop) == 0:
				phylop = numpy.NaN
			out.append(str(numpy.nanpercentile(phylop, 75)))
			out.append("{:10.2f}".format(phylop_coverage))
			if len(phastCons) == 0:
				phastCons = numpy.NaN
			out.append(str(numpy.nanpercentile(phastCons, 75)))
			out.append("{:10.2f}".format(phastCons_coverage))
			out.append(','.join(features))
			out.append(','.join(featuretypes))
			out.append("{:10.2f}".format(REGcoverage))
			outf.write('\t'.join(out))
			outf.write('\t' + '\t'.join(regulatory) + '\n')
	outf.closed




def getGenomicCoordinates_fimo_v4_11_4(promoter_info):
	genomic_coords_promoter = promoter_info.split(';')[2].split(':')
	ch = genomic_coords_promoter[0].replace('chr', '')
	start = int(genomic_coords_promoter[1].split('-')[0])
	rel_start = int(s[3]) - 1
	rel_end = int(s[4]) - 1
	return ch, start + rel_start, start + rel_end





def getGenomicCoordinates_fimo_v4_12(promoter_info, content):
	genomic_coords_promoter = promoter_info.split(';')[2].split(':')
	ch = genomic_coords_promoter[0].replace('chr', '')
	start = int(genomic_coords_promoter[1].split('-')[0])
	end = int(genomic_coords_promoter[1].split('-')[1])	
	start_match = int(content.split('\t')[3])
	end_match = int(content.split('\t')[4])
	if( start_match < start or end_match > end ):
		print 'Warning! Coordinates error at line: ' + content + '\n'
	return ch, start_match, end_match


