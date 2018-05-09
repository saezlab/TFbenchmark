#!/usr/bin/env python2
# coding: utf-8



def MatchENSR (region, RegulatoryFeatures_coord,RegulatoryFeatures_name):
	positions = [str(region[0]) + '_' + str(i) for i in range(region[1], region[2], 1) ]
	features = [RegulatoryFeatures_coord[position] for position in positions if position in RegulatoryFeatures_coord ]
	u_features =  list(set(features))
	u_featuretype = [ RegulatoryFeatures_name[fe] for fe in u_features if fe in RegulatoryFeatures_name  ]
	coverage = float(len(features))/float(len(positions))
	return u_features,u_featuretype,coverage




def loadRegulatoryFeatures():
	RegulatoryFeatures_coord = {}
	RegulatoryFeatures_name = {}
	regulation_file = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/data/genomic_regulatory_regions/ensembl/human_regulatory_features_GRCh38.p10.txt'
	with open(regulation_file) as f:
		next(f)
		for line in f:
			(feature_id, ch, start, end, feature_type) = line.strip().split('\t')
			RegulatoryFeatures_name[feature_id] = feature_type
			for i in range(int(start), int(end), 1):
				RegulatoryFeatures_coord[str(ch)+'_'+str(i)] = feature_id
	return RegulatoryFeatures_coord,RegulatoryFeatures_name


