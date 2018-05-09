


def loadPromotersXRef():
	print 'Annotating promoters only'
	promoter_xref = dict()
	promoters_file = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/data/TF_target_sources/TFBS_scanning/promoter_sequences/Hsapiens.UCSC.hg38.knownGene.1000up.250down.hgncsymbol.fasta'
	with open(promoters_file) as f:
		for line in f:
			if '>' in line:
				full_name = line.strip().split('>')[1]
				transcript_id = full_name.split(';')[0]
				promoter_xref[transcript_id] = full_name
	return promoter_xref

