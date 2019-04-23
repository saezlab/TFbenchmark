
library(GenomicFeatures)
library(org.Hs.eg.db)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(Biostrings)
library(biomaRt)
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
annot = getBM(attributes=c("hgnc_symbol", "entrezgene", 'ensembl_gene_id', 'gene_biotype', 'ensembl_transcript_id', 'uniprotswissprot'), mart = ensembl)
home = '~/Google Drive/projects/TFbenchmark/'
setwd(home)


download_human_TSS = function(txdb, genome, upstream=1000, downstream=200, out.file){
  # Get promoter sequences
  prom = suppressWarnings(promoters(txdb, upstream=upstream, downstream=downstream, columns=c("tx_name", "TXCHROM", "TXSTART", "TXEND", "TXSTRAND", "TXTYPE", "gene_id")))
  prom = trim(prom)
  genome = getBSgenome(genome)
  prom_seqs = getSeq(genome, prom)
  idx = which( ! is.na(annot$hgnc_symbol[ match(as.character(mcols(prom)$gene_id), annot$entrezgene) ]) & ! is.na(as.character(mcols(prom)$gene_id))   )
  write.table(paste('>', as.character(mcols(prom)$tx_name)[idx], ';', 
                    as.character(mcols(prom)$gene_id)[idx],  ';', 
                    as.character(mcols(prom)$TXCHROM)[idx],  ':', 
                    as.character(mcols(prom)$TXSTART)[idx],  '-', 
                    as.character(mcols(prom)$TXEND)[idx],  ';', 
                    as.character(mcols(prom)$TXSTRAND)[idx],  ';', 
                    annot$hgnc_symbol[ match(as.character(mcols(prom)$gene_id), annot$entrezgene) ][idx], '\n',
                    unlist(sapply(prom_seqs[idx], as.character)), sep = ""), 
              file = out.file, quote = F, col.names = F, row.names = F)
}




download_human_TSS(txdb = TxDb.Hsapiens.UCSC.hg38.knownGene, genome = "hg38", upstream=1000, downstream=250, out.file = 'data/TF_target_sources/TFBS_scanning/promoter_sequences/Hsapiens.UCSC.hg38.knownGene.1000up.250down.hgncsymbol.fasta')

