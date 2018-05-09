rm(list = ls())
home = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/'
setwd(home)

source('code/TFperturbations/old/download_and_DGE_lib_experimental_designs.r')
library(pheatmap)


my_GEOid = 'GSE31534'#'GSE31912'



# Create directory and download dataset information from GEO
study_files_directory = paste('/Users/luzgaral/tmp/GEO/', my_GEOid, '/', sep = '')
dir.create(study_files_directory, showWarnings = F)
GSEDataobj = getGEO(my_GEOid, GSEMatrix = F, destdir = study_files_directory)
GSEid = Meta(GSEDataobj)$geo_accession
platform_id = Meta(GSEDataobj)$platform_id



# Format phenotipic data
sAnnot = t(sapply(Meta(GSEDataobj)$sample_id, function(x)  c(x, Meta( getGEO(x, destdir = study_files_directory))$title) ))
sAnnot = as.data.frame(sAnnot)
names(sAnnot)[1] = 'sample'
names(sAnnot)[2] = 'phenotype'
sAnnot$phenotype = format_phenoname(pheno_name = sAnnot$phenotype)
sAnnot$phenotype =  gsub('CONTORL3', 'CONTROL', gsub('A375_', '', gsub('_48H', '', gsub('UNTREATED[1-9]', 'UNTREATED', gsub('CONTROL[1-9]', 'CONTROL', gsub('A375_SIRNA_', '', sAnnot$phenotype))))))
sAnnot$phenotype[ sAnnot$phenotype == 'CDC2']  = 'CDK1'
sAnnot$phenotype[ sAnnot$phenotype == 'FOXO3A']  = 'FOXO3'
sAnnot$treatment = sAnnot$phenotype
sAnnot$treatment[ ! sAnnot$treatment %in% c('UNTREATED', 'CONTROL') ] = 'shRNA'
rownames(sAnnot) = sAnnot$sample
sAnnot = sAnnot[ sAnnot$sample != 'GSM782746', ] # problematic sample



# download cels, build eset and rma normalise
getGEOSuppFiles(GEO = GSEid, baseDir = study_files_directory)
cel_tarfile = list.files(study_files_directory, pattern = '_RAW.tar', full.names = T, ignore.case = T, recursive = T)
untar(cel_tarfile[1], exdir = study_files_directory, extras = '-vz')
cel_files = list.files(study_files_directory, pattern = 'cel', full.names = T, ignore.case = T)
sampleNames = sapply(list.files(study_files_directory, pattern = 'cel', ignore.case = T), function(x) unlist(strsplit(unlist(strsplit(x, split = '_')), split = '\\.'))[1] )
my_cels = ReadAffy(filenames = cel_files[ sampleNames %in% rownames(sAnnot) ], sampleNames = sampleNames[ sampleNames %in% rownames(sAnnot) ], verbose = F, compress = T)
eset = affy::rma(my_cels)



# identify the 3 hidden batches
clus = pheatmap(cor(exprs(eset)))
batch = cutree(hclust(dist(cor(exprs(eset)))), k = 3)
pc = prcomp( t ( exprs( eset ) ) )
plot( pc$x[ , 1:2 ], col = c('skyblue3', 'coral3', 'yellowgreen')[ match(batch, unique(batch))], pch = batch)



# correct batch effect
library(sva)
sAnnot$batch = batch
modcombat = model.matrix(~as.factor(phenotype), data=sAnnot)
combat_eset = ComBat(dat=exprs(eset), batch=batch, mod=modcombat, par.prior=T, prior.plots = T)




# obtain z-scores
Esymbol = combat_eset
gene_annot = find.geneannot(platform_id)
rownames(Esymbol) = getSYMBOL(rownames(exprs(eset)), data = gene_annot)
colnames(Esymbol) = paste(colnames(Esymbol), sAnnot[colnames(Esymbol),]$phenotype , sep = '_')
Esymbol = Esymbol[ ! is.na(rownames(Esymbol)), ]
mycor = cor(Esymbol)
pheatmap(mycor)
Esymbol_zscore = t(scale(t(Esymbol)))
save(Esymbol_zscore, file = 'data/regulons_QC/pertubations/GSE31534/expression_signatures/A375_zscores.rdata')



# check quality of shRNA assays
genes = intersect(sAnnot$phenotype, rownames(Esymbol_zscore))
ranks = sapply(genes, function(gene){
  ra = sort(Esymbol_zscore[gene,])
  mean(grep(paste(gene, '$', sep = ''), names(ra)))
})
df = melt(ranks)
df$gene =  rownames(df)
names(df)[1] = 'rank_position'
df$label = df$gene
df$label[ df$rank_position <= 5 ] = ''
df$color = 'gray'
df$color[ df$rank_position > 5 ] = 'skyblue3'
df$color[ df$rank_position > 10 ] = 'purple3'
df$color[ df$rank_position > 20 ] = 'red'
library(ggplot2)
library(ggrepel)
ggplot(df, aes(x=rank_position)) + geom_histogram() + 
  theme_light()
ggplot(df, aes(y=rank_position, x=gene, label=label, color = color)) + geom_point() + geom_text_repel() + 
  geom_hline(yintercept = 20, size=.1)+ 
  ylab('rank\'s position of knock out gene across samples') + 
  theme_light() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_color_identity()
ggsave(filename = 'data/regulons_QC/pertubations/GSE31534/expression_signatures/A375_zscores_ranks.png', width = 12, height = 5)

