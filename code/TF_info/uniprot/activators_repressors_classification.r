rm(list = ls())
home = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/'
setwd(home)
source('code/utils.r')

data = read.delim('data/TF_info/regulation_type/uniprot/raw/uniprot_transcriptionFactor_functionCC.tab', header = T, stringsAsFactors = F)
data$Gene.names = sapply(strsplit(data$Gene.names, split = ' '), head, 1)

TFs = load_TFs_census()
data = subset(data, Gene.names %in% TFs )
data$is_Trepressor = F
data$is_Tactivator = F


# Activators
data$is_Tactivator[ grep('activates transcription', data$Function..CC., ignore.case = T) ] = T
data$is_Tactivator[ grep('activate transcription', data$Function..CC., ignore.case = T) ] = T
data$is_Tactivator[ grep('activates the transcription', data$Function..CC., ignore.case = T) ] = T
data$is_Tactivator[ grep('activate the transcription', data$Function..CC., ignore.case = T) ] = T
data$is_Tactivator[ grep('activates expression', data$Function..CC., ignore.case = T) ] = T
data$is_Tactivator[ grep('activate expression', data$Function..CC., ignore.case = T) ] = T
data$is_Tactivator[ grep('activates the expression', data$Function..CC., ignore.case = T) ] = T
data$is_Tactivator[ grep('activate the expression', data$Function..CC., ignore.case = T) ] = T
data$is_Tactivator[ grep('induces the expression', data$Function..CC., ignore.case = T) ] = T
data$is_Tactivator[ grep('induce the expression', data$Function..CC., ignore.case = T) ] = T
data$is_Tactivator[ grep('activates gene expression', data$Function..CC., ignore.case = T) ] = T
data$is_Tactivator[ grep('activate gene expression', data$Function..CC., ignore.case = T) ] = T
data$is_Tactivator[ grep('transcriptional activation', data$Function..CC., ignore.case = T) ] = T
data$is_Tactivator[ grep('transcriptional activator', data$Function..CC., ignore.case = T) ] = T
data$is_Tactivator[ grep('transcription activator', data$Function..CC., ignore.case = T) ] = T
data$is_Tactivator[ grep('ligand-activated transcription', data$Function..CC., ignore.case = T) ] = T
table(data$is_Tactivator)



# Both
data$is_Tactivator[ grep('activate or repress transcription', data$Function..CC., ignore.case = T) ] = T
data$is_Trepressor[ grep('activate or repress transcription', data$Function..CC., ignore.case = T) ] = T
data$is_Tactivator[ grep('activate or repress the transcription', data$Function..CC., ignore.case = T) ] = T
data$is_Trepressor[ grep('activate or repress the transcription', data$Function..CC., ignore.case = T) ] = T
data$is_Tactivator[ grep('activate and repress transcription', data$Function..CC., ignore.case = T) ] = T
data$is_Trepressor[ grep('activate and repress transcription', data$Function..CC., ignore.case = T) ] = T
data$is_Tactivator[ grep('activate and repress the transcription', data$Function..CC., ignore.case = T) ] = T
data$is_Trepressor[ grep('activate and repress the transcription', data$Function..CC., ignore.case = T) ] = T
data$is_Tactivator[ grep('transcriptional activator or repress', data$Function..CC., ignore.case = T) ] = T
data$is_Trepressor[ grep('transcriptional activator or repress', data$Function..CC., ignore.case = T) ] = T
data$is_Tactivator[ grep('transcriptional activator and repress', data$Function..CC., ignore.case = T) ] = T
data$is_Trepressor[ grep('transcriptional activator and repress', data$Function..CC., ignore.case = T) ] = T
data$is_Tactivator[ grep('transcriptional activation and repress', data$Function..CC., ignore.case = T) ] = T
data$is_Trepressor[ grep('transcriptional activation and repress', data$Function..CC., ignore.case = T) ] = T
data$is_Tactivator[ grep('transcriptional activation or repress', data$Function..CC., ignore.case = T) ] = T
data$is_Trepressor[ grep('transcriptional activation or repress', data$Function..CC., ignore.case = T) ] = T
data$is_Tactivator[ grep('transcription repressor and activ', data$Function..CC., ignore.case = T) ] = T
data$is_Trepressor[ grep('transcription repressor and activ', data$Function..CC., ignore.case = T) ] = T
data$is_Tactivator[ grep('stimulates or represses gene transcription', data$Function..CC., ignore.case = T) ] = T
data$is_Trepressor[ grep('stimulates or represses gene transcription', data$Function..CC., ignore.case = T) ] = T
table(data$is_Trepressor)
table(data$is_Tactivator)



# Repressor (Not co-repressors)
data$is_Trepressor[ grep('represses transcription', data$Function..CC., ignore.case = T) ] = T
data$is_Trepressor[ grep('represses its transcription', data$Function..CC., ignore.case = T) ] = T
data$is_Trepressor[ grep('repress transcription', data$Function..CC., ignore.case = T) ] = T
data$is_Trepressor[ grep('represses the transcription', data$Function..CC., ignore.case = T) ] = T
data$is_Trepressor[ grep('repress the transcription', data$Function..CC., ignore.case = T) ] = T
data$is_Trepressor[ grep('represses expression', data$Function..CC., ignore.case = T) ] = T
data$is_Trepressor[ grep('repress expression', data$Function..CC., ignore.case = T) ] = T
data$is_Trepressor[ grep('represses the expression', data$Function..CC., ignore.case = T) ] = T
data$is_Trepressor[ grep('repress the expression', data$Function..CC., ignore.case = T) ] = T
data$is_Trepressor[ grep('inhibits the expression', data$Function..CC., ignore.case = T) ] = T
data$is_Trepressor[ grep('inhibits the expression', data$Function..CC., ignore.case = T) ] = T
data$is_Trepressor[ grep('represses gene expression', data$Function..CC., ignore.case = T) ] = T
data$is_Trepressor[ grep('repress gene expression', data$Function..CC., ignore.case = T) ] = T
data$is_Trepressor[ grep('transcriptional repressor', data$Function..CC., ignore.case = T) ] = T
data$is_Trepressor[ grep('transcription repressor', data$Function..CC., ignore.case = T) ] = T
table(data$is_Trepressor)



data$Veredict_noncurated = 'none'
data$Veredict_noncurated[ data$is_Tactivator ] = 'A'
data$Veredict_noncurated[ data$is_Trepressor ] = 'R'
data$Veredict_noncurated[ data$is_Tactivator & data$is_Trepressor ] = 'A,R'
data = data[ order(data$Gene.names), ]
write.csv(data[, c('Entry', 'Gene.names', 'Veredict_noncurated', 'Function..CC.') ], file = 'data/TF_info/regulation_type/uniprot/activators_repressors_noncurated.csv', row.names = F)
