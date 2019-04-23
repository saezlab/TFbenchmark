rm(list = ls())
home = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/'
setwd(home)
source('code/lib/data.r')
source('code/lib/utils.r')




save_network = function(sdata, filename){
  sdata = Filter(Negate(is.null), sdata)
  sdata = melt(sdata, id.vars = names(sdata[[1]]) )
  network = unique(sdata[, c('motif_alt_id', 'gene') ])
  names(network) = c('TF', 'target')
  network$TF = toupper(network$TF)
  network = subset(network, target != '')
  network = subset(network, TF %in% TF_census)
  # summarize_network(network)
  write.table(network, file = filename, sep = '\t', col.names = T, row.names = F, quote = F)
}



uniprot_annot = load_proteins_annotation(type = 'uniprot')
uniprot_annot = subset(uniprot_annot, Gene.name.primary != '')

TF_census = load_TFs_census()


for (dataset in c('jaspar_2018', 'hocomoco_v11')){
  message(rep('#', 10), ' ',  dataset, ' ', rep('#', 10))
  path = paste('data/TF_target_sources/TFBS_scanning/', dataset, '/scanning_output/', sep = '')
  data = read.delim(paste(path, 'fimo/fimo_annotated.txt', sep = ''), header = T, stringsAsFactors = F)
  if( dataset == 'hocomoco_v11'){
    # data[, 3:ncol(data) ] = data[, 3:ncol(data)-1 ]
    hocomoco_mof_annot = data.frame(X..motif_id = unique(data$X..motif_id), stringsAsFactors = F)
    hocomoco_mof_annot$Uniprot.name = sapply(strsplit(hocomoco_mof_annot$X..motif_id, split = '\\.'), head, 1)
    hocomoco_mof_annot$Uniprot.acc = uniprot_annot$Entry[ match(hocomoco_mof_annot$Uniprot.name, uniprot_annot$Entry.name) ]
    hocomoco_mof_annot$name = uniprot_annot$Gene.name.primary[ match(hocomoco_mof_annot$Uniprot.name, uniprot_annot$Entry.name) ]
    hocomoco_mof_annot$name[ is.na(hocomoco_mof_annot$name) ] = gsub('_HUMAN', '', hocomoco_mof_annot$Uniprot.name[ is.na(hocomoco_mof_annot$name) ]) 
    data$motif_alt_id = hocomoco_mof_annot$name[ match(data$X..motif_id, hocomoco_mof_annot$X..motif_id) ]
  }
  # Remove exact match duplicates
  data = unique(data[, ! colnames(data) %in% c('sequence_name', 'start', 'stop') ])
  # see global stats
  network = unique(data[, c('motif_alt_id', 'gene') ])
  names(network) = c('TF', 'target')
  for (i in c(100, 200, 500, 1000) ){
    message(rep('-', 10), ' ', i , ' ', rep('-', 10))
    # select all
    filename = paste(path, 'network_n', i, '.sif', sep = '')
    sdata = lapply(unique(data$X..motif_id), function(tf) head(subset(data, X..motif_id == tf ), i)  )
    save_network(sdata, filename)
    # select "open chromatin"
    filename = paste(path, 'network_n', i, '_opRegion.sif', sep = '')
    OP_idx = unique(c(grep('open_chromatin_region', data$SWEmbl_R0005_IDR), grep('ChIP_seq_region', data$SWEmbl_R0005_IDR)))
    # sdata = lapply(unique(data$X..motif_id[OP_idx]), function(tf) head(subset(data[OP_idx,], X..motif_id == tf ), i)  )
    sdata = lapply(unique(data$X..motif_id[OP_idx]), function(tf) data[ intersect(which(data$X..motif_id == tf)[1:i], OP_idx) , ]  )
    save_network(sdata, filename)
    # select conserved
    filename = paste(path, 'network_n', i, '_con.sif', sep = '')
    Con_idx = which( as.numeric(data$Phylop_scores) >= 3 | as.numeric(data$PhastCons_scores) >= 0.95 )
    # sdata = lapply(unique(data$X..motif_id[Con_idx]), function(tf) head(subset(data[Con_idx,], X..motif_id == tf ), i)  )
    sdata = lapply(unique(data$X..motif_id[Con_idx]), function(tf) data[ intersect(which(data$X..motif_id == tf)[1:i], Con_idx) , ]  )
    save_network(sdata, filename)
    # select "open chromatin" & conserved
    filename = paste(path, 'network_n', i, '_opRegion_con.sif', sep = '')
    both_idx = intersect(OP_idx, Con_idx) 
    # sdata = lapply(unique(data$X..motif_id[both_idx]), function(tf) head(subset(data[both_idx,], X..motif_id == tf ), i)  )
    sdata = lapply(unique(data$X..motif_id[both_idx]), function(tf) data[ intersect(which(data$X..motif_id == tf)[1:i], both_idx) , ]  )
    save_network(sdata, filename)
  }
}




