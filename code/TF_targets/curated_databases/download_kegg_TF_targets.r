rm(list = ls())
home = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/'
setwd(home)
library(KEGGgraph)
library(KEGGREST)
library(KEGG.db)
library(org.Hs.eg.db)
library(reshape2)
options(error = recover)
source('code/utils.r')


get_GErel_edges = function(path_name, path_id){
  # path_name = pName[116]
  # path_id = pId[116]
  
  # path_name = pName[i]
  # path_id = pId[i]
  message(path_name)
  tmp = paste(tempfile(tmpdir = '~/tmp/'), '.xml', sep = '')
  my_kgml = retrieveKGML(pathwayid = path_id, organism="hsa", destfile = tmp)
  mapkG = try(parseKGML2Graph(tmp, expandGenes=T), silent = T)
  if(class(mapkG) == "try-error")
    return(NULL)    
  pathedges =  getKEGGedgeData(mapkG)
  pathedges_type = sapply(pathedges, getType)
  pathedges_GErel = pathedges[ pathedges_type == 'GErel' ]
  pathedges_entryID = sapply(pathedges_GErel, getEntryID)
  pathedges_subtype_name = sapply(pathedges_GErel, function(x) getName(getSubtype(x)[[1]]) )
  pathedges_subtype_value = sapply(pathedges_GErel, function(x) getValue(getSubtype(x)[[1]]) )
  if( is.list(pathedges_entryID) )
    return(NULL)
  df = as.data.frame(cbind(t(pathedges_entryID), pathedges_subtype_name, pathedges_subtype_value, path_name = path_name, path_id = path_id), stringsAsFactors = F)
  df$Entry1Name = sapply(mget(translateKEGGID2GeneID(df$Entry1ID), org.Hs.egSYMBOL, ifnotfound=NA), "[[",1)
  df$Entry2Name = sapply(mget(translateKEGGID2GeneID(df$Entry2ID), org.Hs.egSYMBOL, ifnotfound=NA), "[[",1)
  return(df)
}



mykegglist = keggList(organism = "hsa", database = 'pathway')
pId = gsub('path:hsa', '', names(mykegglist))
pName = mykegglist
GErel_edges = mapply(get_GErel_edges, pName, pId)
xGErel_edges = Filter(Negate(is.null), GErel_edges)
mGErel_edges = melt(xGErel_edges, id.vars = names(xGErel_edges[[1]]))
mGErel_edges = mGErel_edges[ order(mGErel_edges$Entry2Name) , ]
mGErel_edges = mGErel_edges[ order(mGErel_edges$Entry1Name) , ]
write.table(mGErel_edges, file = 'data/TF_target_sources/curated_databases/kegg/GErel_hsa.txt', col.names = T, row.names = F, quote = F, sep = '\t')



mGErel_edges = read.delim(file = 'data/TF_target_sources/curated_databases/kegg/GErel_hsa.txt', header = T, stringsAsFactors = F)
network = unique(mGErel_edges[, 7:8])
names(network)[1:2] = c('TF', 'target')
network = network[ network$TF %in% load_TFs_census(), ]
network = network[ order(network$TF), ]
summarize_network(network)
write.table(network, file = 'data/TF_target_sources/curated_databases/kegg/network.sif', sep = '\t', row.names = F, quote = F)
