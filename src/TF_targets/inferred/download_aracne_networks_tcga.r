rm(list = ls())
home = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/'
setwd(home)
source('code/lib/utils.r')
source('code/lib/data.r')
source('code/TF_targets/inferred/lib_download_aracne_networks_tcga.r')



ensemblgenes_annot = load_hugo_entrez_annotation()
TFcensus = load_TFs_census()



require(aracne.networks)
for (regulon_name in data(package="aracne.networks")$results[, "Item"]){
  print(regulon_name)
  regulon = retieve_aracne_network(regulon_name)
  reg_name = gsub('regulon', '', regulon_name)
  names(regulon) = paste(names(regulon), ' tcga_', reg_name, sep = '' )
  save(regulon, file = paste('data/TF_target_sources/inferred/TCGA/cancer_specific/', reg_name, '_viperRegulon.rdata' , sep = '') )
}
