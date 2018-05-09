rm(list = ls())
home = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/'
setwd(home)
source('code/utils.r')
library(fgsea)




perturbation_files = list.files(path = 'data/regulons_QC/pertubations/GEO/contrasts', pattern = 'rdata', full.names = T)
TFs =  unique(sapply(strsplit(sapply(strsplit(perturbation_files, split = "/"), tail, 1), split = "\\."), head, 1) )


# Load networks
network_files = grep('sif', list.files('data/TF_target_sources/', recursive = T, pattern = 'network'), value=T)
network_files = setdiff(network_files, "reverse_engenieered/network.sif")
NETWORKS = list()
for (n in network_files){
  message(n)
  # Load network and fromat regulons
  net = read.delim(paste('data/TF_target_sources/', n, sep = ''), stringsAsFactors = F, sep = '\t')
  net = subset(net, TF %in% TFs)
  if( ! is.null(net$effect) )
    net = subset(net, effect != -1 )
  NETWORKS[[n]] = net
}

# Define background targets
targets = unlist(sapply(NETWORKS, function(net) net$target))
u_targets = sort(unique(targets))
f_targets = table(targets)


# Format regulons as gene-sets
networks = list()
for (n in network_files){
  message(n)
  net = NETWORKS[[n]]
  for (tf in unique(net$TF) ){
    tf_name = paste(tf, n, sep = ' - ')
    networks[[tf_name]] = unique(subset(net, TF == tf & target != "")$target)
    # random network
    for (ri in 1:10){
      si = length(networks[[tf_name]])
      ntargets = sample(u_targets, replace = F, prob = f_targets, size = si)
      networks[[paste(tf , ' - random_balanced ', ri, ' - ', n, sep = '')]] = ntargets
      ntargets = sample(u_targets, replace = F, size = si)
      networks[[paste(tf , ' - random_unbalanced ', ri, ' - ', n, sep = '')]] = ntargets
    }
  }
}




# Load perturbations and run GSEA
for (myFile in perturbation_files){
  message('Loading file: ', myFile)
  DEGs = get(load(myFile))
  DEGs = subset(DEGs, Symbol != "" )
  message('Formating data for analysis')
  myRanks = DEGs$logFC
  names(myRanks) = DEGs$Symbol
  myRanks = myRanks[ order(myRanks) ]
  
  message('------ GSEA ------')
  results = fgsea(pathways = networks,  stats = myRanks, minSize=4, maxSize=2000, nperm=1000)
  results = results[ order(results$NES, decreasing = T) , ]
  
  analysis_name = gsub('contrasts/', 'gsea/fgsea_', myFile)
  save(results, file = analysis_name)
  message('-----------------')
  message('')
}





