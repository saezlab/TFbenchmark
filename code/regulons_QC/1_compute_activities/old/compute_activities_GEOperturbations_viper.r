rm(list = ls())
home = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/'
setwd(home)
source('code/utils.r')
library(viper)





# List perturbation files and extract perturbed TFd
perturbation_files = list.files(path = 'data/regulons_QC/pertubations/GEO/contrasts', pattern = 'rdata', full.names = T)
TFs_perturbed = sapply(strsplit(perturbation_files, split = 'contrasts/'), tail, 1)
TFs_perturbed = unique(sapply(strsplit(TFs_perturbed, split = '\\.'), head, 1))


# Load TF census and repressors
TFs =  unique(sapply(strsplit(sapply(strsplit(perturbation_files, split = "/"), tail, 1), split = "\\."), head, 1) )
repressors = load_repressors()



# Set minimum regulon size
N = 4



# Load REAL networks + filter TFs_perturbed + format for viper input
network_files = list.files('data/TF_target_sources/', recursive = T, pattern = 'network') %>% grep('sif', ., value=T)  %>% grep('RANDOM', ., value=T, invert= T)
network_files = setdiff(network_files, "reverse_engenieered/network.sif")
networks = list()
for (n in network_files){
  message(n)
  # Load network and fromat regulons
  net = read.delim(paste('data/TF_target_sources/', n, sep = ''), stringsAsFactors = F, sep = '\t')
  net = subset(net, TF %in% TFs_perturbed)
  # Remove self regulation
  net = subset(net, TF != target)
  # Remove empty targets
  net = subset(net, target != "")
  tfs_of_interest = names(which(table(net$TF) >= N))
  for (tf in tfs_of_interest ){
    tf_name = paste(tf, n, sep = ' - ')
    tf_targets = subset(net, TF == tf)$target
    if( ! is.null(net$effect) ){
      tf_targets_sign = subset(net, TF == tf)$effect
      tf_targets_sign = tf_targets_sign[ tf_targets != tf ]
      tf_targets = tf_targets[ tf_targets != tf ]
      tf_targets_sign[ tf_targets_sign == 0 ] = 1
    }else{
      tf_targets = tf_targets[ tf_targets != tf ]
      tf_targets_sign = rep(1, length(tf_targets))
      if( tf %in% repressors )
        tf_targets_sign = -tf_targets_sign
    }
    if( length(tf_targets_sign) == 0 )
      next()
    si = length(tf_targets)
    # build real regulon
    networks[[tf_name]]$tfmode = tf_targets_sign # NOTE: ttmode has to be the first element of the list. Otherwise msviper will fail.
    names(networks[[tf_name]]$tfmode) = tf_targets
    networks[[tf_name]]$likelihood = rep(1, length(tf_targets))
  }
}


# Load RANDOM networks + filter those with the sizes of interest + format for viper input
network_files = list.files('data/TF_target_sources/RANDOM_NETWORKS', recursive = T, pattern = 'network', full.names = T)
for (n in network_files){
  message(n)
  # Load network and fromat regulons
  net = read.delim(n, stringsAsFactors = F, sep = '\t')
  for (tf in unique(net$TF) ){
    tf_name = tf
    tf_targets = subset(net, TF == tf)$target
    tf_targets = tf_targets[ tf_targets != tf ]
    tf_targets_sign = rep(1, length(tf_targets))
    if( tf %in% repressors )
      tf_targets_sign = -tf_targets_sign
    si = length(tf_targets)
    # build real regulon
    networks[[tf_name]]$tfmode = tf_targets_sign # NOTE: ttmode has to be the first element of the list. Otherwise msviper will fail.
    names(networks[[tf_name]]$tfmode) = tf_targets
    networks[[tf_name]]$likelihood = rep(1, length(tf_targets))
  }
}



# Load perturbations and run GSEA
for (myFile in perturbation_files){
  message('Loading file: ', myFile)
  DEGs = get(load(myFile))
  DEGs = subset(DEGs, Symbol != "" )
  
  message('Formating data for analysis')
  myStatistics = matrix(DEGs$logFC, dimnames = list(DEGs$Symbol, 'logFC') )
  myPvalue = matrix(DEGs$P.Value, dimnames = list(DEGs$Symbol, 'P.Value') )
  # Although Gene Expression Signature (GES) can be defined by the t-statistic, to be consistent with the z-score based null model for msVIPER (see section 6.2 in its manual), 
  # we will estimate z-score values for the GES as indicated in the manual:
  mySignature = (qnorm(myPvalue/2, lower.tail = FALSE) * sign(myStatistics))[, 1]
  mySignature = mySignature[order(mySignature, decreasing = T)]

  message('------ VIPER ------')
  mrs = msviper(mySignature, networks, verbose = T, minsize = N, pleiotropy = F, ges.filter = F)
  results = data.frame(Regulon = names(mrs$es$nes),
                       Size = mrs$es$size[ names(mrs$es$nes) ], 
                       NES = mrs$es$nes, 
                       p.value = mrs$es$p.value, 
                       FDR = p.adjust(mrs$es$p.value, method = 'fdr') )
  results = results[ order(results$NES, decreasing = T) , ]
  
  analysis_name = gsub('contrasts/', 'gsea/viper_', myFile)
  save(results, file = analysis_name)
  message('\n-----------------')
  message('')
}





