rm(list = ls())
home = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/'
setwd(home)
source('code/utils.r')
source('code/TFperturbations/contrast_lib.r')
library(viper)




# Set minimum regulon size
Nmax = 10000
Nmin = 4



# find perturbed TFs
perturbed_TFs = list.dirs(path = 'data/regulons_QC/pertubations/GEO/design/', full.names = F) %>% setdiff(., '')
message('\n\n\n\n\nLoading regulons for ', length(perturbed_TFs), ' Transcription Factors')



# Load REAL networks + filter TFs_perturbed
network_files = list.files('data/TF_target_sources', recursive = T, pattern = 'viperRegulon.rdata')
network_files = grep("RANDOM_NETWORKS", network_files, invert = T, value = T)
network_files = grep("aracne", network_files, invert = T, value = T)
networks = list()
for (n in network_files){
  # Load network and fromat regulons
  net = get(load(paste('data/TF_target_sources/', n, sep = '')))
  if( length(net) < 2 )
    next()
  net_tfs = sapply(strsplit(names(net), split = ' '), head, 1)
  net_tfs = sapply(strsplit(net_tfs, split = '_'), head, 1)
  net = net[ net_tfs %in% perturbed_TFs ]
  if( length(net) > 1 ){
    networks = append(networks, net)
    message(n)
  }
}




# Load perturbations and run GSEA
perturbation_files = list.files(path = 'data/regulons_QC/pertubations/GEO/contrasts', pattern = 'rdata', full.names = T)
for (myFile in rev(perturbation_files)){
  message('Loading file: ', myFile)
  DEGs = get(load(myFile))
  DEGs = subset(DEGs, Symbol != "" )
  DEGs = subset(DEGs, ! duplicated(Symbol))
  
  message('Reading design file to extract tissue')
  perturbation_id = tail(unlist(strsplit(myFile, split = '/')),1) %>% gsub('.rdata', '', .) %>% strsplit(., split = '\\.') %>% unlist(.)
  design_file = list.files('data/regulons_QC/pertubations/GEO/design', recursive = T, full.names = T) %>% grep(perturbation_id[1], ., value = T) %>% grep(perturbation_id[2], ., value = T)
  tissue = read_desigfile(design_file)$GTEx_tissue
  
  message('Formating data for analysis')
  myStatistics = matrix(DEGs$logFC, dimnames = list(DEGs$Symbol, 'logFC') )
  myPvalue = matrix(DEGs$P.Value, dimnames = list(DEGs$Symbol, 'P.Value') )
  # Although Gene Expression Signature (GES) can be defined by the t-statistic, to be consistent with the z-score based null model for msVIPER (see section 6.2 in its manual), 
  # we will estimate z-score values for the GES as indicated in the manual:
  mySignature = (qnorm(myPvalue/2, lower.tail = FALSE) * sign(myStatistics))[, 1]
  mySignature = mySignature[order(mySignature, decreasing = T)]

  message('Running VIPER ...')
  mrs = msviper(mySignature, networks, verbose = T, minsize = Nmin, ges.filter = F)
  results = data.frame(Regulon = names(mrs$es$nes),
                       Size = mrs$es$size[ names(mrs$es$nes) ], 
                       Size_adaptative = sapply(mrs$regulon[ names(mrs$es$nes) ], function(x) sum((x$likelihood/max(x$likelihood))^2) ), 
                       NES = mrs$es$nes, 
                       p.value = mrs$es$p.value, 
                       FDR = p.adjust(mrs$es$p.value, method = 'fdr'),
                       perturbation_tissue = tissue)
  results = results[ order(results$p.value, decreasing = F) , ]
  
  analysis_name = gsub('contrasts/', 'activities/real/viper_', myFile)
  save(results, file = analysis_name)
  message('\nDone!\n\n')
  message('')
}





