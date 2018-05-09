rm(list = ls())
home = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/'
setwd(home)
source('/Volumes/GoogleDrive/My Drive/projects/pathway_activities/project_data_loads.r')
source('/Volumes/GoogleDrive/My Drive/projects/pathway_activities/DoRothEA/CODE/TF_scores_methods_v2/lib_enrichment_scores.r')
source('code/utils.r')



# Load data
E_cell_lines = t(scale(t(load_expression_gdsc(RNASEQ = T, norm_method = 'none'))))
E_pert_GSE31534 = get(load(file = 'data/regulons_QC/pertubations/GSE31534/GSE31534_A375_TFknockdown_zscores.rdata'))
E_pert_GSE31912 = get(load(file = 'data/regulons_QC/pertubations/GSE31912/GSE31912_MCF7_TFknockdown_zscores.rdata'))



# Set min gene set size
N = 4


# Load networks
networks = grep('sif', list.files('data/TF_target_sources/', recursive = T, pattern = 'network'), value=T)
networks = setdiff(networks, 'curated_databases/NFIRegulomeDB/network.sif')
networks = setdiff(networks, 'RANDOM_NETWORKS/3_network.sif')
for (n in networks){
  message(n)
  # Load network and fromat regulons
  net = read.delim(paste('data/TF_target_sources/', n, sep = ''), stringsAsFactors = F, sep = '\t')
  regulon = list()
  regulon$NAME = unique(net$TF)
  regulon$GENES = lapply(regulon$NAME, function(tf){ 
    targets = subset(net, TF == tf)$target
    x = rep(1, length(targets))
    names(x) = targets
    return(x)
    })
  outfile = n %>% gsub('/', '.', .) %>% gsub('network_', '', .) %>%  gsub('network', '', .) %>% gsub('sif', 'activities.rdata', .)  %>% gsub('\\.\\.', '.', .)
  # Cell lines activities
  activities = SLEA(E_cell_lines, genesets = regulon, method = 'VIPER', min_geneset_size = N, filter_E = F, max_shared = Inf)$NES
  save(activities, file = paste('data/TF_activities/cell_lines/', outfile, sep = ''))
  # Perturbation 1 activities
  activities = SLEA(E_pert_GSE31534, genesets = regulon, method = 'VIPER', min_geneset_size = N, filter_E = F, max_shared = Inf)$NES
  save(activities, file = paste('data/TF_activities/perturbations/GSE31534.', outfile, sep = ''))
  # Perturbation 2 activities
  activities = SLEA(E_pert_GSE31912, genesets = regulon, method = 'VIPER', min_geneset_size = N, filter_E = F, max_shared = Inf)$NES
  save(activities, file = paste('data/TF_activities/perturbations/GSE31912.', outfile, sep = ''))
}



