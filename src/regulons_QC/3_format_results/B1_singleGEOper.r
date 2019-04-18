#' Merging tissue and non tissue specifi reults to generate the accurancy plots
#' corresponding to the B1 benchmarck datasets 
#' described in: https://www.biorxiv.org/content/biorxiv/early/2018/06/18/337915.full.pdf 
#'
#' This code merges all the information into a data frame 
#' that is the input to compute the accurancy values and plot the results.







rm(list = ls())
home = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/'
setwd(home)
source('code/lib/functions.r')
source('code/lib/utils.r')
# source('code/lib/data.r')



# Prepare data frame to plot accuracy
df = merge_data2plot_accuracy(file = 'data/regulons_QC/B1_perturbations/results/aggregated_activities.rdata')

# Define perturbed TFs
df$is_TF_perturbed = df$TF == df$perturbed_TF


# add experiment information
experiments_info = read.csv('data/regulons_QC/B1_perturbations/contrasts/TF_perturbation_confidence.csv', stringsAsFactors = F)
experiments_info$perturbation_GEOid = strsplit(as.character(experiments_info$id), '\\.') %>% sapply(., tail, 1) %>% strsplit(., '-') %>%  sapply(., head, 1)
df$perturbation_effect = experiments_info$perturbation_effect[ match(df$experiment, experiments_info$id) ]
df$perturbation_treatment = experiments_info$perturbation_treatment[ match(df$experiment, experiments_info$id) ]
df$perturbation_GEOid = experiments_info$perturbation_GEOid[ match(df$experiment, experiments_info$id) ]
# NOTE: we checked if the TF expression changed according the perturbation (i.e. higher expression after overexpression & lower expression after inhibition)
df$is_TF_DE_opposite = experiments_info$is_TF_DE_opposite[ match(df$experiment, experiments_info$id) ] # means the TF is differentially expressed in the opposite direction as expected
df$is_TF_DE_expected = experiments_info$is_TF_DE_expected[ match(df$experiment, experiments_info$id) ] # means the TF is differentially expressed as expected



# # add TF's extra information
# activators = load_activators()
# repressors = load_repressors()
# dual = load_dualActRep()
# df$TF_type = 'unknown'
# df$TF_type[ df$TF %in% activators ] = 'activators'
# df$TF_type[ df$TF %in% repressors ] = 'repressors'
# df$TF_type[ df$TF %in% dual ] = 'dual'


# Selecte TFs to benchmark
df = subset(df, ! is_TF_DE_opposite)
df = subset(df, is_TF_DE_expected | ! perturbation_treatment %in% c("shRNA", "siRNA", "overexpression" ) )
tfs_of_interest =  intersect(df$TF , df$perturbed_TF)
df = subset(df, TF %in% tfs_of_interest & perturbed_TF %in% tfs_of_interest) # Filter TFs: only test TFs that are perturbed at least in one experiment



# Invert rank
df$rank_nes = 1 - df$rank_nes
df$NES = df$NES * -(1)


# save
save(df, file = 'data/regulons_QC/B1_perturbations/results/aggregated_activities_formated.rdata')


