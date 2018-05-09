rm(list = ls())
home = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/'
setwd(home)
source('code/regulons_QC/3_format_results/lib_format_results.r')



# LOAD ##################################################################################################
#### UNSIGNED
# Load merged data
df = build_mdf(file = 'data/regulons_QC/pertubations/GEO/results/viper_aggregated_results_coexpressionPriors.rdata')
df_te = build_mdf(file = 'data/regulons_QC/pertubations/GEO/results/viper_aggregated_results_coexpressionPriors_tissuespecific.rdata')



# MERGE and FILTER ##################################################################################################
# merge
df = rbind(df, df_te)
df$is_TF_perturbed = df$TF == df$perturbed_TF
df$regulon_group = manage_datasets_names(df$regulon_dataset)


# PROCESS ##################################################################################################
# add experiment information
experiments_info = read.csv('data/regulons_QC/pertubations/GEO/contrasts/TF_perturbation_confidence.csv', stringsAsFactors = F)
experiments_info$perturbation_GEOid = strsplit(as.character(experiments_info$id), '\\.') %>% sapply(., tail, 1) %>% strsplit(., '-') %>%  sapply(., head, 1)
df$perturbation_effect = experiments_info$perturbation_effect[ match(df$experiment, experiments_info$id) ]
df$perturbation_treatment = experiments_info$perturbation_treatment[ match(df$experiment, experiments_info$id) ]
df$perturbation_GEOid = experiments_info$perturbation_GEOid[ match(df$experiment, experiments_info$id) ]
df$is_TF_DE_opposite = experiments_info$is_TF_DE_opposite[ match(df$experiment, experiments_info$id) ]
df$is_TF_DE_expected = experiments_info$is_TF_DE_expected[ match(df$experiment, experiments_info$id) ]



# Selecte TFs to benchmark
tfs_of_interest =  intersect(df$TF , df$perturbed_TF)
df = subset(df, TF %in% c(tfs_of_interest, 'random'))
df = subset(df, ! (perturbation_treatment == 'overexpression' & is_TF_DE_opposite) )
length(tfs_of_interest)
length(unique(df$experiment))




# add extra information
activators = load_activators()
repressors = load_repressors()
dual = load_dualActRep()
df$TF_type = 'unknown'
df$TF_type[ df$TF %in% activators ] = 'activators'
df$TF_type[ df$TF %in% repressors ] = 'repressors'
df$TF_type[ df$TF %in% dual ] = 'dual'



# remove is_TF_DE_opposite
ggplot(subset(df, is_TF_perturbed &  regulon_evidence == 'old_consensus' ), aes(x = is_TF_DE_opposite, y = rank_nes) ) + geom_boxplot()
df = subset(df, ! is_TF_DE_opposite)


# save
save(df, file = 'data/regulons_QC/pertubations/GEO/results/viper_aggregated_results_coexpressionPriors_formated.rdata')


