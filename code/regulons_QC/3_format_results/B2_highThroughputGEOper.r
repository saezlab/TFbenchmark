#' Merging tissue and non tissue specifi reults to generate the accurancy plots
#' corresponding to the B2 benchmarck datasets
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
# source('code/regulons_QC/4_plot_accuracy/lib_plots_regulonsQC.r')



# Prepare data frame to plot accuracy
df1 = merge_data2plot_accuracy('data/regulons_QC/B2_perturbations/GSE31534/activity_comparison_results.rdata')
df2 = merge_data2plot_accuracy('data/regulons_QC/B2_perturbations/GSE31912/activity_comparison_results.rdata')
df1$expreiment_GEOid = 'GSE31534'
df2$expreiment_GEOid = 'GSE31912'
df = rbind(df1, df2)





# # add information
# activators = load_activators()
# repressors = load_repressors()
# dual = load_dualActRep()
# df$TF_type = 'unknown'
# df$TF_type[ df$TF %in% activators ] = 'activators'
# df$TF_type[ df$TF %in% repressors ] = 'repressors'
# df$TF_type[ df$TF %in% dual ] = 'dual'




# Remove low quality perturbations: i.e. gene not downregulated after knock-down 
low_quality_perturbations = c("MITF",
                              "ETV7", "FOXJ2",  "FOXP1" , "JUND", "ID1", "PAX3", "STAT5B",
                              "JUN") 
df$is_highquality_perturbation = ! df$perturbed_TF %in% low_quality_perturbations
# Selecte TFs to benchmark
tfs_of_interest = setdiff(intersect(df$TF , df$perturbed_TF), low_quality_perturbations)
df = subset(df, perturbed_TF %in% tfs_of_interest )


# Save
save(df, file = 'data/regulons_QC/B2_perturbations/aggregated_activities_formated.rdata')


