rm(list = ls())
home = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/'
setwd(home)
source('code/regulons_QC/3_format_results/lib_format_results.r')
source('code/regulons_QC/4_plot_benchmark_results/lib_plots_regulonsQC.r')


# Load and merge activity comparisons
df1 = build_mdf('data/regulons_QC/pertubations/GSE31534/activity_comparison_results.rdata')
df2 = build_mdf('data/regulons_QC/pertubations/GSE31912/activity_comparison_results.rdata')
df1$expreiment_GEOid = 'GSE31534'
df2$expreiment_GEOid = 'GSE31912'
df = rbind(df1, df2)
table(df$is_TF_perturbed)
table(is.na(df$is_TF_perturbed))




# add information
activators = load_activators()
repressors = load_repressors()
dual = load_dualActRep()
df$TF_type = 'unknown'
df$TF_type[ df$TF %in% activators ] = 'activators'
df$TF_type[ df$TF %in% repressors ] = 'repressors'
df$TF_type[ df$TF %in% dual ] = 'dual'





# Check low quality perturbations
low_quality_perturbations = c("MITF",
                              "ETV7", "FOXJ2",  "FOXP1" , "JUND", "ID1", "PAX3", "STAT5B",
                              "JUN") 
df$is_highquality_perturbation = ! df$perturbed_TF %in% low_quality_perturbations
p = plot_boxplots(subset(df, regulon_dataset == 'TOP'), wrap_variable = 'is_highquality_perturbation', p2_with_names = T)
p[[2]]



# Selecte TFs to benchmark
tfs_of_interest = setdiff(intersect(df$TF , df$perturbed_TF), low_quality_perturbations)
df = subset(df, perturbed_TF %in% tfs_of_interest )
save(df, file = 'data/regulons_QC/pertubations/GSE_activity_comparison_results_formated.rdata')


