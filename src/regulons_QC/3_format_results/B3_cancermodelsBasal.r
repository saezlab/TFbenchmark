#' Merging tissue and non tissue specifi reults to generate the accurancy plots
#' corresponding to the B3 benchmarck datasets (Garcia-Alonso et al. Cancer Research 2018)
#' described in: https://www.biorxiv.org/content/biorxiv/early/2018/06/18/337915.full.pdf 
#'
#' This code merges all the information into a data frame 
#' that is the input to compute the accurancy values and plot the results.
#' 
#' NOTE: For B3, perturbed TFs are defined from independent essentilaity and CNA experiments on cancer cell lines
#' positive TFs (~perturbed) are defined as those showing cell line-specific essentiality (sd>4; essentiality score)
#' negative TFs are defined as non-essential (sd<-4; essentiality score) or CNA deleted TFs.




rm(list = ls())
home = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/'
setwd(home)
source('code/lib/functions.r')
source('code/lib/utils.r')
source('code/lib/data.r')




# Load and merge activity comparisons
df = merge_data2plot_accuracy(file = 'data/regulons_QC/B3_cell_lines/aggregated_activities.rdata')




## Add ACHILLES essentiality scores
essentiality = load_ACHILLES()
essentiality = t(scale(t(essentiality)))
essentiality = essentiality[ rownames(essentiality) %in% df$TF, ]
df$TF_essentiality_achilles = NA
idx = df$TF %in% rownames(essentiality)  & df$experiment %in% colnames(essentiality)
df$TF_essentiality_achilles[idx] = essentiality[ cbind(df$TF[idx], df$experiment[idx]) ]



## Add _DRIVEproject essentiality scores
essentiality = load_DRIVEproject()
essentiality = t(scale(t(essentiality)))
essentiality = essentiality[ rownames(essentiality) %in% df$TF, ]
df$TF_essentiality_DRIVEproject = NA
idx = df$TF %in% rownames(essentiality)  & df$experiment %in% colnames(essentiality)
df$TF_essentiality_DRIVEproject[idx] = essentiality[ cbind(df$TF[idx], df$experiment[idx]) ]



# Define essential genes
df$is_TF_essential = F
df$is_TF_essential[ which(df$TF_essentiality_achilles < -4 | df$TF_essentiality_DRIVEproject < -4 ) ] = T
df$is_TF_nonessential = F
df$is_TF_nonessential[  which(df$TF_essentiality_achilles > 4 | df$TF_essentiality_DRIVEproject > 4 ) ] = T



# Add homoDeletions
cell_lines_homoDel = load_homDel()
cell_lines_homoDel = cell_lines_homoDel[ rownames(cell_lines_homoDel) %in% df$TF, ]
df$is_TF_deleted = NA
idx = df$TF %in% rownames(cell_lines_homoDel)  & df$experiment %in% colnames(cell_lines_homoDel)
df$is_TF_deleted[idx] = cell_lines_homoDel[ cbind(df$TF[idx], df$experiment[idx]) ] == 1




# Define positive and negative samples
df$is_TF_inactive = F
df$is_TF_active = F
df$is_TF_inactive[ which(df$is_TF_deleted) ] = T
df$is_TF_inactive[ which(df$is_TF_nonessential) ] = T
df$is_TF_active[ which(df$is_TF_essential) ] = T

df$is_TF_perturbed = df$is_TF_active



# Filter TFs of interest: i.e. TFs in the active/inactive group
df$is_TF_of_interest = F
df$is_TF_of_interest [ df$is_TF_inactive | df$is_TF_active ] = T
df = subset(df, is_TF_of_interest)



# Invert rank
df$rank_nes = 1 - df$rank_nes
df$NES = df$NES * -(1)


save(df, file = 'data/regulons_QC/B3_cell_lines/aggregated_activities_formated.rdata' )
