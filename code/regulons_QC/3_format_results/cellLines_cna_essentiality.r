rm(list = ls())
home = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/'
setwd(home)
source('code/regulons_QC/3_format_results/lib_format_results.r')




# Load and merge activity comparisons
df = build_mdf(file = 'data/regulons_QC/cell_lines/activity_comparison_results.rdata')



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
# df$is_TF_essential[ which(df$TF_essentiality_achilles <= -4 ) ] = T
df$is_TF_nonessential = F
df$is_TF_nonessential[  which(df$TF_essentiality_achilles > 4 | df$TF_essentiality_DRIVEproject > 4 ) ] = T
# df$is_TF_nonessential[  which(df$TF_essentiality_achilles >= 4 ) ] = T
table(df$is_TF_essential)
table(df$is_TF_nonessential)





# Add homoDeletions
cell_lines_homoDel = load_homDel()
cell_lines_homoDel = cell_lines_homoDel[ rownames(cell_lines_homoDel) %in% df$TF, ]
# cell_lines_homoDel[ rowSums(cell_lines_homoDel) < 500, ]
df$is_TF_deleted = NA
idx = df$TF %in% rownames(cell_lines_homoDel)  & df$experiment %in% colnames(cell_lines_homoDel)
df$is_TF_deleted[idx] = cell_lines_homoDel[ cbind(df$TF[idx], df$experiment[idx]) ] == 1
table(df$is_TF_deleted)



# Define positive and negative samples
df$is_TF_inactive = F
df$is_TF_inactive[ which(df$is_TF_deleted) ] = T
df$is_TF_inactive[ which(df$is_TF_nonessential) ] = T
table(df$is_TF_inactive)

df$is_TF_active = F
df$is_TF_active[ which(df$is_TF_essential) ] = T
table(df$is_TF_active)


df$is_TF_perturbed = df$is_TF_active
df$is_TF_of_interest = F
df$is_TF_of_interest [ df$is_TF_inactive | df$is_TF_active ] = T


df = subset(df, is_TF_of_interest)
save(df, file = 'data/regulons_QC/cell_lines/activity_comparison_results_formated.rdata' )