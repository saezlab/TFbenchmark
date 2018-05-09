rm(list = ls())
home = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/'
setwd(home)
source('code/applications/drug_response/associations_drugresponse_lib.r')
source('code/utils.r')




#---------------------- TF activities
message('Loading activities')
activities = get(load('data/regulons_QC/cell_lines/activities/omnipath_scores.TOP.network.rdata'))
rownames(activities) = sapply(strsplit(rownames(activities), split = ' - '), head, 1)

# Remove low confidence regulons
activities = activities[ grep('_E$', rownames(activities), invert = T) , ]
activities = activities[ grep('_D$', rownames(activities), invert = T) , ]
# activities = activities[ grep('_C$', rownames(activities), invert = T) , ]
# activities = activities[ grep('_B$', rownames(activities), invert = T) , ]
# activities = activities[ grep('_A$', rownames(activities), invert = T) , ]
dim(activities)
#----------------------




# # #---------------------- Pancancer ASSOS
# associations_df = analyse_TF_drug_associations(tf_activities = activities, tissue = 'pancancer')
# plot_volcano(associations_df, th = 0.01, main = 'TF-gene drug interactions\n(pancancer)')
# ggsave(filename = 'results/applications/drugresponse_assos/pancancer.png', dpi = 300, width = 6, height = 6)
# save(associations_df, file = 'results/applications/drugresponse_assos/pancancer.rdata')
# # #----------------------




#---------------------- Cancer-type ASSOS
tissues_of_interest = names(which(table(covariates.annot$gdsc_desc_2)  > 10) )
for (tissue in rev(sort(tissues_of_interest))){
  message(tissue)
  associations_df = analyse_TF_drug_associations(activities, tissue = tissue)
  plot_volcano(associations_df)
  ggsave(filename = paste('results/applications/drugresponse_assos/', tissue, '.png', sep = ''), dpi = 300, width = 10, height = 10)
  save(associations_df, file = paste('results/applications/drugresponse_assos/', tissue, '.rdata', sep = ''))
}
#----------------------

