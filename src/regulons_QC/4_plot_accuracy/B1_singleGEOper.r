#' Computing benchmark accuracy and plots
#' corresponding to the B1 benchmarck datasets 
#' described in: https://www.biorxiv.org/content/biorxiv/early/2018/06/18/337915.full.pdf 
#'
#' This code calculates the accuracy (PR AUC) for each regulon dataset
#' and generates the corresponding 
#' 1) Covergae vs accuracy plots
#' 2) Accuracy plots on the overlapping TFs




rm(list = ls())
home = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/'
setwd(home)
source('code/lib/utils.r')
source('code/lib/functions.r')
source('code/lib/plots.r')
source('code/regulons_QC/4_plot_accuracy/lib_accuracy.r')




# Load data 
load(file = 'data/regulons_QC/B1_perturbations/results/aggregated_activities_formated.rdata')


# Filter regulons to plot for publication
df = subset(df, ! regulon_dataset %in% c('cnio', 'TOP', 'TOP_PANCANCER') )
# group tissue specific regulons
df$regulon_dataset[ df$regulon_group == 'GTEx tissue-specific'] = 'GTEx tissue-specific'
df$regulon_dataset[ df$regulon_group == "GTEx tissue-specific noCOMBAT"] = "GTEx tissue-specific noCOMBAT"
df$regulon_dataset[ df$regulon_group == "cancer-specific"] = "cancer-specific"




# Define test and plot features
balanced_accuracy = F
plot_title = 'B1'
outdir = 'results/regulons_QC/B1_perturbations/'



# TFs
TFs = unique(df$TF[ df$is_TF_perturbed ])
num_TFs = length(unique(df$TF[ df$is_TF_perturbed ]))
num_TFs
write.table(TFs, file = 'data/regulons_QC/benchmark1_TFs.txt', quote = F, row.names = F, col.names = F)




#############################
# Coverage vs Accuracy plots
#############################


# Compute accuracy
accuracy = activities2accuracy(df, balanced_accuracy = balanced_accuracy)
save(accuracy, file = paste(outdir, 'regulons_accuracy_PRauc.rdata', sep = ''))


# Define plot range
x_scale = c(min(accuracy$AUC), signif(max(accuracy$AUC),1))
if( min(accuracy$AUC) > 0.475 )
  x_scale = c(0.47, max(accuracy$AUC)+0.01)
y_scale = c(0, max(accuracy$coverage))


# Plot individual resources
idx = Reduce(intersect, list(grep('_n100_', accuracy$regulon_dataset,  invert = T), 
                             grep('consensus', accuracy$regulon_evidence,  invert = T),
                             grep('cancer', accuracy$regulon_evidence,  invert = T),
                             grep('COMBAT', accuracy$regulon_dataset,  invert = T),
                             grep('omnipath_scores', accuracy$regulon_evidence,  invert = T)))
P = auc2coverage_plot(activities2accuracy_format(accuracy[idx,]), 
                      y_scale = y_scale, x_scale = x_scale, legend_position = 'bottom') + 
  ggtitle(plot_title)
P
ggsave(filename = paste(outdir, 'singe_evidences.png', sep = ''), width = 4.5, height = 5, dpi = 300)



# Plot combined resources
idx = grep('^consensus', accuracy$regulon_evidence,  invert = F)
P = auc2coverage_plot(activities2accuracy_format(accuracy[ idx , ]), 
                      y_scale = y_scale, x_scale = x_scale, legend_position = 'bottom') + 
  ggtitle(plot_title)
P
ggsave(filename = paste(outdir, 'combined_evidences.png', sep = ''), width = 4.5, height = 5, dpi = 300)



# Plot Omnipath scores
idx = grep('omnipath_scores', accuracy$regulon_evidence)
P = auc2coverage_plot(activities2accuracy_format(accuracy[ idx , ]), 
                      y_scale = y_scale, x_scale = x_scale, legend_position = 'bottom') + 
  ggtitle(plot_title)
P
ggsave(filename = paste(outdir, 'omnipathscores_evidences.png', sep = ''), width = 4.5, height = 5, dpi = 300)



# Plot GTEx vs TCGA specific inferred regulons
idx = setdiff(grep('inferred', accuracy$regulon_evidence), grep('noCOMBAT', accuracy$regulon_dataset) )
P = auc2coverage_plot_oncology(activities2accuracy_format(accuracy[ idx , ]),
                               y_scale = y_scale, x_scale = x_scale, legend_position = 'bottom') + 
  ggtitle(plot_title)
P
ggsave(filename = paste(outdir, 'combined_inferred_cancer.png', sep = ''), width = 4.5, height = 5, dpi = 300)
save(P, file = paste(outdir, 'combined_inferred_cancer.rdata', sep = ''))






##########################################################
# Accuracy boxplots on overlapping TFs
##########################################################
# Plot main lines of evidence
sdf = subset(df, regulon_dataset %in% 
               c('ReMap_n500', 'hocomoco_v11_n500', 'jaspar_2018_n200', 'GTEx_pantissue_i3', '1curateddatabases', '2curateddatabases') )
TFs_filter = Reduce(intersect, dlply(sdf, 'regulon_dataset',  function(x) unique(x$TF[x$is_TF_perturbed]) ))
line_colors = c(my_color_palette$EMBL[3], my_color_palette$EMBL[3], my_color_palette$EMBL[4], my_color_palette$EMBL[1], my_color_palette$EMBL[1], 'coral')
plot_PR(sdf, balanced_accuracy = balanced_accuracy, line_colors = line_colors, TFs_filter = TFs_filter) 
ggsave(filename = paste(outdir, 'overlapping_TFs.png', sep = ''), width = 5.5, height = 3.5, dpi = 300)



# Plot TFBS
TFs_filter = intersect(df$TF[grep('hocomoco', df$regulon_dataset)], df$TF[grep('jaspar', df$regulon_dataset)])
line_colors = c(brewer.pal(5, 'Blues')[2:5], brewer.pal(5, 'Reds')[2:5])
p = plot_rocs(df[grep('0$', df$regulon_dataset), ], 
              balanced_accuracy = balanced_accuracy,
              regulon_evidence_filter = 'TFBS_scanning',
              TFs_filter = TFs_filter, 
              line_colors =  my_color_palette$EMBL[1])
p + ggtitle(plot_title)
save(p, file = paste(outdir, 'ROC_TFBS_size.rdata', sep = ''))
p = plot_PR(df[grep('0$', df$regulon_dataset), ], 
            balanced_accuracy = balanced_accuracy,
            regulon_evidence_filter = 'TFBS_scanning',
            TFs_filter = TFs_filter, 
            line_colors =  my_color_palette$EMBL[1] )
p + ggtitle(plot_title)
save(p, file = paste(outdir, 'PR_TFBS_size.rdata', sep = ''))

p = plot_rocs(df[grep('200', df$regulon_dataset), ], 
              balanced_accuracy = balanced_accuracy,
              regulon_evidence_filter = 'TFBS_scanning',
              TFs_filter = TFs_filter, 
              line_colors =  my_color_palette$EMBL[1])
p + ggtitle(plot_title)
save(p, file = paste(outdir, 'ROC_TFBS_200filters.rdata', sep = ''))
p = plot_PR(df[grep('200', df$regulon_dataset), ], 
            balanced_accuracy = balanced_accuracy,
            regulon_evidence_filter = 'TFBS_scanning',
            TFs_filter = TFs_filter, 
            line_colors =  my_color_palette$EMBL[1])
p + ggtitle(plot_title)
save(p, file = paste(outdir, 'PR_TFBS_200filters.rdata', sep = ''))



# plot inferred
# GTEx only
sdf = df[ intersect(grep(df$regulon_dataset, pattern = '^GTEx'), 
                    grep(df$regulon_dataset, pattern = 'noCOMBAT', invert = T)) , ]
TFs_filter = Reduce(intersect, dlply(sdf, 'regulon_dataset',  function(x) unique(x$TF[x$is_TF_perturbed]) ))
p = plot_rocs(sdf,
              balanced_accuracy = balanced_accuracy,
              regulon_evidence_filter = 'inferred',
              TFs_filter = TFs_filter, 
              line_colors = my_color_palette$EMBL[4])
p + ggtitle(plot_title)
save(p, file = paste(outdir, 'ROC_inferred.rdata', sep = ''))
p = plot_PR(sdf,
              balanced_accuracy = balanced_accuracy,
              regulon_evidence_filter = 'inferred',
              TFs_filter = TFs_filter, 
              line_colors = my_color_palette$EMBL[4])
p + ggtitle(plot_title)
save(p, file = paste(outdir, 'PR_inferred.rdata', sep = ''))

# GTEx: COMBAT ves nonCOMBAT
p = plot_PR_compareCOMBAT(df[ grep('GTEx', df$regulon_dataset),],
                          balanced_accuracy = balanced_accuracy,
                          TFs_filter = TFs_filter, 
                          line_colors =  my_color_palette$EMBL[4])
p + ggtitle(plot_title)
p
ggsave(filename = paste(outdir, 'PR_inferred_ttest.png', sep = ''), width = 7, height = 7, dpi = 300)
save(p, file = paste(outdir, 'PR_inferred_ttest.rdata', sep = ''))

# GTEx vs cancer-TCGA
sdf = df[ grep(df$regulon_dataset, pattern = 'noCOMBAT', invert = T) , ]
sdf = sdf[ grep('inferred', sdf$regulon_evidence),]
TFs_filter = Reduce(intersect, dlply(sdf, 'regulon_dataset',  function(x) unique(x$TF) ))
p = plot_PR_compareCANCER(sdf,
                          balanced_accuracy = balanced_accuracy,
                          TFs_filter = TFs_filter, 
                          line_colors =  my_color_palette$EMBL[4])
p + ggtitle(plot_title)
save(p, file = paste(outdir, 'PR_inferred_cancer.rdata', sep = ''))




# Plot ChIP-Seq
sdf = df[ df$regulon_evidence == 'ChIP_Seq',]
TFs_filter = Reduce(intersect, dlply(sdf, 'regulon_dataset',  function(x) unique(x$TF) ))
p = plot_rocs(sdf, 
              TFs_filter = TFs_filter,
              balanced_accuracy = balanced_accuracy,
              line_colors =  'coral')
p + ggtitle(plot_title)
save(p, file = paste(outdir, 'ROC_chipSeq.rdata', sep = ''))
p = plot_PR(sdf, 
            TFs_filter = TFs_filter,
            balanced_accuracy = balanced_accuracy,
            line_colors =  'coral')
p + ggtitle(plot_title)
save(p, file = paste(outdir, 'PR_chipSeq.rdata', sep = ''))




# plot curated_databases
sdf = df[ df$regulon_evidence == 'curated_databases',]
TFs_filter = Reduce(intersect, dlply(sdf, 'regulon_dataset',  function(x) unique(x$TF) ))
p = plot_rocs(sdf,
              balanced_accuracy = balanced_accuracy,
              line_colors =  my_color_palette$EMBL[3])
p + ggtitle(plot_title)
save(p, file = paste(outdir, 'ROC_curated_databases.rdata', sep = ''))
p = plot_PR(sdf,
            balanced_accuracy = balanced_accuracy,
            line_colors =  my_color_palette$EMBL[3])
p + ggtitle(plot_title)
save(p, file = paste(outdir, 'PR_curated_databases.rdata', sep = ''))
