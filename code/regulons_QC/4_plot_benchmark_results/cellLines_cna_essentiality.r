rm(list = ls())
home = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/'
setwd(home)
source('code/regulons_QC/4_plot_benchmark_results/lib_plots_regulonsQC.r')




# Load and merge activity comparisons
load(file = 'data/regulons_QC/cell_lines/activity_comparison_results_formated.rdata' )
df = subset(df, ! regulon_dataset %in% c('all_weighted', 'cnio', 'inferredGTEx_dynamic', 'evidence4', 'TOP') )
df$rank_nes = 1 - df$rank_nes
df$NES = df$NES * -(1)


# Counts
nrow(unique(df[ df$is_TF_active, c('experiment', 'TF') ]))
num_TFs = length(unique(df$TF[  df$is_TF_active ]))
num_TFs


# Define test and plot features
balanced_accuracy = F
plot_title = 'B3'
outdir = 'results/regulons_QC/cell_lines/'






# PLOTS MAIN ##################################################################################################
x_scale = c(0.47, 1)
y_scale = c(0, num_TFs)



idx = Reduce(intersect, list(grep('_n100_', df$regulon_dataset,  invert = T), 
                             grep('consensus', df$regulon_evidence,  invert = T),
                             grep('omnipath_scores', df$regulon_evidence,  invert = T)))
P = auc2coverage_plot(df[ idx , ], balanced_accuracy = balanced_accuracy, y_scale = y_scale, x_scale = x_scale, legend_position = 'bottom')
P$aucs
P$P + ggtitle(plot_title)
ggsave(filename = paste(outdir, 'singe_evidences.png', sep = ''), width = 4.5, height = 5, dpi = 300)


idx = grep('^consensus', df$regulon_evidence,  invert = F)
P = auc2coverage_plot(df[ idx , ], balanced_accuracy = balanced_accuracy, y_scale = y_scale, x_scale = x_scale, legend_position = 'bottom')
P$aucs
P$P + ggtitle(plot_title)
ggsave(filename = paste(outdir, 'combined_evidences.png', sep = ''), width = 4.5, height = 5, dpi = 300)

idx = grep('evidence2_', df$regulon_dataset,  invert = F)
P = auc2coverage_plot(df[ idx , ], balanced_accuracy = balanced_accuracy, y_scale = y_scale, x_scale = x_scale, filter_datasets = F, legend_position = 'bottom')
P$aucs
P$P + ggtitle(plot_title)
ggsave(filename = paste(outdir, 'combined_evidences_signed.png', sep = ''), width = 4.5, height = 5, dpi = 300)






# plot OMNIPATH SCORES
TFs_filter = Reduce(intersect, lapply(LETTERS[1], function(let) subset(df, regulon_dataset == let)$TF ) )
line_colors = c(brewer.pal(5, 'Paired')[c(2,1,4,3,5)], 'gray20', 'purple')
plot_rocs(df, 
          balanced_accuracy = balanced_accuracy,
          regulon_evidence_filter = 'omnipath_scores',
          # TFs_filter = TFs_filter,
          line_colors = line_colors )
ggsave(filename = paste(outdir, 'ROC_omnipath_scores.png', sep = ''), width = 4, height = 4, dpi = 300)
plot_PR(df, 
        balanced_accuracy = balanced_accuracy,
        regulon_evidence_filter = 'omnipath_scores',
        # TFs_filter = TFs_filter,
        line_colors = line_colors)
ggsave(filename = paste(outdir, 'PR_omnipath_scores.png', sep = ''), width = 7, height = 4, dpi = 300)


idx = grep('omnipath_scores', df$regulon_evidence)
P = auc2coverage_plot(df[ idx , ], balanced_accuracy = balanced_accuracy, y_scale = y_scale, x_scale = x_scale, legend_position = 'bottom')
P$aucs
P$P
ggsave(filename = paste(outdir, 'omnipathscores_evidences.png', sep = ''), width = 4, height = 4, dpi = 300)



# plot TFBS
TFs_filter = intersect(df$TF[grep('hocomoco', df$regulon_dataset)], df$TF[grep('jaspar', df$regulon_dataset)])
line_colors = c(brewer.pal(5, 'Blues')[2:5], brewer.pal(5, 'Reds')[2:5])
p = plot_rocs(df[grep('0$', df$regulon_dataset), ], 
              balanced_accuracy = balanced_accuracy,
              regulon_evidence_filter = 'TFBS_scanning',
              TFs_filter = TFs_filter, 
              line_colors = line_colors)
p
ggsave(filename = paste(outdir, 'ROC_TFBS_size.png', sep = ''), width = 4, height = 4, dpi = 300)
save(p, file = paste(outdir, 'ROC_TFBS_size.rdata', sep = ''))
p = plot_PR(df[grep('0$', df$regulon_dataset), ], 
            balanced_accuracy = balanced_accuracy,
            regulon_evidence_filter = 'TFBS_scanning',
            TFs_filter = TFs_filter, 
            line_colors =  my_color_palette$EMBL[1] )
p
ggsave(filename = paste(outdir, 'PR_TFBS_size.png', sep = ''), width = 7, height = 4, dpi = 300)
save(p, file = paste(outdir, 'PR_TFBS_size.rdata', sep = ''))

p = plot_rocs(df[grep('200', df$regulon_dataset), ], 
              balanced_accuracy = balanced_accuracy,
              regulon_evidence_filter = 'TFBS_scanning',
              TFs_filter = TFs_filter, 
              line_colors = rev(line_colors)[ c(5:8,1:4) ])
p
ggsave(filename = paste(outdir, 'ROC_TFBS_200filters.png', sep = ''), width = 4, height = 4, dpi = 300)
save(p, file = paste(outdir, 'ROC_TFBS_200filters.rdata', sep = ''))
p = plot_PR(df[grep('200', df$regulon_dataset), ], 
            balanced_accuracy = balanced_accuracy,
            regulon_evidence_filter = 'TFBS_scanning',
            TFs_filter = TFs_filter, 
            line_colors =  my_color_palette$EMBL[1])
p
ggsave(filename = paste(outdir, 'PR_TFBS_200filters.png', sep = ''), width = 7, height = 4, dpi = 300)
save(p, file = paste(outdir, 'PR_TFBS_200filters.rdata', sep = ''))


idx = Reduce(intersect, list(c(setdiff(grep('TFBS', df$regulon_evidence), grep('_n100_', df$regulon_dataset)),
                               grep('TFBS', df$regulon_dataset))))
P = auc2coverage_plot(df[ idx , ], balanced_accuracy = balanced_accuracy, y_scale = y_scale, x_scale = x_scale)
P$aucs
P$P + ggtitle(plot_title)
ggsave(filename = paste(outdir, 'combined_TFBS.png', sep = ''), width = 4, height = 4, dpi = 300)







# plot inferredGTEx
inferredGTEx_common_TFs = Reduce(intersect, list(df$TF[intersect(grep('inferredGTEx_2', df$regulon_dataset), which(df$is_TF_perturbed))],
                                                 df$TF[intersect(grep('inferredGTEx_3', df$regulon_dataset), which(df$is_TF_perturbed))],
                                                 df$TF[intersect(grep('inferredGTEx_5', df$regulon_dataset), which(df$is_TF_perturbed))],
                                                 df$TF[intersect(grep('inferredGTEx_10', df$regulon_dataset), which(df$is_TF_perturbed))]) )
p = plot_rocs(df[ grep('inferredGTEx_', df$regulon_dataset) ,],
              balanced_accuracy = balanced_accuracy,
              regulon_evidence_filter = 'inferredGTEx',
              TFs_filter = inferredGTEx_common_TFs, 
              line_colors = my_color_palette$shiny)
p
ggsave(filename = paste(outdir, 'ROC_inferredGTEx.png', sep = ''), width = 4, height = 4, dpi = 300)
save(p, file = paste(outdir, 'ROC_inferredGTEx.rdata', sep = ''))
p = plot_PR(df[ grep('inferredGTEx_', df$regulon_dataset) ,],
            balanced_accuracy = balanced_accuracy,
            regulon_evidence_filter = 'inferredGTEx',
            TFs_filter = inferredGTEx_common_TFs, 
            line_colors =  my_color_palette$EMBL[4])
p
ggsave(filename = paste(outdir, 'PR_inferredGTEx.png', sep = ''), width = 7, height = 4, dpi = 300)
save(p, file = paste(outdir, 'PR_inferredGTEx.rdata', sep = ''))


idx = unique( c(grep('inferredGTEx', df$regulon_evidence),
                grep('inferredGTEx', df$regulon_dataset)))
P = auc2coverage_plot(df[ idx , ], balanced_accuracy = balanced_accuracy, y_scale = y_scale, x_scale = x_scale)
P$aucs
P$P + ggtitle(plot_title)
ggsave(filename = paste(outdir, 'combined_inferredGTEx.png', sep = ''), width = 4, height = 4, dpi = 300)








# plot chipX
p = plot_rocs(df,
              balanced_accuracy = balanced_accuracy,
              regulon_evidence_filter = 'chipSeq',
              line_colors = my_color_palette$shiny)
p
ggsave(filename = paste(outdir, 'ROC_chipSeq.png', sep = ''), width = 4, height = 4, dpi = 300)
save(p, file = paste(outdir, 'ROC_chipSeq.rdata', sep = ''))
p = plot_PR(df,
            balanced_accuracy = balanced_accuracy,
            regulon_evidence_filter = 'chipSeq',
            line_colors =  'coral')
p
ggsave(filename = paste(outdir, 'PR_chipSeq.png', sep = ''), width = 7, height = 4, dpi = 300)
save(p, file = paste(outdir, 'PR_chipSeq.rdata', sep = ''))



idx = unique( c(grep('chipSeq', df$regulon_evidence),
                grep('chipSeq', df$regulon_dataset)))
P = auc2coverage_plot(df[ idx , ], balanced_accuracy = balanced_accuracy, y_scale = y_scale, x_scale = x_scale)
P$aucs
P$P + ggtitle(plot_title)
ggsave(filename = paste(outdir, 'combined_chipSeq.png', sep = ''), width = 4, height = 4, dpi = 300)








# plot curated_databases
plot_rocs(subset(df, ! regulon_dataset %in% c('TFe', 'tfact', 'trrust', 'oreganno', 'trrd_via_tfact')  ),
          balanced_accuracy = balanced_accuracy,
          regulon_evidence_filter = 'curated_databases',
          line_colors = c(brewer.pal(11, 'Paired'), my_color_palette$EMBL))
ggsave(filename = paste(outdir, 'ROC_curated_databases.png', sep = ''), width = 4, height = 4, dpi = 300)
p = plot_PR(df,
            balanced_accuracy = balanced_accuracy,
            regulon_evidence_filter = 'curated_databases',
            line_colors =  my_color_palette$EMBL[3])
p
ggsave(filename = paste(outdir, 'PR_curated_databases.png', sep = ''), width = 7, height = 4, dpi = 300)
save(p, file = paste(outdir, 'PR_curated_databases.rdata', sep = ''))

idx = unique( c(grep('curated', df$regulon_evidence),
                grep('curated', df$regulon_dataset)))
P = auc2coverage_plot(df[ idx , ], balanced_accuracy = balanced_accuracy, y_scale = y_scale, x_scale = x_scale)
P$aucs
P$P + ggtitle(plot_title)
ggsave(filename = paste(outdir, 'combined_curated.png', sep = ''), width = 4, height = 4, dpi = 300)

idx = Reduce(intersect, list(which(df$regulon_dataset %in% c("tfact_signed", "trrust_signed", "TFe_signed", "trrd_via_tfact_signed", "oreganno_signed",
                                                             "tfact", "trrust", "TFe", "trrd_via_tfact", "oreganno") )))
P = auc2coverage_plot(df[ idx , ], balanced_accuracy = balanced_accuracy, y_scale = y_scale, x_scale = x_scale, filter_datasets = F)
P$aucs
P$P + ggtitle(plot_title)
ggsave(filename = paste(outdir, 'combined_curateddatabases_signed.png', sep = ''), width = 4, height = 4, dpi = 300)
