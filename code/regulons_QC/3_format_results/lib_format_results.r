
source('code/utils.r')
library(pROC)
library(plotROC)





build_mdf = function(file){
  load(file)
  mdf = melt(activity_ranks)
  names(mdf) = c('regulon', 'experiment', 'rank_nes')
  mdf$NES = melt(activity_nes)$value
  # mdf$pvalue1tailed = melt(activity_pvalue1tailed)$value
  mdf$perturbed_TF = columns_perturbation_annot$perturbed_TF[ match(mdf$experiment, columns_perturbation_annot$perturbation_id) ]
  mdf$TF = rows_regulons_annot$TF[ match(mdf$regulon, rows_regulons_annot$regulon_id) ]
  mdf$regulon_evidence = rows_regulons_annot$regulon_evidence[ match(mdf$regulon, rows_regulons_annot$regulon_id) ]
  mdf$regulon_dataset = rows_regulons_annot$regulon_dataset[ match(mdf$regulon, rows_regulons_annot$regulon_id) ]
  mdf$regulon_group = rows_regulons_annot$regulon_group[ match(mdf$regulon, rows_regulons_annot$regulon_id) ]
  mdf$is_TF_perturbed = mdf$TF == mdf$perturbed_TF
  mdf$experiment = as.character(mdf$experiment)
  # remove NA activities
  mdf = subset(mdf, ! is.na(NES) )
  return(mdf)
}

