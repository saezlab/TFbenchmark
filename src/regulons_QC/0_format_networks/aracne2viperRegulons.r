#' Generate and format TF regulons for the TF resources benchmark 
#' available at: https://www.biorxiv.org/content/biorxiv/early/2018/06/18/337915.full.pdf 
#'
#' This code loads the individual TF-target networks from running ARACNE (network.txt file)
#' and transforms them into R objects with a format compatible with VIPER regulons (www.bioconductor.org/packages/release/bioc/html/viper.html)
#' Each row indicates an individual TF-target interaction
#' Columns indicate Regulator,	Target,	MI and pvalue respectively
#' The supplied expression file has to contain genes in HGNC symbols
#' 
#' 
#' The processing includes the steps to convert aracne networks to viper regulons as indicated in VIPER documentation
#'    (package vignette section "5 Generating the regulon object). Please, cite (Alvarez et al. Nature Genetics 2016)



### Set enviroment
rm(list = ls())
home = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/'
setwd(home)
source('code/lib/utils.r')
source('code/lib/functions.r')
source('code/lib/data.r')



### Variables
regulons_folder = 'data/TF_target_sources/inferred/GTEx/tissue_specific/'
expression_folder = 'data/TF_target_sources/inferred/GTEx/expression_matrixes/'



# Generate ARACNe tissue specific regulons - GTEx
tissues = list.dirs(regulons_folder, full.names = F)[-1]
for (tissue in tissues){
  aracne_file = paste(regulons_folder, tissue, '/network.txt', sep = '')
  expression_file = paste(expression_folder, tissue, '.txt', sep = '')
  tissue_regulon = aracne2viper_tissueSpecificNet(aracne_file, expression_file)
  if( ! is.null(tissue_regulon) ){
    names(tissue_regulon) = paste(names(tissue_regulon), gsub(' ', '_', tissue) )
    save(tissue_regulon, file = paste(regulons_folder, tissue, '/', tissue, '_viperRegulon.rdata' , sep = '') )
  }
}

