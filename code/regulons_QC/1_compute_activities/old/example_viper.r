rm(list = ls())

load('/Volumes/GoogleDrive/My Drive/borrar/regulons.rdata')
load('/Volumes/GoogleDrive/My Drive/borrar/expression_matrix.rdata')
library(viper)
activities = viper(eset = E, regulon = regulon, minsize = 4, nes = T, method = 'none', eset.filter = F, pleiotropy = F, verbose = T)
dim(activities)
sessionInfo()

source('/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/code/regulons_QC/compute_activities/lib_viper/viper_1.12.r') # from viper_1.12.0.tar.gz
activities = viper(eset = E, regulon = regulon, minsize = 4, nes = T, method = 'none', eset.filter = F, pleiotropy = F, verbose = T)
dim(activities)