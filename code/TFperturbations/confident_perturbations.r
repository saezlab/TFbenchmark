rm(list = ls())
home = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/'
setwd(home)
source('code/utils.r')


signature_files = list.files(path = 'data/regulons_QC/pertubations/GEO/contrasts', pattern = 'rdata', full.names = T)
signature_tfs = strsplit(signature_files, '\\.') %>% sapply(., head, 1) %>% strsplit(., '/') %>% sapply(., tail, 1)
signature_ids = strsplit(signature_files, '/') %>% sapply(., tail, 1) %>% gsub('.rdata', '', .)



is_TF_DE_expected = mapply(function(tf, fi){
            DEsignature = get(load(fi))
            TF_DEsignature = subset(DEsignature, Symbol == tf)
            any(TF_DEsignature$P.Value < 0.05 & TF_DEsignature$logFC > 0)
          }, signature_tfs, signature_files)
table(is_TF_DE_expected)
which( ! is_TF_DE_expected)
is_TF_DE_opposite = mapply(function(tf, fi){
  DEsignature = get(load(fi))
  TF_DEsignature = subset(DEsignature, Symbol == tf)
  any(TF_DEsignature$P.Value < 0.05 & TF_DEsignature$logFC < 0) & ! any(TF_DEsignature$P.Value < 0.05 & TF_DEsignature$logFC > 0)
}, signature_tfs, signature_files)
table(is_TF_DE_opposite)
which( is_TF_DE_opposite )
is_TF_in_array = mapply(function(tf, fi){
  DEsignature = get(load(fi))
  tf %in% DEsignature$Symbol
}, signature_tfs, signature_files)
table(is_TF_in_array)
which( ! is_TF_in_array)



source('code/TFperturbations/contrast_lib.r')
design_files = list.files('data/regulons_QC/pertubations/GEO/design', recursive = T, pattern = 'G', full.names = T)
designs = lapply(design_files, read_desigfile)
design_id = sapply(designs, function(x) x$id)
design_effect = sapply(designs, function(x) x$effect)
design_treatment = sapply(designs, function(x) x$treatment)
design_GTEx_tissue = sapply(designs, function(x) x$GTEx_tissue)



df = data.frame(id = signature_ids, tf = signature_tfs, 
                is_TF_DE_expected = is_TF_DE_expected, 
                is_TF_DE_opposite = is_TF_DE_opposite,
                is_TF_in_array = is_TF_in_array, 
                perturbation_effect = design_effect[ match(signature_ids, design_id) ],
                perturbation_treatment = design_treatment[ match(signature_ids, design_id) ],
                perturbation_GTExtissue = design_GTEx_tissue[ match(signature_ids, design_id) ],
                stringsAsFactors = F)
write.csv(x = df, file = 'data/regulons_QC/pertubations/GEO/contrasts/TF_perturbation_confidence.csv', row.names = F)
