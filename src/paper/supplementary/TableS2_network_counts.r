rm(list = ls())
home = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/'
setwd(home)
source('code/utils.r')


nets_files = list.files(path = 'data/TF_target_sources/', full.names = T, recursive = T, pattern = 'viper') %>% grep('network', ., value = T) 
nets_files = grep('omni', nets_files, invert = T, value = T) %>% grep('old', ., invert = T, value = T) %>% grep('evidence', ., invert = T, value = T) %>% grep('cnio', ., invert = T, value = T) %>% grep('namic', ., invert = T, value = T) 
nets_id = lapply(strsplit(nets_files, split = "/"), tail, 3) %>% sapply(., paste, collapse = '/') %>% gsub('_viperRegulon.rdata', '', .)



out = data.frame(regulon_id = nets_id, 
                 regulon_dataset = manage_datasets_names(nets_id, what = 'regulon_dataset_full2dataset'),
                 regulon_source = sapply(strsplit(manage_datasets_names(nets_id, what = 'regulon_dataset_full2dataset'), '_'), head, 1),
                 regulon_evidence = manage_datasets_names(nets_id, what = 'regulon_dataset_full2evidence'),
                 TFs = 0 , targets = 0,  TF_target_interactions = 0 , stringsAsFactors = F)
rownames(out) = nets_files
head(out)
all_regulons = list()
for( n in nets_files ){
  net = get(load(n))
  names(net) = sapply(strsplit(names(net), ' '), head, 1)
  for( tf in names(net) )
    all_regulons[[tf]] = unique(c(all_regulons[[tf]], names(net[[tf]]$tfmode)))
  out[n,]$TFs = length(net)
  out[n,]$targets = lapply(net, function(x) x$tfmode) %>% lapply(., names) %>% unlist(.) %>% unique(.) %>% length(.)
  out[n,]$TF_target_interactions = lapply(net, function(x) x$tfmode) %>% lapply(., names) %>% unlist(.) %>% length(.)
}
out$regulon_dataset = gsub('hocomocov11_n', 'hocomocov11_top', out$regulon_dataset ) %>% gsub('jaspar2018_n', 'jaspar2018_top', .) %>% gsub('inferredGTEx_', 'inferredGTEx_Ntissues', .) 
head(out)
tail(out)
write.csv(out[,-1], row.names = F, file = 'paper/GarciaAlonso_tableS2_unformated.csv', quote = F)


total = melt(all_regulons)
nrow(unique(total))
length(unique(total$L1))
length(unique(total$value))
