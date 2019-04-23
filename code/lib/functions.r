  
  
# --------------  
# -------------- Generate regulon object
# --------------
# Generate regulon objects from table file
#   n = file containing network information. Columns: TF (mandatory), target (mandatory), effect (optional), weigth (optional)
#   TFcensus = vector containing the TF names
#   repressors = vector containing the TF with repressors role. 
#              These are used to assing effect to their targets 
#              if it is not specified in the network "n" file
#   returns a viper regulon object
sif2viper = function(n, regulons_folder, TFcensus, repressors){ 
  message(n)
  viper_regulon = list()
  # Load network and format regulons
  net = read.delim(paste(regulons_folder, '/', n, sep = ''), stringsAsFactors = F, sep = '\t')
  # Remove duplicated interactions
  net = unique(net)
  # Remove self regulation
  net = subset(net, TF != target)
  # Remove empty targets
  net = subset(net, target != "")
  # Check number of TFs
  tfs = unique(net$TF)
  message('TFs: ', length(tfs))
  # Clean TF names with additional information, such as dorothea score (example MYC_A; where "A is the score")
  if( length(grep('_', net$TF)) > 0 ){ # This step is only needed to format Dorothea regulons where the TF has a score attached by a "_"
    net$regulon_id = net$TF
    net$TF = sapply(strsplit(net$TF, split = '_'), head, 1) # This separates the TF name from the score
  }
  # Remove overrepresented targets
  if( remove_overrepresented_targets & length(tfs) > 20 ){
    overrepresented_targets = which( table(net$target) / length(tfs) > overrepresented_targets_th) %>% names(.)
    net = subset(net, ! target %in% overrepresented_targets)
  }
  # Exclude TFs not in the TFcensus
  net = subset(net, TF %in% TFcensus)
  # Exclude TFs with less than Nmin targets
  tfs_of_interest = names(which(table(net$TF) >= Nmin))
  message('TFs with (n>=', Nmin,'): ', length(tfs_of_interest))
  if( length(tfs_of_interest) < 2 )
    next()
  #  Iterate each TF to build the reguon object and ...
  for (tf in tfs_of_interest ){
    tf_targets = subset(net, TF == tf)$target
    # Define TF mode
    if( ! is.null(net$effect) ){ # if known, use predefined
      tf_targets_sign = subset(net, TF == tf)$effect
      if( tf %in% repressors ) # unless effect equals 0, check if tf is repressor
        tf_targets_sign[ tf_targets_sign == 0 ] = - 1
      if( ! tf %in% repressors ) # unless effect equals 0, check if tf is repressor
        tf_targets_sign[ tf_targets_sign == 0  ] = 1
    }else{ # if unknown, check if tf is repressor
      tf_targets_sign = rep(1, length(tf_targets))
      if( tf %in% repressors )
        tf_targets_sign = -tf_targets_sign
    }
    if( length(tf_targets_sign) == 0 )
      next()
    # Define regulation weight
    if( ! is.null(net$weight) ){
      tf_targets_weight = subset(net, TF == tf)$weight
    }else{
      tf_targets_weight = rep(1, length(tf_targets))
    }
    si = length(tf_targets)
    # Build VIPER regulon
    tf_name = paste(tf, n, sep = ' - ')
    if( ! is.null(net$regulon_id) )
      tf_name = paste(subset(net, TF == tf)$regulon_id[1], n, sep = ' - ')
    viper_regulon[[tf_name]]$tfmode = tf_targets_sign # NOTE: ttmode has to be the first element of the list. Otherwise msviper will fail.
    names(viper_regulon[[tf_name]]$tfmode) = tf_targets
    viper_regulon[[tf_name]]$likelihood = tf_targets_weight
  }
  return(viper_regulon)
}


# Generate regulon objects from aracne and expression file
#   aracne_file = file containing aracne information
#   expression_file = file containing the target gene extression (used to infer MoR by aracne2regulon)
#   returns a viper regulon object
aracne2viper_tissueSpecificNet = function(aracne_file, expression_file){
  # Load ARACNE network file
  A = read.delim(aracne_file, stringsAsFactors = F)
  if( nrow(A) == 0 )
    return(NULL)
  # Write temporary file
  write.table(A, file = '~/tmp/network2viper.txt', sep='\t', quote = F, col.names = F, row.names = F)
  # Load expression matrix
  E = read.delim(file = expression_file, stringsAsFactors = F)
  # Format gene names
  rownames(E) = E$gene
  E = as.matrix(E[,-1])
  # Run aracne2regulon function from VIPER to infer tfmode and likelihood
  require(viper)
  viper_regulon = aracne2regulon(afile = '~/tmp/network2viper.txt', eset = E, verbose = F, format = '3col')
  return(viper_regulon)
}



# --------------  
# -------------- Load and aggregate regulon objects into VIPER Regulon object
# --------------
load_and_merge_regulons = function(regulon_files, regulons_folder, filter_TFs = NULL){
  aggregated_networks = list()
  for (n in network_files){
    message(' - ', n)
    # Load network and fromat regulons
    net = get(load(paste(regulons_folder, '/', n, sep = '')))
    # Clean TF regulon names (such as the Dorothea ones containing the score in their label, example MYC_A)
    net_tfs = sapply(strsplit(names(net), split = ' '), head, 1)
    net_tfs = sapply(strsplit(net_tfs, split = '_'), head, 1)
    # Filter TFs to includo only those perturbed in the experiment
    if (! is.null(filter_TFs))
      net = net[ net_tfs %in% filter_TFs ]
    if( length(net) < 2 ) # At least 2 TF regulons required
      next()
    aggregated_networks = append(aggregated_networks, net)
  }
  return(aggregated_networks)
}







# --------------  
# -------------- Load and aggregate regulon objects into a DATA FRAME
# --------------
load_and_merge_regulons_asdf = function(networks2merge){
  regulon = list()
  for (n in networks2merge){
    message(n)
    # Load network and fromat regulons
    net = read.delim(n, stringsAsFactors = F, sep = '\t')
    evidence = unlist(strsplit(n, '/'))[3]
    database = unlist(strsplit(n, '/'))[4]
    if( is.null(regulon[[evidence]]) )
      regulon[[evidence]] = list()
    if( is.null(net$effect) )
      net$effect = 0
    for( tf in unique(net$TF) ){
      regulon[[evidence]][[database]][[tf]] = unique(rbind(regulon[[evidence]][[database]][[tf]],
                                                           subset(net, TF == tf)[, c('target', 'effect') ] ))
    }
  }
  # build merged data.frame ~ database
  df = unique(melt(regulon, id.vars = names(regulon[[1]][[1]][[1]]) ))
  df = df[, c(3,1,2,5,4)]
  names(df) = c('TF', 'target', 'effect', 'evidence', 'database')
  return(df)
}






# --------------  
# -------------- Managing TF-target effect in a consensus network df
# --------------
fix_sign_conflict = function (ddf){
  # Identify TF-target interactions whith effect sign conflict (i.e. 1 and -1)
  x = unique(subset(ddf, effect != 0 & evidence == 'curated_databases')[, c('TF', 'target', 'effect', 'TFtarget')])
  TFtargte_with_effect_confict = unique(x$TFtarget[ duplicated(x$TFtarget) ])
  message('NOTE: ', length(TFtargte_with_effect_confict), ' whith effect sign conflict.\nPrioritizing those with matching inferredGTEx effect sign')
  print(TFtargte_with_effect_confict)
  # If in inferredGTEx, then use this effect
  with_coexpression = subset(ddf, evidence == 'inferredGTEx' & TFtarget %in% intersect(TFtargte_with_effect_confict, ddf$TFtarget[ ddf$evidence == 'inferredGTEx' ]) )
  # If NOT in inferredGTEx, then prioritize according database
  without_coexpression = subset(ddf, evidence == 'curated_databases' & TFtarget %in% setdiff(TFtargte_with_effect_confict, ddf$TFtarget[ ddf$evidence == 'inferredGTEx' ]) )
  without_coexpression = rbind(subset(without_coexpression, database == 'TFe_signed'),
                               subset(without_coexpression, database == 'tfact_signed'),
                               subset(without_coexpression, database == 'trrust_signed'),
                               subset(without_coexpression, database == 'trrd_via_tfact_signed'),
                               subset(without_coexpression, database == 'NFIRegulomeDB'),
                               subset(without_coexpression, database == 'oreganno_signed')  )
  # generate hash-maps
  with_coexpression = unique(with_coexpression[, c('TF', 'target', 'effect', 'TFtarget')])
  without_coexpression = unique(without_coexpression[ ! duplicated(without_coexpression$TFtarget) , c('TF', 'target', 'effect', 'TFtarget')])
  # fix effect sign
  idx = ddf$TFtarget %in% with_coexpression$TFtarget & ddf$evidence == 'curated_databases' & ddf$effect != 0
  ddf$effect[ idx ] = with_coexpression$effect[ match(ddf$TFtarget[ idx ], with_coexpression$TFtarget) ]
  idx = ddf$TFtarget %in% without_coexpression$TFtarget & ddf$evidence == 'curated_databases' & ddf$effect != 0
  ddf$effect[ idx ] = without_coexpression$effect[ match(ddf$TFtarget[ idx ], without_coexpression$TFtarget) ]
  return(ddf)
}


assign_TFTG_sign = function(ddf){
  ddf$TFtarget = paste(ddf$TF, ddf$target)
  # order: prioritize evidence curated (curated > inferred) and signed effect (1/-1 > 0)
  ddf = ddf[ order(ddf$database, decreasing = T), ]
  ddf = ddf[ order(ddf$evidence), ]
  ddf = ddf[ order(abs(ddf$effect), decreasing = T), ]
  ddf = ddf[ order(ddf$TFtarget), ]
  new_ddf = ddf[ ! duplicated(ddf$TFtarget) , 1:3]
  return(new_ddf)
}






# --------------  
# -------------- Load matching TCGA regulon, if exists
# --------------
find_matching_TCGA_regulon = function(design_file, TCGA_folder, filter_TFs = NULL){
  source('code/lib/contrast.r')
  message('- Reading perturbation DE design file to extract sample\'s cancer type')
  tissue = read_desigfile(design_file)$TCGA_label
  tissue
  if (tissue == 'none' ){
    message('\tNOTE: non-cancer samples used in this perturbation\n')
    return(NULL)
  }
  
  regulon_file = grep(paste(tissue, '_viperRegulon', sep = ''), network_files, value = T)
  if ( length(regulon_file) == 0){
    message('\tNOTE: no matching regulon for the cancer-type: ', toupper(tissue), '\n')
    return(NULL)
  }
  
  message('- Regulon found for the cancer type: ', toupper(tissue) )
  message('- Loading cancer-specific regulon and filter TFs in perturbations')
  networks = get(load(paste(TCGA_folder, regulon_file, sep = '')))
  if( ! is.null(filter_TFs) ){
    message('- Filter TFs in perturbations')
    networks = networks[ sapply(strsplit(names(networks), ' '), head, 1) %in% filter_TFs ]
  }
}





# --------------  
# -------------- Load matching GTEx regulon, if exists
# --------------
find_matching_GTEx_regulon = function(design_file, GTEx_folder, filter_TFs = NULL){
  source('code/lib/contrast.r')
  message('- Reading perturbation DE design file to extract sample\'s tissue type')
  tissue = read_desigfile(design_file)$GTEx_tissue
  
  regulon_file = grep(paste(tissue, '_viperRegulon', sep = ''), network_files, value = T)
  if ( length(regulon_file) == 0){
    message('\tNOTE: no matching regulon for the GTEx tissue: ', toupper(tissue), '\n')
    return(NULL)
  }
  
  message('- Regulon found for the tissue: ', toupper(tissue) )
  message('- Loading tissue-specific regulon')
  networks = get(load(paste(GTEx_folder, regulon_file, sep = '')))
  if( ! is.null(filter_TFs) ){
    message('- Filter TFs in perturbations')
    networks = networks[ sapply(strsplit(names(networks), ' '), head, 1) %in% filter_TFs ]
  }
}






# --------------  
# -------------- Compute TF activities from DE signatures
# --------------
DEs2activities = function(DE_file, design_folder, networks){
  source('code/lib/contrast.r')
  
  message('- Loading DE file: ', DE_file)
  DEGs = get(load(DE_file))
  DEGs = subset(DEGs, Symbol != "" )
  DEGs = subset(DEGs, ! duplicated(Symbol))
  
  message('- Reading design file to extract tissue')
  perturbation_id = tail(unlist(strsplit(DE_file, split = '/')),1) %>% gsub('.rdata', '', .) %>% strsplit(., split = '\\.') %>% unlist(.)
  design_file = list.files(design_folder, recursive = T, full.names = T) %>% grep(perturbation_id[1], ., value = T) %>% grep(perturbation_id[2], ., value = T)
  tissue = read_desigfile(design_file)$GTEx_tissue
  
  message('- Formating data for analysis')
  myStatistics = matrix(DEGs$logFC, dimnames = list(DEGs$Symbol, 'logFC') )
  myPvalue = matrix(DEGs$P.Value, dimnames = list(DEGs$Symbol, 'P.Value') )
  # Although Gene Expression Signature (GES) can be defined by the t-statistic, 
  # to be consistent with the z-score based null model for msVIPER (see section 6.2 in its vignette), 
  # we will estimate z-score values for the GES as indicated in the vignette.
  mySignature = (qnorm(myPvalue/2, lower.tail = FALSE) * sign(myStatistics))[, 1]
  mySignature = mySignature[order(mySignature, decreasing = T)]
  
  message('- Running VIPER ...')
  mrs = msviper(mySignature, networks, verbose = T, minsize = Nmin, ges.filter = F)
  activities = data.frame(Regulon = names(mrs$es$nes),
                       Size = mrs$es$size[ names(mrs$es$nes) ], 
                       Size_adaptative = sapply(mrs$regulon[ names(mrs$es$nes) ], function(x) sum((x$likelihood/max(x$likelihood))^2) ), # relevant only in weighted regulons  
                       NES = mrs$es$nes, 
                       p.value = mrs$es$p.value, 
                       FDR = p.adjust(mrs$es$p.value, method = 'fdr'),
                       perturbation_tissue = tissue)
  
  activities = activities[ order(activities$p.value, decreasing = F) , ]
  message('')
  return(activities)
}



# --------------  
# -------------- Load and filter differential activities from msVIPER
# --------------
load_activity_msVIPERfile = function(f, Nmin, Nmax){
  # Load
  viper_results = get(load(f))
  # Filter according to the regulon size
  viper_results = subset(viper_results, Size >= Nmin & Size <= Nmax)
  viper_results = viper_results[ order(viper_results$Regulon), ]
  return(viper_results)
}


# --------------  
# -------------- Load and filter activity matrixes from VIPER
# --------------
load_activity_VIPERfile = function(f, Nmin, Nmax){
  nes = get(load(f))
  if( length(grep('specific', f)) == 0 )
    nes = nes[ nes[,'Size'] <= N_max & nes[,'Size'] >= N_min , ]
  if( length(grep('noCOMBAT', f)) == 1 ) #noCOMBAT analysis requested by a referee. Ignore it otherwise
    rownames(nes) = paste(rownames(nes), '_noCOMBAT', sep = '') # This adds the "_noCOMBAT" label to the regulon name.
  if( ! is.matrix(nes) )
    return(NULL)
  if( nrow(nes) == 0)
    return(NULL)
   viper_results = melt(nes[,-1])
  names(viper_results) = c('Regulon', 'Sample', 'NES')
  return(viper_results)
}


# --------------  
# -------------- Generate the aggregated NES matrix including all regulons
# --------------
aggregate_activities_NESmatrix = function(activities_list){
  df = melt(activities_list, id.vars = names(activities_list[[1]]))
  df$value = df$NES
  activity_nes = acast(df, Regulon ~ Sample, fill = NA)
  return(activity_nes)
}



# --------------  
# -------------- Rank activity_nes 
# --------------
activity_nes2ranks = function(activity_nes){
  activity_ranks = t(apply(activity_nes, 1, rank, ties.method = 'min'))
  activity_ranks = activity_ranks / apply(activity_ranks, 1, max)
  activity_ranks[ is.na(activity_nes) ] = NA
  return(activity_ranks)
}







# --------------  
# -------------- Merge activity_ranks, columns_perturbation_annot and rows_regulons_annot to plot accuracy
# --------------
merge_data2plot_accuracy = function(file){
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





# --------------  
# -------------- Extract regulons information from regulons name
# --------------
annotate_regulons = function(regulons){
  rows_regulons_annot = data.frame(regulon_id = regulons,
                                   TF = manage_datasets_names(regulons, what = 'regulon_id2TF'),
                                   regulon_dataset_full = manage_datasets_names(regulons, what = 'regulon_id2regulon_dataset_full'),
                                   stringsAsFactors = F)
  rows_regulons_annot$regulon_evidence = manage_datasets_names(rows_regulons_annot$regulon_dataset_full, what = 'regulon_dataset_full2evidence')
  rows_regulons_annot$regulon_dataset = manage_datasets_names(rows_regulons_annot$regulon_dataset_full, what = 'regulon_dataset_full2dataset')
  tissue_specific_index = grep('/', rows_regulons_annot$regulon_dataset_full, invert = T)
  rows_regulons_annot$regulon_dataset[ tissue_specific_index ] =  paste('inferred', rows_regulons_annot$regulon_dataset_full[tissue_specific_index] ) 
  rows_regulons_annot$regulon_evidence[ tissue_specific_index ] = 'inferred'
  rows_regulons_annot$regulon_group = manage_datasets_names(rows_regulons_annot$regulon_dataset, what = 'dataset2group')
  rows_regulons_annot$regulon_evidence[ rows_regulons_annot$regulon_group == 'cancer-specific' ] = 'inferred_cancer'
  unique(rows_regulons_annot[, c('regulon_dataset', 'regulon_evidence', 'regulon_group') ])
  return(rows_regulons_annot)
}


# --------------  
# -------------- Manage regulons information
# --------------
manage_datasets_names = function(x, what = 'dataset2group'){
  # message('Options for what are:')
  # message('- regulon_id2TF')
  # message('- regulon_id2regulon_dataset_full')
  # message('- regulon_dataset_full2evidence')
  # message('- regulon_dataset_full2dataset')
  # message('- dataset2group')
  x_out = x
  if(what == 'regulon_id2TF'){
    x_out = sapply(strsplit(x, ' '), head, 1)
    if( length(grep('_[A-E]$', x_out)) > 0 )
      x_out[ grep('_[A-E]$', x_out) ] = sapply(strsplit(x_out[ grep('_[A-E]$', x_out) ], '_'), head, 1)
  } 
  if(what == 'regulon_id2regulon_dataset_full'){
    x_out = sapply(strsplit(x, ' '), tail, 1) %>% gsub('.sif', '', .)
  }
  if(what == 'regulon_dataset_full2evidence'){
    x_out = sapply(strsplit(x, '/'), head, 1)
    x_out[ intersect(grep('PANCANCER', x, ignore.case = T), grep('omnipath', x)) ] = 'omnipath_scores_cancer'
    x_out[ intersect(grep('pancancer', x, ignore.case = T), grep('inferred', x)) ] = 'inferred_cancer'
  }
  if(what == 'regulon_dataset_full2dataset'){
    x_out = sapply(strsplit(x, '/'), function(xx) paste(xx[2:length(xx)], collapse = '_') ) %>%
      gsub('network', '', .) %>%
      gsub('__', '_', .) %>% 
      gsub('_signed_signed', '_signed', .)  
    x_out = sapply(strsplit(x_out, 'via'), head, 1)
    x_out = gsub('_$', '', x_out)  %>% 
      gsub('^_', '', .) 
  }
  if(what == 'dataset2group'){
    x_out[ grep('pantissue', x) ] = 'GTEx pantissue'
    x_out[ intersect(grep('noCOMBAT', x), grep('pantissue', x)) ] = 'GTEx pantissue_noCOMBAT'
    x_out[ grep('pancancer', x) ] = 'pancancer'
    x_out[ grep('inferred ', x) ] = 'GTEx tissue-specific'
    x_out[ intersect(grep('inferred ', x), grep('noCOMBAT', x)) ] = 'GTEx tissue-specific noCOMBAT'
    x_out[ grep('tcga_', x) ] = 'cancer-specific'
    # x_out[ x_out %in% c(LETTERS[1:5], 'TOP') ] = 'omnipath'
    # x_out[ x_out %in% paste(c(LETTERS[1:5], 'TOP'), 'PANCANCER', sep ='_') ] = 'omnipath_cancer'
    x_out[ grep('hocomoco', x_out) ] = 'hocomoco'
    x_out[ grep('ReMap', x_out) ] = 'ReMap'
    x_out[ grep('jaspar', x_out) ] = 'jaspar'
    # x_out[ x_out %in% c(LETTERS[1:5], 'TOP')  ] = 'omnipath_scores'
    x_out = gsub('_signed_signed', '', x_out)
    x_out = gsub('e2_', '', x_out)
  }
  if(what == 'evidence2color'){
    # x_out = c('coral', my_color_palette$EMBL[c(5,5,8,3,4,2,6,1,7)])[
    x_out = c('coral', my_color_palette$EMBL[c(5,5,3,3,4,2,6,1)], brewer.pal(n = 6, name = 'Paired')[6])[
      match(x_out,
            c('ChIP_Seq', 'consensus', 'consensus_all', 'consensus_curated', 'curated_databases', 'inferred',  'old_consensus', 'TFBS_scanning', 'omnipath_scores', 'omnipath_scores_cancer') )]
  }
  if(what == 'evidence2shape'){
    x_out = c(15, 6, 6, 6, 16:17, 5, 18, 16, 16, 1)[
      match(x_out,
            c('ChIP_Seq', 'consensus', 'consensus_all', 'consensus_curated', 'curated_databases', 'inferred',  'old_consensus', 'TFBS_scanning', 'omnipath_scores', 'omnipath_scores_cancer') )]
  }
  # x_out = gsub('_databases', '', x_out) %>% gsub('databases', '', .) %>% gsub('_scanning', '', .) 
  return(x_out)
}


