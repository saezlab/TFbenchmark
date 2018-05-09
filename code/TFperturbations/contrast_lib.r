library(limma)
library(affy)
library(hgu133plus2.db) # GPL570
library(hgu133a.db)  # GPL96
library(hgu133a2.db) # GPL571
library(illuminaHumanv4.db) # GPL10558
library(pd.hugene.1.0.st.v1) # GPL6244
library(illuminaHumanv4.db) # GPL10558
library(hgfocus.db) #GPL201
library(hgu95av2.db)
library(annotate)
library(GEOquery)
library(ggplot2)
library(ggrepel)


read_desigfile = function(design_file){
  x = t(read.delim(design_file, header = F, sep = '\t'))
  colnames(x) = x[1,]
  my_design = as.list(x[2,])
  my_design$positive_samples = unlist(strsplit(my_design$positive_samples, split = ','))
  my_design$negative_samples = unlist(strsplit(my_design$negative_samples, split = ','))
  return(my_design)
}




find.geneannot = function(platform){
  gene_annot = NULL
  if( platform == 'GPL570')
    gene_annot = 'hgu133plus2.db'
  if( platform == 'GPL96')
    gene_annot = 'hgu133a.db'
  if( platform == 'GPL571')
    gene_annot = 'hgu133a2.db'
  if( platform == 'GPL10558')
    gene_annot = 'illuminaHumanv4.db'
  if( platform == 'GPL6244')
    gene_annot = 'pd.hugene.1.0.st.v1'
  if( platform == 'GPL10558')
    gene_annot = 'illuminaHumanv4.db'
  if( platform == 'GPL8300')
    gene_annot = 'hgu95av2.db'
  if( platform == 'GPL201')
    gene_annot = 'hgfocus.db'
  return(gene_annot)
}

run_differential_expression = function(norm_exp, de, path_plots = NULL){
  # Map probes to gene ids
  ID = featureNames(norm_exp)
  Symbol = NULL
  gene_annot = find.geneannot(de$platform)
  if(! is.null(gene_annot) )
    try(Symbol <- getSYMBOL(ID, data = gene_annot), silent = T)
  if( is.null(Symbol) & class(norm_exp$GEOobj) == "GDS" ){
    annot = fData(GDS2eSet(norm_exp$GEOobj))
    annot$Symbol = as.character(annot$`Gene symbol`)
    Symbol = annot$Symbol[ match(ID, annot$ID) ]
  }
  if( is.null(Symbol) & class(norm_exp$GEOobj) == "GSE" ){
    x = GPLList(norm_exp$GEOobj)[[1]]
    annot = Table(dataTable(x))
    symbol_column = grep('symbol', grep('^gene', names(annot), value = T, ignore.case = T), value = T, ignore.case = T)
    if( length(symbol_column) == 1 ){
      annot$Symbol = as.character(annot[,symbol_column])
    }else{
      try(annot$Symbol <- sapply(strsplit(as.character(annot$gene_assignment), split = " \\/\\/ "), function(st) st[2]), silent = T)
      try(annot$Symbol <- as.character(annot$Symbol), silent = T)
    }
    if ( is.null(annot$Symbol) )
      cat('NOTE: Unknown filed name containing a "gene symbol"\n', 'Check:', paste(names(annot), collapse = '\n\t') ,'\n')
    Symbol = annot$Symbol[ match(ID, annot$ID) ]
  }
  fData(norm_exp$eset) = data.frame(ID=ID,Symbol=Symbol)
  # Design matrix
  design = model.matrix(~ 0+c(rep(0, length(de$negative_samples)), 
                              rep(1, length(de$positive_samples)) ) )
  design = cbind(design, (design-1)*-1 )
  colnames(design) = c("Positive", "Negative")
  rownames(design) = c(de$negative_samples, de$positive_samples)
  design = design[ rownames(design) %in% colnames(norm_exp$eset) , ]
  # Remove genes with NA
  idx = apply(norm_exp$eset[, rownames(design) ], 1, function(x) any(is.na(x)) )
  exp_m = norm_exp$eset[ ! idx , rownames(design) ]
  # Fit model with contrasts
  fit = NULL
  try(fit <- lmFit(object = exp_m, design = design))
  if(is.null(fit))
    return(NULL)
  cont.matrix = makeContrasts(TFstatus=Positive-Negative, levels=design)
  fit2 = contrasts.fit(fit, contrasts = cont.matrix)
  ebayes = eBayes(fit2, trend=TRUE)
  results = topTable(ebayes, adjust="fdr", number = nrow(norm_exp$eset))
  if( ! is.null(path_plots) ){
    try(pheatmap::pheatmap(exp_m[ results$ID[1:200], ], scale = 'row', annotation_col = as.data.frame(design), filename = paste(path_plots, de$id, '_heatmap.png', sep = ''), show_rownames = F), silent = T)
    results$Status = ifelse(results$adj.P.Val < 0.05, "Sig. FDR < 0.05", "Not Sig")
    p = ggplot(data=results, aes(x=logFC, y=-log10(P.Value), colour=Status)) +
      geom_point(alpha=0.5, size=1.75) + 
      scale_color_manual(values = c("gray", "steelblue") ) + 
      geom_text_repel(data = subset(results, Symbol == de$pathway), aes(label = Symbol), color = 'black', size = 6) + 
      geom_point(data= subset(results, Symbol == de$pathway), alpha=0.8, size=2.5, color = 'red') + 
      theme_light(18) + xlab("log2 fold change") + ylab("-log10 p-value") + ggtitle(paste(de$id, de$effect, de$treatment, sep = ' - '))
    ggsave(p, filename = paste(path_plots, de$id, '_volcano.png', sep = ''))
  }
  return(results)
}

