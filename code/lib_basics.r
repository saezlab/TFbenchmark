library(reshape2)
library(ggplot2)
library(Matrix)
library(plyr)
library(scales)
library(igraph)
library(gplots)
library(devtools)
library(RColorBrewer)
library(ggfortify)
library(qvalue)
library(mixtools)
library(ggrepel)
library(cowplot)
library(pheatmap)
library(UpSetR)
library(magrittr)
#source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")


find.cutoff <- function(x, proba=0.5) {
  x = x[! is.na(x)]
  if( length(x) > 10){
    k = as.numeric(names(sort(-table(ncomp[,5])))[1] ) # get the mode (median, mean, mode)
    model = normalmixEM(x, k = k)
    i = which.min(model$mu)  # Index of component with lower mean
    ## Cutoff such that Pr[drawn from bad component] == proba
    f <- function(x) {
      proba - (model$lambda[i]*dnorm(x, model$mu[i], model$sigma[i]) /
                 (model$lambda[1]*dnorm(x, model$mu[1], model$sigma[1]) + model$lambda[2]*dnorm(x, model$mu[2], model$sigma[2])))
    }
    if( k ==3 ){
      if(k == 2){
        plot(model, which=2, breaks = 50, col2= c('brown2', 'steelblue3'))
      }else{
        plot(model, which=2, breaks = 50, col2= ggcolors(k))
      }
      lines(density(x), lty=2, lwd=2)
    }
    if( min(x) >= 0 | max(x) < 0 ){
      return(NA)   
    }else{
      value = uniroot(f=f, lower=min(x), upper=max(x))$root 
    }
    abline(v = value, col='gray40')
    return(value)  # Careful with division by zero if changing lower and upper
  }
}



wt = function(x, file, quote=F, sep="\t", row.names=F, col.names=F, ...){
  write.table(x, file, quote=quote, sep=sep, row.names=row.names, col.names=col.names, ...)
}

rt = function(file, header=F, ...){
  read.delim(file=file, header=header, stringsAsFactors=F, ...)
}

spl = function(x, split, ...){
  unlist(strsplit(x, split=split, ...))
}

past = function(...){
  paste(..., sep="")
}

asdf = function(...){
  as.data.frame(..., stringsAsFactors=F)
}
  

DF = function(...){
  data.frame(..., stringsAsFactors=F)
}

intw = function(x, y){
  intersect(which(x), which(y))
}

capitalizeFirst = function(x) {
  s = strsplit(x, " ")[[1]]
  past(toupper(substring(s, 1,1)), substring(s, 2), collapse=" ")
}


intersectSeveral = function(...) { 
  Reduce(intersect, ...)
} 


multigrep = function(patterns, x, ...){
  unlist(lapply(patterns, grep, x = x , ...))
}


list2binaryMat = function(li, li.l, thresold=1){
  cat("- melt\n")
  dat.melt = melt(unlist(li, recursive=F)[ - which((li.l == 0) ) ])
  dat.melt$value = as.character(dat.melt$value)
  tvalue = asdf(table(dat.melt$value))
  cat("- filter\n")
  dat.melt = dat.melt[ dat.melt$value %in% tvalue$Var1[ tvalue$Freq >=  thresold ] , ]
  cat("- cast\n")
  values = unique(dat.melt$value)
  cases = unique(dat.melt$L1)
  if( length(values) > 100 ){ # if more than 500 values then split matrix and aggregate as dcast is sensible to values length
    cat("\taggregating: ")
    AggSequence = mapply(function(i,j) i:j, seq(1, length(values), 100), c(seq(1, length(values), 100)[-1]-1, length(values))   )
    dat.cast = data.frame(cases=cases) # this is to make sure dcast will not remove any case (as cases without values are removed)
    for(i in AggSequence){
      cat(i[1], ", ")
      X = rbind(DF(L1=cases, value="all")  , dat.melt[ dat.melt$value %in% values[i] ,])
      d = dcast(X, L1~value, function(x) 1, fill = 0)
      rnames =  d$L1
      dat.cast = cbind(dat.cast, d[-1])
      rownames(dat.cast) = rnames
    }
    cat("\n")
  }else{
    dat.cast = dcast(dat.melt, L1~value, function(x) 1, fill = 0)
    rownames(dat.cast) =  dat.cast$L1
  }
  dat.cast = dat.cast[ , ! colnames(dat.cast) %in% c("all", "cases")]
  return(dat.cast)
}


trans_pval = function(pval, asterisk=F){
  if( asterisk ){
    if( pval < 0.001 ){
      pval = "***"
    } else if( pval < 0.01 ){
      pval = "**"
    }else if( pval < 0.05 ){
      pval = "*"
    }
  }else{
    pval = signif(pval, 3)
    if( pval < 0.001 ){ pval = past("p < 0.001") }else{ pval = past("p < ", pval) }
  }
  return(pval)
}


ggcolors = function(n=6, h=c(0, 360) +15){
  if ((diff(h)%%360) < 1) h[2] = h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot = function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots = c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout = matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx = as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

g_legend<-function(p){
  tmp = ggplotGrob(p)
  leg = which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend = tmp$grobs[[leg]]
  return(legend)
}


get_rgb = function(colorname){
  colors = as.list(col2rgb(colorname))
  names(colors) = rownames(col2rgb(colorname))
  colors = as.data.frame(colors)
  return(colors)
}


plot_igraph = function(g, l=NULL){
  if ( is.null(l) ) {
    l <- layout.fruchterman.reingold(g,niter=1000,area=30*vcount(g)^2,repulserad=vcount(g)^4)
  }
  #l = layout.fruchterman.reingold(g,niter=500,area=vcount(g)^2,repulserad=vcount(g)^2.9)
  if( ! is.null(V(g)$color) ){
    vertex.color=V(g)$color
  }else{
    vertex.color=rgb(get_rgb(ggcolors(20)[11]), alpha=151, maxColorValue=256)
  }
  plot(g, 
       layout=l,
       edge.width=1,
       edge.arrow.size=0.25,
       vertex.size=10,
       vertex.label.color="gray10",
       vertex.label.dist=0, 
       vertex.label.cex=0.75,
       vertex.label.font=1,
       vertex.label.family='sans',
       vertex.frame.color='white',
       vertex.color=vertex.color)
}


get.igraph.all.shortest.paths = function(g, nodes, directed = F, onlymin=F){
   sh = lapply(nodes, function(n) get.all.shortest.paths(g, n, setdiff(nodes, n), mode = 'all')$res )
   if ( onlymin ){
     sh = sh[ sapply(sh, length) > 1 ]
     sh = lapply(sh, function(x) x[ sapply(x, length) == min(sapply(x, length)) ] )
   } 
   ssh = lapply(sh, function(ss) lapply(ss, function(s) rep(s, times=1, each = 2)[ -c(1, length(s)*2) ] )   )
   ids = unique(get.edge.ids(g, unlist(ssh)))
   gg = graph.data.frame(d = get.edges(g, unique(get.edge.ids(g, unlist(ssh)))), directed = directed)
   return(gg)
}



makeTransparent = function(..., alpha=0.5) {
  if(alpha<0 | alpha>1) stop("alpha must be between 0 and 1")
  alpha = floor(255*alpha)  
  newColor = col2rgb(col=unlist(list(...)), alpha=FALSE)
  .makeTransparent = function(col, alpha) {
    rgb(red=col[1], green=col[2], blue=col[3], alpha=alpha, maxColorValue=255)
  }
  newColor = apply(newColor, 2, .makeTransparent, alpha=alpha)
  return(newColor)
  
}




effectsize_cohensD  =  function(x, y) { # Effect size
  lx  =  length(x)- 1
  ly  =  length(y)- 1
  # md   =  abs(mean(x) - mean(y)) ## mean difference (numerator)
  md   =  mean(x) - mean(y) ## mean difference (numerator)
  csd  =  lx * var(x) + ly * var(y)
  csd  =  csd/(lx + ly)
  csd  =  sqrt(csd)                     ## common sd computation
  cd   =  md/csd                        ## cohen's d
}
