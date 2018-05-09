# ======================= 
# User-defined variables
# ======================= 
options(warn = -1, scipen = 1000) #R will not report any warnings (warn = -1), R will not use scientific notation (scipen = 1000) 
#=======================> DONE! 




# =====================
# Parsing Arguments
# =====================
args = commandArgs(trailingOnly = TRUE) # Read Arguments from command line

#Parsing arguments and storing values
for (each.arg in args) {
	#bed file names
	if (grepl('^-filename=',each.arg)) {
		arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] # split on =
		if (! is.na(arg.split[2]) ) {
			filename <- arg.split[2]
		} else {
			stop('No bed file names')
		} 
	}
}
options(stringsAsFactors = FALSE)
#=======================> DONE! 


standardize = function(x,a,b) {
	x = (b-a) * ( (x-min(x)) / (max(x)-min(x)) ) + a
	return(x)
}





#filename = "ReMap2_bach2_nrPeaks.bed.closest"
l = read.table(filename, sep = "\t", comment.char = "", quote = "")
name=basename(filename)

a = pdf(paste0(unique(l$V4), "-distance.pdf"))
if (length(l$V1) > 1) {
	a = hist(log10(l$V8 + 1), breaks = "scott", xlim = c(0,8), xlab = "log10[Distance to Closest TSS + 1]", main = l$V4[1], col = "cornflowerblue", xaxp = c(0,8,8))
}
a = dev.off()

if (median(l$V8) == 0) {
	med = 100 
} else {
	med = median(l$V8)
}
dist_score = exp( -( l$V8 / ( (med) * 10) ) )
a = pdf(paste0(unique(l$V4), "-distanceScore.pdf"))
a = plot(log10(l$V8+1), dist_score, xlab = "log10[Distance to Closest TSS + 1]", ylab = "Distance Score", xlim = c(0,8), xaxp = c(0,8,8), ylim = c(0,1), yaxp = c(0,1,10), main = paste0(l$V4[1], " - median distance: ", median(l$V8)))
a = dev.off()


peak_score = l$V5 * dist_score


all = data.frame(name = l$V7, score = peak_score)
all = aggregate(score ~ name, data = all, sum)
all$score = round(all$score, 3)
all = all[order(all$score, decreasing = T),]
no = which(all$score == 0)



if (length(all$score) == 1) {

	write(paste(all$name, l$V4[1], ".", 1000, sep = "\t"), ncolumns = 1, file = paste0(l$V4[1], ".assign.txt"))

} else {

	if(length(no) >= 1) {
		write(paste(all$name[-no], l$V4[1], ".", round(standardize(all$score[-no], 1, 1000),3), sep = "\t"), ncolumns = 1, file = paste0(l$V4[1], ".assign.txt"))
	} else {
		write(paste(all$name, l$V4[1], ".", round(standardize(all$score, 1, 1000),3), sep = "\t"), ncolumns = 1, file = paste0(l$V4[1], ".assign.txt"))	
	}
	
}

