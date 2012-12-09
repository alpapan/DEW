# TODO: test
# 
# Author: brian haas and alexie papanicolaou
###############################################################################

normalize = function(x){
	if(length(which(is.na(x)))==0){
		(x-mean(x))/apply(x, 2, sd)
	}
	else {
		(x-mean(x,na.rm=T))/apply(x,2,sd,na.rm=T)
	}
}

TMM_normalize = function(targetsFile, outfile='',minSumTagCounts = 10) {
	library(edgeR,quietly=T,warn.conflicts=F,verbose=F)	
	# format should be like so, including header:
#   files	group	description
#   fileA.txt  grpA    my grpA replicate
#   fileB.txt  grpB    my grpB replicate
#   ...
	
	targets = read.delim(file=targetsFile, stringsAsFactors = FALSE)
	myDGEList = readDGE(targets)
	
	# only use those rows that have a minimum number of tags (min to be DE)
	myDGEList = myDGEList[rowSums(myDGEList$counts) >= minSumTagCounts, ];
	
	# reset library sizes
	myDGEList$samples$lib.size = colSums(myDGEList$counts)
	
	# run TMM normalization
	myDGEList = calcNormFactors(myDGEList)
	
	myDGEList$samples$eff.lib.size = myDGEList$samples$lib.size * myDGEList$samples$norm.factors
	if (!is.null(outfile)){
		unlink(outfile)
		write.table(myDGEList$samples, file=outfile, quote=F, sep="\t", row.names=F)
	}
	return (myDGEList) 
}


edgeR_DE_analysis = function (targetsFile, datafile='tmp', dispersion='auto', FDR=0.05,kclusters = 10) {
	library(gdata,quietly=T,warn.conflicts=F,verbose=F)
	library(gplots,quietly=T,warn.conflicts=F,verbose=F)
	library(cluster,quietly=T,warn.conflicts=F,verbose=F)
	library(Biobase,quietly=T,warn.conflicts=F,verbose=F)
	library(rjson,quietly=T,warn.conflicts=F,verbose=F)
	library(ape,quietly=T,warn.conflicts=F,verbose=F)
	library(edgeR,quietly=T,warn.conflicts=F,verbose=F)
	library(fdrtool,quietly=T,warn.conflicts=F,verbose=F)
	
	DGEList = TMM_normalize(targetsFile)
	
	if (dispersion == 'auto'){
		DGEList <- estimateCommonDisp(DGEList)
		DGEList <- estimateTagwiseDisp(DGEList)
	}
	# set dispersion to zero for Poisson-equivalent results.    
	
	de.tgw = exactTest(DGEList, dispersion=dispersion)
	
	# qvalues with fdrtool; correlates very well with FDR but more stringent; lfdr essentially bins values to compute fdr and is even more stringent
	postscript(file=paste(datafile,'.qvalues.ps',sep=''))
	fdrresults<-fdrtool(de.tgw$table$PValue,statistic='pvalue',plot=T)
	dev.off()

	# $table: data frame containing columns for the log2-fold-change, 'logFC',
	# the average log2-counts-per-million, 'logCPM', the two-sided p-value 'PValue'
	## so log2 of 0.1375035 and -0.1520031 is 1.10 and 0.90 FC respectively
	#
	# genes with positive log-fold change are up-regulated in group B compared
	# with group A (and vice versa for genes with negative log-fold change).
	de.tgw$table$QValue <- fdrresults$qval
	de.tgw$table$lfdr <- fdrresults$lfdr
	
	# potentially housekeeping genes
	housekeeping<-rownames(subset(de.tgw$table, PValue > 0.95))
	unlink(paste(datafile,'.housekeeping.txt',sep=''))
	lapply(housekeeping, write, paste(datafile,'.housekeeping.txt',sep=''), append=T, ncolumns=1,sep ="\n")
	
	
	topTagCount = sum(de.tgw$table$QValue <= FDR); # edgeR R-2.13+ 
	toptgw = topTags(de.tgw, n = topTagCount);
	detags_toptgw = toptgw$table[toptgw$table$FDR <= FDR,];
	detags_names = rownames(detags_toptgw)
	
	postscript(file=paste(datafile,'.ps',sep=''))
	plotSmear(DGEList, de.tags=detags_names)
	dev.off()
	
	# table with entries of interest
	unlink(paste(datafile,'.results.txt',sep=''))
	write.table(de.tgw$table, quote=F, file=paste(datafile,'.all.txt',sep=''), sep="\t",col.names=NA)
	write.table(detags_toptgw, quote=F, file=paste(datafile,'.results.txt',sep=''), sep="\t",col.names=NA)

	ordered.data <- as.data.frame(matrix(t(DGEList$counts[rownames(DGEList$counts) %in% detags_names]),ncol=2))
	rownames(ordered.data) <- detags_names
	colnames(ordered.data)<-as.vector(DGEList$samples$group)
	gene_names <- detags_names
	sample_names <- colnames(ordered.data)
	
	# gene plots for 2 treatments
	pdf(file=paste(datafile,'.per_gene_plots.pdf',sep=''))
	par(mfrow=c(3, 2))
	#remember [row,column]
	for (i in 1:length(ordered.data[,1])) {
		d <- ordered.data[i,]
		ymin <- min(d);
		ymax <- max(d);
		plot(as.numeric(d), type="l", ylim=c(ymin,ymax), main=paste(gene_names[i],paste('FDR:',formatC(ordered.data[i,length(ordered.data[i,])],format='e',digits=2),sep=' ') ,sep="\n"), col="blue", xaxt="n", xlab="", ylab="TMM-norm. eff. counts")
		axis(side=1, at=1:length(d), labels=sample_names, las=0)
	}  
	dev.off()
	rm(d,i,ymin,ymax)
	
	#transform data
	logged.data <- log2(1 + as.matrix(ordered.data)) 
	if (length(logged.data[,1]) == 1){
		return
	}else{
		centered.data <- normalize(logged.data) # normalize
#		centered.data <- logged.data # dont normalize
	}
	
	# cluster data
	hc_genes <- agnes(centered.data, diss=FALSE, metric="euclidean") # cluster genes
	hc_samples <- hclust(as.dist(1-cor(centered.data, method="spearman")), method="complete") # cluster conditions
	myheatcol <- redgreen(75)
	gene_partition_assignments <- cutree(as.hclust(hc_genes), k=kclusters)
	partition_colors <- rainbow(length(unique(gene_partition_assignments)), start=0.4, end=0.95)
	gene_colors <- partition_colors[gene_partition_assignments]
	
	# heatmap
	postscript(file=paste(datafile,'.heatmap.ps',sep=''), horizontal=FALSE, width=8, height=18, paper="special")
	figure.heatmap1<-heatmap.2(centered.data, dendrogram="both", Rowv=as.dendrogram(hc_genes), Colv=as.dendrogram(hc_samples), col=myheatcol, RowSideColors=gene_colors, scale="none", density.info="none", trace="none", key=TRUE, keysize=1.2, cexCol=2.5, margins=c(15,15), lhei=c(0.4,2), lwid=c(2.5,4))
	dev.off();
	
	## data for JS interactive map
	fdr.data<-as.data.frame(detags_toptgw$FDR)
	rownames(fdr.data)<-detags_names
	fdr.data.json<-toJSON(as.data.frame(t(fdr.data)))
	figure.heatmap1.data.json<-toJSON(as.data.frame(figure.heatmap1$carpet))
	figure.heatmap1.samples.json<-toJSON(rownames(figure.heatmap1$carpet))
	figure.heatmap1.genes.json<-toJSON(colnames(figure.heatmap1$carpet))
	cat(fdr.data.json,file=paste(datafile,'.fdr.json',sep=''))
	cat(figure.heatmap1.data.json,file=paste(datafile,'.data.json',sep=''))
	cat(figure.heatmap1.samples.json,file=paste(datafile,'.samples.json',sep=''))
	cat(figure.heatmap1.genes.json,file=paste(datafile,'.genes.json',sep=''))
	hc_genes$labels <-hc_genes$order.lab # add labels
	write.tree( as.phylo(as.hclust(hc_genes)),paste(datafile,'.hcgenes.newick',sep=''))
	write.tree( as.phylo(as.hclust(hc_samples)),paste(datafile,'.hcsamples.newick',sep=''))
}


edgeR_heatmaps_all = function(datafile='tmp',kclusters = 10){
	library(gdata,quietly=T,warn.conflicts=F,verbose=F)
	library(gplots,quietly=T,warn.conflicts=F,verbose=F)
	library(cluster,quietly=T,warn.conflicts=F,verbose=F)
	library(Biobase,quietly=T,warn.conflicts=F,verbose=F)
	library(rjson,quietly=T,warn.conflicts=F,verbose=F)
	library(ape,quietly=T,warn.conflicts=F,verbose=F)
	
	# variables
	data <- read.table(file=datafile, header=T, com="", sep="\t",row.names=1)
	
	ordered.data <- data[with(data, order(data[,length(data[1,])])), ]
	## last column is now FDR column
	gene_names <- rownames(ordered.data)
	sample_names <- colnames(ordered.data[,-length(ordered.data[1,])])
	
	# gene plots
	pdf(file=paste(datafile,'.per_gene_plots.pdf',sep=''))
	par(mfrow=c(3, 2))
	#remember [row,column]
	for (i in 1:length(ordered.data[,1])) {
		d <- ordered.data[i,-length(ordered.data[1,])]
		ymin <- min(d);
		ymax <- max(d);
		plot(as.numeric(d), type="l", ylim=c(ymin,ymax), main=paste(gene_names[i],paste('FDR:',formatC(ordered.data[i,length(ordered.data[i,])],format='e',digits=2),sep=' ') ,sep="\n"), col="blue", xaxt="n", xlab="", ylab="TMM-norm. eff. counts")
		axis(side=1, at=1:length(d), labels=sample_names, las=2)
	}  
	dev.off()
	rm(d,i,ymin,ymax)
	
	# transform data
	logged.data <- log2(1 + as.matrix(ordered.data[,-length(ordered.data[1,])])) # convert to matrix remove FDR; convert to log2
	# centered.data <- normalize(logged.data) # center rows, mean substracted
	centered.data <- logged.data
	
	# cluster data
	hc_genes <- agnes(centered.data, diss=FALSE, metric="euclidean") # cluster genes
	hc_samples <- hclust(as.dist(1-cor(centered.data, method="spearman")), method="complete") # cluster conditions
	myheatcol <- redgreen(75)
	gene_partition_assignments <- cutree(as.hclust(hc_genes), k=kclusters)
	partition_colors <- rainbow(length(unique(gene_partition_assignments)), start=0.4, end=0.95)
	gene_colors <- partition_colors[gene_partition_assignments]
	
	# heatmap
	postscript(file=paste(datafile,'.heatmap.ps',sep=''), horizontal=FALSE, width=8, height=18, paper="special")
	figure.heatmap1<-heatmap.2(centered.data, dendrogram="both", Rowv=as.dendrogram(hc_genes), Colv=as.dendrogram(hc_samples), col=myheatcol, RowSideColors=gene_colors, scale="none", density.info="none", trace="none", key=TRUE, keysize=1.2, cexCol=2.5, margins=c(15,15), lhei=c(0.4,2), lwid=c(2.5,4))
	dev.off();
	
	## data for JS interactive map
	fdr.data<-as.data.frame(ordered.data$FDR)
	rownames(fdr.data)<-rownames(ordered.data)
	fdr.data.json<-toJSON(as.data.frame(t(fdr.data)))
	figure.heatmap1.data.json<-toJSON(as.data.frame(figure.heatmap1$carpet))
	figure.heatmap1.samples.json<-toJSON(rownames(figure.heatmap1$carpet))
	figure.heatmap1.genes.json<-toJSON(colnames(figure.heatmap1$carpet))
	cat(fdr.data.json,file=paste(datafile,'.fdr.json',sep=''))
	cat(figure.heatmap1.data.json,file=paste(datafile,'.data.json',sep=''))
	cat(figure.heatmap1.samples.json,file=paste(datafile,'.samples.json',sep=''))
	cat(figure.heatmap1.genes.json,file=paste(datafile,'.genes.json',sep=''))
	hc_genes$labels <-hc_genes$order.lab # add labels
	write.tree( as.phylo(as.hclust(hc_genes)),paste(datafile,'.hcgenes.newick',sep=''))
	write.tree( as.phylo(as.hclust(hc_samples)),paste(datafile,'.hcsamples.newick',sep=''))
}


