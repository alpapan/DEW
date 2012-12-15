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

TMM_normalize = function(targetsFile, minCPM = 2,minLibs = 1,outfile='') {
	library(edgeR,quietly=T,warn.conflicts=F,verbose=F)	
	# format should be like so, including header:
#   files	group	description
#   fileA.txt  grpA    my grpA replicate
#   fileB.txt  grpB    my grpB replicate
#   ...
	
	targets = read.delim(file=targetsFile, stringsAsFactors = FALSE)
	DGEList = readDGE(targets)
	# how many factors are we going to look at ?
	DGEList$factors_to_do<-c()
	for(i in colnames(DGEList$samples)) { 
		if (i %in% c("files","description","lib.size","norm.factors","eff.lib.size", "group")){next}; 
		DGEList$factors_to_do<-append(DGEList$factors_to_do,i)
	}
	
	# only use those rows that have a minimum number of tags (min to be DE)
	#TODO actually use minLibs in perl script.
	DGEList<-DGEList[rowSums(cpm(DGEList) > minCPM) >= minLibs,]
	
	# reset library sizes
	DGEList$samples$lib.size = colSums(DGEList$counts)
	
	# run TMM normalization
	DGEList = calcNormFactors(DGEList)
	rownames(DGEList$samples)<-DGEList$samples$description
	colnames(DGEList$counts)<-DGEList$samples$description
	
	DGEList$samples$eff.lib.size = DGEList$samples$lib.size * DGEList$samples$norm.factors
	# normalized CPM
	DGEList$CPM <- cpm(DGEList,normalized.lib.sizes=T)
	DGEList$logCPM <- log2(1+DGEList$CPM)
	
	# just in case normalization really changed things
	DGEList<-DGEList[rowSums(DGEList$CPM >= minCPM) >= minLibs,]
	
	if (!is.null(outfile)){
		unlink(outfile)
		write.table(DGEList$samples, file=outfile, quote=F, sep="\t", row.names=F)
	}
	return (DGEList) 
}

# TODO split between those with reps and without reps
edgeR_DE_explore = function (targetsFile, baseout='tmp', dispersion='auto', FDR=0.05,kclusters = 10,minCPM = 2,minLibs = 1){
	library(edgeR,quietly=T,warn.conflicts=F,verbose=F)
	DGEList = TMM_normalize(targetsFile,minCPM,minLibs )
	if (length(DGEList)==0){return}
	
	if (length(DGEList$factors_to_do) == 0){
		pdf(file=paste(baseout,'.qc_graphs.pdf',sep=''))
		plotMDS(DGEList, main="Biological Co-efficient of Variation distance MDS plot");
		plotMDS(DGEList$logCPM, main="Counts per million distance MDS plot")
		if (dispersion == 'auto'){
			DGEList <- estimateCommonDisp(DGEList)
			DGEList <- estimateTrendedDisp(DGEList)
			DGEList <- estimateTagwiseDisp(DGEList)
			plotBCV(DGEList, main="Dispersion trends")
			# set dispersion to zero for Poisson-equivalent results.	
		}
		dev.off()
			cat ("Processing one factor using CML\n");
			de.tgw = exactTest(DGEList, dispersion=dispersion)
			de.tgw <- edgeR_DE_postanalysis(DGEList, de.tgw, baseout, dispersion, FDR,kclusters)
	}else{
		cat ("Processing multiple factors using additive GLMs\n");
		pdf(file=paste(baseout,'.qc_graphs.pdf',sep=''))
		plotMDS(DGEList, main="Biological Co-efficient of Variation distance MDS plot");
		plotMDS(DGEList$logCPM, main="Counts per million distance MDS plot")
		group <- as.factor(DGEList$samples[, 'group'])
		# many factors, do glm
		form <- as.formula(
				paste("~ group + ", paste(
								paste('as.factor(DGEList$samples[, DGEList$factors_to_do[', 1:length(DGEList$factors_to_do),']])', sep="")
								, collapse= "+")
				)
		)
		design<-model.matrix(form)
		rownames(design) <- colnames(DGEList)
		DGEList <- estimateGLMCommonDisp(DGEList, design, verbose=TRUE)
		DGEList <- estimateGLMTrendedDisp(DGEList, design)
		DGEList <- estimateGLMTagwiseDisp(DGEList, design, prior.n=3) # inspect BCV ?!
		plotBCV(DGEList, main="Dispersion trends")
		dev.off()
		
		fit <- glmFit(DGEList, design)
		if (length(fit$fail) >0){return} 
		
		# Explore the effect of every factor
		#coeff 1 is intercept;
		# coeff 2 is group1 (starts with 'group')
		# then the other coefficients could be other groups or other factors.
		coefficients <-colnames(fit$coefficients)
		positives <-c()
		names_positives<-c()
		# we will use this later
		de.tgw <- glmLRT(DGEList, fit, coef=2)
		positives <-append(positives,sum(p.adjust(de.tgw$table$PValue, method="BH")< FDR))
		names_positives<- append(names_positives,(paste('Coefficient group vs others:',sep="")))
		
		for (i in 3:length(coefficients)){
			de.tgw.t <- glmLRT(DGEList, fit, coef=i)
			positives <-append(positives,sum(p.adjust(de.tgw.t$table$PValue, method="BH")< FDR))
			names_positives<- append(names_positives,(paste("Coeffs ",i,' vs others:',sep="")))
		}
		if (length(coefficients) > 3){
			for (i in 3:length(coefficients)){
				for (k in 4:length(coefficients)){
					if (k == i){next}
					de.tgw.t <- glmLRT(DGEList, fit, coef=i:k)
					positives <-append(positives,sum(p.adjust(de.tgw.t$table$PValue, method="BH")< FDR))
					names_positives <- append(names_positives,(paste("Coeffs ",i,':',k,' vs others:',sep="")))
				}
			}
		}
		rm(de.tgw.t)
		positives <-as.data.frame(positives)
		rownames(positives)<-names_positives
		write.table(positives,row.names=T,quote=F,sep="\t",file=paste(baseout,".factors.positives",sep=""))
		write(c("factors:"),file=paste(baseout,".factors",sep=""))
		write.table(DGEList$factors_to_do,col.names=F,row.names=T,quote=F,sep="\t",file=paste(baseout,".factors",sep=""),append=T)
		write(c("\ncoefficients:"),file=paste(baseout,".factors",sep=""),append=T)
		write.table(coefficients,col.names=F,row.names=T,quote=F,sep="\t",file=paste(baseout,".factors",sep=""),append=T)
		write(c("\nTests:"),file=paste(baseout,".factors",sep=""),append=T)
		write.table(col.names=F,quote=F,row.names=F,length(de.tgw$table$PValue),file=paste(baseout,".factors",sep=""),append=T)
		write(c("\nPositives past ",FDR,':'),file=paste(baseout,".factors",sep=""),append=T)
		write.table(positives,col.names=F,row.names=T,quote=F,sep="\t",file=paste(baseout,".factors",sep=""),append=T)
		
		# give us the results for the group factor
		de.tgw <- edgeR_DE_postanalysis(DGEList, de.tgw, baseout, dispersion, FDR,kclusters)
	}
#	return (de.tgw)
}

edgeR_DE_analysis_no_reps_housekeeping = function (DGEList, baseout='tmp', dispersion='auto', FDR=0.05,kclusters = 10) {
	#TODO use housekeeping if no reps
	#If there exist a sizeable number of control transcripts that should not be DE, the the
	#dispersion could be estimated from them. For example, suppose that housekeeping is
	#an index variable identifying housekeeping genes that do not respond to the treatment
	#used in the experiment. First create a copy of the data object with only one treatment
	#group:
	#	> d1 <- d
	#> d1$samples$group <- 1
	#Then estimate the dispersion from the housekeeping genes and all the libraries as one
	#group:
	#		> d0 <- estimateCommonDisp(d1[housekeeping,])
	#Then insert this into the full data object and proceed:
	#		> d$common.dispersion <- d0$common.dispersion
	#> et <- exactTest(d)
}

edgeR_DE_postanalysis = function (DGEList, de.tgw, baseout='tmp', dispersion='auto', FDR=0.05,kclusters = 10) {
	library(gdata,quietly=T,warn.conflicts=F,verbose=F)
	library(gplots,quietly=T,warn.conflicts=F,verbose=F)
	library(cluster,quietly=T,warn.conflicts=F,verbose=F)
	library(Biobase,quietly=T,warn.conflicts=F,verbose=F)
	library(rjson,quietly=T,warn.conflicts=F,verbose=F)
	library(ape,quietly=T,warn.conflicts=F,verbose=F)
	library(edgeR,quietly=T,warn.conflicts=F,verbose=F)
	library(fdrtool,quietly=T,warn.conflicts=F,verbose=F)
	
	
	
	# qvalues with fdrtool; correlates very well with FDR but more stringent; lfdr essentially bins values to compute fdr and is even more stringent
	# $table: data frame containing columns for the log2-fold-change, 'logFC',
	# the average log2-counts-per-million, 'logCPM', the two-sided p-value 'PValue'
	## so log2 of 0.1375035 and -0.1520031 is 1.10 and 0.90 FC respectively
	#
	# genes with positive log-fold change are up-regulated in group B compared
	# with group A (and vice versa for genes with negative log-fold change).
	pdf(file=paste(baseout,'.qc_graphs2.pdf',sep=''))
	logFC.box <-boxplot(de.tgw$table$logFC,main="log2 fold-change boxplot")
	logFC.dev <- abs(0-logFC.box$stats[3])
	fdrresults<-fdrtool(de.tgw$table$PValue,statistic='pvalue',plot=T)
	de.tgw$table$QValue <- fdrresults$qval
	de.tgw$table$lfdr <- fdrresults$lfdr
	qval.box <- boxplot(de.tgw$table$QValue,main="Q-values (bayesian, tail of area) boxplot")
	dev.off()
	
	# potentially housekeeping genes
	#	housekeeping<-rownames(subset(de.tgw$table, PValue >= (1-FDR)))
	# housekeeping<-rownames(subset(de.tgw$table, logFC >=logFC.box$stats[4] ))
	housekeeping<-rownames(de.tgw$table[de.tgw$table$logFC > -logFC.dev & de.tgw$table$logFC < logFC.dev ,])
	unlink(paste(baseout,'.housekeeping.txt',sep=''))
	lapply(housekeeping, write, paste(baseout,'.housekeeping.txt',sep=''), append=T, ncolumns=1,sep ="\n")
	
	
	topTagCount = sum(de.tgw$table$QValue < FDR) 
	toptgw = topTags(de.tgw, n = topTagCount);
	detags_toptgw = toptgw$table[toptgw$table$FDR < FDR,];
	detags_names = rownames(detags_toptgw)
	
	postscript(file=paste(baseout,'.ps',sep=''))
	plotSmear(DGEList, de.tags=detags_names)
	abline(h=c(-1,1), col="blue")
	dev.off()
	
	# table with entries of interest
	unlink(paste(baseout,'.results.txt',sep=''))
	write.table(de.tgw$table, quote=F, file=paste(baseout,'.all.txt',sep=''), sep="\t",col.names=NA)
	write.table(detags_toptgw, quote=F, file=paste(baseout,'.results.txt',sep=''), sep="\t",col.names=NA)
print ("Comparing groups for significant genes...")
#	ordered.data <- as.data.frame(matrix(t(DGEList$counts[rownames(DGEList$counts) %in% detags_names]),ncol=length(as.vector(DGEList$samples$description))))
#	ordered.data <- as.data.frame(matrix(t(de.tgw$CPM[rownames(de.tgw$CPM) %in% detags_names]),ncol=length(colnames(de.tgw$CPM))))
#	rownames(ordered.data) <- detags_names
#	colnames(ordered.data)<-colnames(de.tgw$CPM)
	ordered.data<-de.tgw$CPM
	sample_names <- colnames(ordered.data)
	
	# gene plots for each library
	pdf(file=paste(baseout,'.per_gene_plots.de.pdf',sep=''))
	par(mfrow=c(3, 2))
	#remember [row,column]
	for (i in 1:length(ordered.data[,1])) {
		gene_name <- rownames(ordered.data)[i]
		if (gene_name %in% detags_names){
			d <- as.numeric(ordered.data[i,])
			ymax <- max(d) 
			#if (is.infinite(min(d))){ymin=0}
			ymin <- min(d) 
			val <- formatC(de.tgw$table[gene_name,'QValue'],format='e',digits=2)
			plot(d, type="l", ylim=c(ymin,ymax), main=paste(gene_name,paste('FDR:',val,sep=' ') ,sep="\n"), col="blue", xaxt="n", xlab="", ylab="CPM")
			axis(side=1, at=1:length(d), labels=sample_names, las=0)
		}
	}  
	dev.off()
	rm(d,i,ymin,ymax)
	
	pdf(file=paste(baseout,'.per_gene_plots.housekeeping.pdf',sep=''))
	par(mfrow=c(3, 2))
	#remember [row,column]
	for (i in 1:length(ordered.data[,1])) {
		gene_name <- rownames(ordered.data)[i]
		if (gene_name %in% housekeeping){
			d <- as.numeric(ordered.data[i,])
			ymax <- max(d) 
			#if (is.infinite(min(d))){ymin=0}
			ymin <- min(d) 
			val <- formatC(de.tgw$table[gene_name,'QValue'],format='e',digits=2)
			plot(d, type="l", ylim=c(ymin,ymax), main=paste(gene_name,paste('FDR:',val,sep=' ') ,sep="\n"), col="blue", xaxt="n", xlab="", ylab="CPM")
			axis(side=1, at=1:length(d), labels=sample_names, las=0)
		}
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
print ("Clustering...")
	hc_genes <- agnes(centered.data, diss=FALSE, metric="euclidean") # cluster genes
	hc_samples <- hclust(as.dist(1-cor(centered.data, method="spearman")), method="complete") # cluster conditions
	myheatcol <- redgreen(75)
	gene_partition_assignments <- cutree(as.hclust(hc_genes), k=kclusters)
	partition_colors <- rainbow(length(unique(gene_partition_assignments)), start=0.4, end=0.95)
	gene_colors <- partition_colors[gene_partition_assignments]

	# heatmap
print ("Heatmap and other graphs")
	postscript(file=paste(baseout,'.heatmap.ps',sep=''), horizontal=FALSE, width=8, height=18, paper="special")
	figure.heatmap1<-heatmap.2(centered.data, dendrogram="both", Rowv=as.dendrogram(hc_genes), Colv=as.dendrogram(hc_samples), col=myheatcol, RowSideColors=gene_colors, scale="none", density.info="none", trace="none", key=TRUE, keysize=1.2, cexCol=2.5, margins=c(15,15), lhei=c(0.4,2), lwid=c(2.5,4))
	dev.off();
	
	## data for JS interactive map
#print ("Preparing JS")
	fdr.data<-as.data.frame(detags_toptgw$FDR)
	rownames(fdr.data)<-detags_names
	fdr.data.json<-toJSON(as.data.frame(t(fdr.data)))
	figure.heatmap1.data.json<-toJSON(as.data.frame(figure.heatmap1$carpet))
	figure.heatmap1.samples.json<-toJSON(rownames(figure.heatmap1$carpet))
	figure.heatmap1.genes.json<-toJSON(colnames(figure.heatmap1$carpet))
	write(fdr.data.json,file=paste(baseout,'.fdr.json',sep=''))
	write(figure.heatmap1.data.json,file=paste(baseout,'.data.json',sep=''))
	write(figure.heatmap1.samples.json,file=paste(baseout,'.samples.json',sep=''))
	write(figure.heatmap1.genes.json,file=paste(baseout,'.genes.json',sep=''))
	hc_genes$labels <-hc_genes$order.lab # add labels
	write.tree( as.phylo(as.hclust(hc_genes)),paste(baseout,'.hcgenes.newick',sep=''))
	write.tree( as.phylo(as.hclust(hc_samples)),paste(baseout,'.hcsamples.newick',sep=''))
	
	return (de.tgw)
}


edgeR_heatmaps_all = function(baseout='tmp',kclusters = 10){
	library(gdata,quietly=T,warn.conflicts=F,verbose=F)
	library(gplots,quietly=T,warn.conflicts=F,verbose=F)
	library(cluster,quietly=T,warn.conflicts=F,verbose=F)
	library(Biobase,quietly=T,warn.conflicts=F,verbose=F)
	library(rjson,quietly=T,warn.conflicts=F,verbose=F)
	library(ape,quietly=T,warn.conflicts=F,verbose=F)
	
	# variables
	data <- read.table(file=baseout, header=T, com="", sep="\t",row.names=1)
	
	ordered.data <- data[with(data, order(data[,length(data[1,])])), ]
	## last column is now FDR column
	gene_names <- rownames(ordered.data)
	sample_names <- colnames(ordered.data[,-length(ordered.data[1,])])
	
	# gene plots
	pdf(file=paste(baseout,'.per_gene_plots.pdf',sep=''))
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
	postscript(file=paste(baseout,'.heatmap.ps',sep=''), horizontal=FALSE, width=8, height=18, paper="special")
	figure.heatmap1<-heatmap.2(centered.data, dendrogram="both", Rowv=as.dendrogram(hc_genes), Colv=as.dendrogram(hc_samples), col=myheatcol, RowSideColors=gene_colors, scale="none", density.info="none", trace="none", key=TRUE, keysize=1.2, cexCol=2.5, margins=c(15,15), lhei=c(0.4,2), lwid=c(2.5,4))
	dev.off();
	
	## data for JS interactive map
	fdr.data<-as.data.frame(ordered.data$FDR)
	rownames(fdr.data)<-rownames(ordered.data)
	fdr.data.json<-toJSON(as.data.frame(t(fdr.data)))
	figure.heatmap1.data.json<-toJSON(as.data.frame(figure.heatmap1$carpet))
	figure.heatmap1.samples.json<-toJSON(rownames(figure.heatmap1$carpet))
	figure.heatmap1.genes.json<-toJSON(colnames(figure.heatmap1$carpet))
	cat(fdr.data.json,file=paste(baseout,'.fdr.json',sep=''))
	cat(figure.heatmap1.data.json,file=paste(baseout,'.data.json',sep=''))
	cat(figure.heatmap1.samples.json,file=paste(baseout,'.samples.json',sep=''))
	cat(figure.heatmap1.genes.json,file=paste(baseout,'.genes.json',sep=''))
	hc_genes$labels <-hc_genes$order.lab # add labels
	write.tree( as.phylo(as.hclust(hc_genes)),paste(baseout,'.hcgenes.newick',sep=''))
	write.tree( as.phylo(as.hclust(hc_samples)),paste(baseout,'.hcsamples.newick',sep=''))
}


