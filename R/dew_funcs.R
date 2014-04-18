# 
# Author: alexie papanicolaou; some original code and inspiration from brian haas
# Complaints to alexie@butterflybase.org; CSIRO Ecosystem Sciences, GPO 1700, Canberra 2601, Australia
#        
###############################################################################

dew_install = function(){
	options("menu.graphics"=FALSE)
	bioconductor_libraries_to_install =c('edgeR','DESeq')
	base_libraries_to_install = c('Biobase','fdrtool','ggplot2','lmtest','lattice','gdata','gplots','cluster','rjson','ape')
	for (i in 1:length(base_libraries_to_install)){
		dew_install_lib(base_libraries_to_install[i])
	}
	
	for (i in 1:length(bioconductor_libraries_to_install)){
		dew_install_lib(bioconductor_libraries_to_install[i],1)
	}
}

dew_install_lib = function(lib,is_bioconductor=0){
	if(require(lib, character.only = TRUE,quietly=T,warn.conflicts=F)){
		cat(paste('Library ',lib," is loaded correctly","\n",sep=''))
	} else {
		cat(paste("Trying to install ",lib,"\n",sep=''))
		if (is.null(is_bioconductor)){
			install.packages(lib)
		}else{
			BiocInstaller::biocLite(lib)
		}
		if(require(lib, character.only = TRUE,quietly=T,warn.conflicts=F )){
			cat(paste('Library ',lib," installed and loaded","\n",sep=''))
		} else {
			stop(c("Could not install ",lib,"\n"))
		}
	}
}

convert_from_uid = function(hash,uids){
	# not in same order: x <- hash[hash[,1] %in% as.vector(uids),2]
	# validated:
	x<-c()
	for (i in 1:length(uids)){
		key <- which(hash[,1] ==  uids[i])
		x<-append(x,hash[key,2])
	}
	return(as.vector(x))
}


scale_with_edgeR = function(edgeR_obj){
	library(edgeR,quietly=T,warn.conflicts=F,verbose=F)
	logCPM_edgeR <- predFC(edgeR_obj, prior.count.total=2*ncol(edgeR_obj))
	return(as.data.frame(logCPM_edgeR))
}

scale_with_DESeq = function(x){
	library(DESeq,quietly=T,warn.conflicts=F,verbose=F)
	library(BiocGenerics,quietly=T,warn.conflicts=F,verbose=F)
	library(locfit,quietly=T,warn.conflicts=F,verbose=F)
	# provide raw counts
	cds <- newCountDataSet( round(x), rep( "dummy", ncol(x) ) )
	cds <- estimateSizeFactors( cds )
	cds <- estimateDispersions(cds,method='blind')
	vsd_stable <- getVarianceStabilizedData( cds )
	return (as.data.frame(vsd_stable))
}

check_distribution = function (edgeR_obj,baseout='tmp'){
	library(ggplot2,quietly=T,warn.conflicts=F,verbose=F)
	library(lattice,quietly=T,warn.conflicts=F,verbose=F)
	data <- edgeR_obj$CPM
	if (is.null(data)){return (FALSE)}
	groups <- colnames( data )
	#genes <- rownames(data)
	
	x_for_density<-data.frame(
			group = factor( 
					rep(groups ,each=length(rownames(data)  ))  
			),
			value = as.vector(as.matrix(data))
	);
	x_for_density_box<-boxplot(x_for_density$value,plot=F)
	
	x_for_density_log2<-data.frame(
			group = factor( 
					rep(groups ,each=length(rownames(data)  ))  
			),
			value = as.vector(log2(1+ as.matrix(data)))
	);
	x_for_density_log2_box<-boxplot(x_for_density_log2$value,plot=F)
	
#	x_for_density_DESeq<-data.frame(			group = factor( 					rep(colnames( data ) ,each=length(rownames(data)  ))  			),			value = as.vector(as.matrix(scale_with_DESeq(data)))	);
#	x_for_density_DESeq_box<-boxplot(x_for_density_DESeq$value,plot=F)
	
	x_for_density_edgeR<-data.frame(
			group = factor( 
					rep(colnames( data ) ,each=length(rownames(data)  ))  
			),
			value = as.vector(as.matrix(scale_with_edgeR(data)))
	);
	x_for_density_edgeR_box<-boxplot(x_for_density_edgeR$value,plot=F)
	
	#	x_for_density_scaled<-data.frame(
#			group = factor( 
#					rep(groups ,each=length(rownames(data)  ))  
#			),
#			value = as.vector(scale(as.matrix(data),center=F))
#	);
#	x_for_density_scaled_box<-boxplot(x_for_density_scaled$value,plot=F)
#	
#	x_for_density_scaled_clustering<-data.frame(
#			group = factor( 
#					rep(groups ,each=length(rownames(data)  ))  
#			),
#			value = as.vector(scale(as.matrix(data)))
#	);
#	x_for_density_scaled_clustering_box<-boxplot(x_for_density_scaled_clustering$value,plot=F)
	
	pdf(file=paste(baseout,'.distrib_qc_graphs_orig_log2_deseq_edgeR.pdf',sep=''))
	#
	legend_style <- guides(fill=guide_legend(ncol=1,override.aes = list(size = 0)))
	if (length(groups) > 10 ){
		legend_style <- guides(fill=guide_legend(ncol=2,override.aes = list(size = 0)))
	}
	print(ggplot(x_for_density[x_for_density$value < x_for_density_box$stats[4], ], aes(x=group, y=value, fill=group)) + geom_boxplot()  + ggtitle(wrap_text('Original data', width = 30))  + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + legend_style )
	print(ggplot(x_for_density, aes(x=value, fill=group)) + geom_density(alpha=.3)  + coord_cartesian(xlim = c(0, x_for_density_box$stats[5])) + ggtitle(wrap_text('Original data', width = 30)) + legend_style )
	
	print(ggplot(x_for_density_log2, aes(x=group, y=value, fill=group)) + geom_boxplot() + ggtitle(wrap_text('Log2 data', width = 30))  + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + legend_style )
	print(ggplot(x_for_density_log2, aes(x=value, fill=group)) + geom_density(alpha=.3) + ggtitle(wrap_text('Log 2data', width = 30)) + legend_style)
	
#	print(ggplot(x_for_density_DESeq, aes(x=group, y=value, fill=group)) + geom_boxplot() + ggtitle(wrap_text('DESeq transformed data', width = 30))  + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + legend_style )
#	print(ggplot(x_for_density_DESeq, aes(x=value, fill=group)) + geom_density(alpha=.3)  + ggtitle(wrap_text('DESeq transformed data', width = 30)) + legend_style)
	
	print(ggplot(x_for_density_edgeR, aes(x=group, y=value, fill=group)) + geom_boxplot() + ggtitle(wrap_text('edgeR transformed data', width = 30))  + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + legend_style )
	print(ggplot(x_for_density_edgeR, aes(x=value, fill=group)) + geom_density(alpha=.3) + ggtitle(wrap_text('edgeR transformed data', width = 30)) + legend_style)
	
	#print(ggplot(x_for_density_scaled, aes(x=value, fill=group)) + geom_density(alpha=.3)  + coord_cartesian(xlim = c(0, 1))
	#print(ggplot(x_for_density_scaled_clustering, aes(x=value, fill=group)) + geom_density(alpha=.3)  + coord_cartesian(xlim = c(0, 1)))
	
	dev.off()
	
	# QQ Plots
#	png(file=paste(baseout,'.distrib_qc_graphs-qq-orig-%d.png',sep=''), width = 1200, height = 800)
#	for (i in 1:(length(groups)-1)) {
#		for (k in (i+1):length(groups)) {
#			print(qq(group ~ value, aspect = 1,data=x_for_density,subset=(group == groups[i] | group == groups[k])))
#		}
#	}
#	dev.off()
#	png(file=paste(baseout,'.distrib_qc_graphs-qq-log2-%d.png',sep=''), width = 1200, height = 800)
#	for (i in 1:(length(groups)-1)) {
#		for (k in (i+1):length(groups)) {
#			print(qq(group ~ value, aspect = 1,data=x_for_density_log2,subset=(group == groups[i] | group == groups[k])))
#		}
#	}
#	dev.off()
#	png(file=paste(baseout,'.distrib_qc_graphs-qq-deseq-%d.png',sep=''), width = 1200, height = 800)
#	for (i in 1:(length(groups)-1)) {
#		for (k in (i+1):length(groups)) {
#			print(qq(group ~ value, aspect = 1,data=x_for_density_DESeq,subset=(group == groups[i] | group == groups[k])))
#		}
#	}
#	dev.off()
#	png(file=paste(baseout,'.distrib_qc_graphs-qq-edgeR-%d.png',sep=''), width = 1200, height = 800)
#	for (i in 1:(length(groups)-1)) {
#		for (k in (i+1):length(groups)) {
#			print(qq(group ~ value, aspect = 1,data=x_for_density_edgeR,subset=(group == groups[i] | group == groups[k])))
#		}
#	}
#	dev.off()	
#	
#	for (i in 1:(length(groups)-1)) {
#		for (k in (i+1):length(groups)) {
#			print(qq(group ~ value, aspect = 1,data=x_for_density_scaled,subset=(group == groups[i] | group == groups[k])))
#		}
#	}
#	
#	for (i in 1:(length(groups)-1)) {
#		for (k in (i+1):length(groups)) {
#			print(qq(group ~ value, aspect = 1,data=x_for_density_scaled_clustering,subset=(group == groups[i] | group == groups[k])))
#		}
#	}
	
}

normalize = function(x){
	if(length(which(is.na(x)))==0){
		(x-mean(x))/apply(x, 2, sd)
	}
	else {
		(x-mean(x,na.rm=T))/apply(x,2,sd,na.rm=T)
	}
}

TMM_normalize = function(targetsFile, minCPM = 2,minLibs = 1,outfile=NULL) {
	library(edgeR,quietly=T,warn.conflicts=F,verbose=F)	
	targets = read.delim(file=targetsFile, stringsAsFactors = FALSE)
	edgeR_obj = readDGE(targets)
	# how many factors are we going to look at ?
	edgeR_obj$factors_to_do<-c()
	for(i in colnames(edgeR_obj$samples)) { 
		if (i %in% c("files","description","lib.size","norm.factors","eff.lib.size", "group")){next}; 
		edgeR_obj$factors_to_do<-append(edgeR_obj$factors_to_do,i)
	}
	
	# only use those rows that have a minimum number of tags (min to be DE)
	edgeR_obj$CPM <- cpm(edgeR_obj,normalized.lib.sizes=F)
	edgeR_obj<-edgeR_obj[rowSums(edgeR_obj$CPM >= minCPM) >= minLibs,]
	# reset library sizes
	edgeR_obj$samples$lib.size = colSums(edgeR_obj$counts)
	
	# run TMM normalization
	edgeR_obj = calcNormFactors(edgeR_obj)
	rownames(edgeR_obj$samples)<-edgeR_obj$samples$description
	colnames(edgeR_obj$counts)<-edgeR_obj$samples$description
	
	edgeR_obj$samples$eff.lib.size = edgeR_obj$samples$lib.size * edgeR_obj$samples$norm.factors
	# normalized CPM
	edgeR_obj$CPM <- cpm(edgeR_obj,normalized.lib.sizes=T)
	edgeR_obj$logCPM <- log2(1+edgeR_obj$CPM)
	
	if (!is.null(outfile)){
		unlink(outfile)
		check_distribution(edgeR_obj,outfile)
		write.table(edgeR_obj$samples, file=outfile, quote=F, sep="\t", row.names=F)
	}
	return (edgeR_obj) 
}


# TODO split between those with reps and without reps
edgeR_DE_explore = function (genesaliases_file,targetsFile, baseout='tmp', dispersion='auto', FDR=0.05,kclusters = 10,minCPM = 2,minLibs = 1){
	#genesaliases_file='/databases/owner/pap056/rna/wei_only_results/wei_only.toalign.md5';targetsFile='/databases/owner/pap056/rna/wei_only_results//edgeR/HaGR01_vs_HaGR16.edgeR_target.files'; baseout='tmp'; dispersion=0.4; FDR=0.05;kclusters = 10;minCPM = 2;minLibs = 1
	
	library(Biobase,quietly=T,warn.conflicts=F,verbose=F)
	library(edgeR,quietly=T,warn.conflicts=F,verbose=F)
	library(fdrtool,quietly=T,warn.conflicts=F,verbose=F)
	
	aliases <- as.matrix(read.table(file=genesaliases_file,sep="\t",header=F,col.names=c('uid','alias')))
	
	edgeR_obj = TMM_normalize(targetsFile,minCPM,minLibs )
	if (length(edgeR_obj)==0){return}
	file_number <- length(edgeR_obj$samples$files)	
	if (length(edgeR_obj$factors_to_do) == 0){
		if (file_number > 2 ){
			pdf(file=paste(baseout,'.qc_graphs.pdf',sep=''))
			plotMDS(edgeR_obj, main="Biological Co-efficient of Variation distance MDS plot");
			plotMDS(edgeR_obj$logCPM, main="Counts per million distance MDS plot")
		}
		if (dispersion == 'auto'){
			edgeR_obj <- estimateCommonDisp(edgeR_obj)
			edgeR_obj <- estimateTrendedDisp(edgeR_obj)
			edgeR_obj <- estimateTagwiseDisp(edgeR_obj)
			plotBCV(edgeR_obj, main="Dispersion trends")
			# set dispersion to zero for Poisson-equivalent results.	
		}
		if (file_number > 2 ) dev.off()
		
		cat ("Processing one factor using CML\n");
		edgeR_obj_de = exactTest(edgeR_obj, dispersion=dispersion)
		
		# pvalues
		pdf(file=paste(baseout,'.qc_graphs2.pdf',sep=''))
		logFC.box <-boxplot(edgeR_obj_de$table$logFC,main="log2 fold-change boxplot")
		logFC.dev <- abs(0-logFC.box$stats[3])
		fdrresults<-fdrtool(edgeR_obj_de$table$PValue,statistic='pvalue',plot=T)
		edgeR_obj_de$table$QValue <- fdrresults$qval
		edgeR_obj_de$table$lfdr <- fdrresults$lfdr
		edgeR_obj_de$logFC.dev <- logFC.dev 
		qval.box <- boxplot(edgeR_obj_de$table$QValue,main="Q-values (bayesian, tail of area) boxplot")
		dev.off()
		
		edgeR_obj_de <- edgeR_DE_postanalysis(aliases,edgeR_obj, edgeR_obj_de, baseout, dispersion, FDR,kclusters)
	}else{
		cat ("Processing multiple factors using additive GLMs\n");
		pdf(file=paste(baseout,'.qc_graphs.pdf',sep=''))
		plotMDS(edgeR_obj, main="Biological Co-efficient of Variation distance MDS plot");
		plotMDS(edgeR_obj$logCPM, main="Counts per million distance MDS plot")
		group <- as.factor(edgeR_obj$samples[, 'group'])
		# many factors, do glm
		form <- as.formula(
				paste("~ group + ", paste(
								paste('as.factor(edgeR_obj$samples[, edgeR_obj$factors_to_do[', 1:length(edgeR_obj$factors_to_do),']])', sep="")
								, collapse= "+")
				)
		)
		design<-model.matrix(form)
		rownames(design) <- colnames(edgeR_obj)
		edgeR_obj <- estimateGLMCommonDisp(edgeR_obj, design, verbose=TRUE)
		edgeR_obj <- estimateGLMTrendedDisp(edgeR_obj, design)
		edgeR_obj <- estimateGLMTagwiseDisp.edgeR_obj(edgeR_obj, design=design) 
		plotBCV(edgeR_obj, main="Dispersion trends")
		dev.off()
		
		fit <- glmFit(edgeR_obj, design)
		if (length(fit$fail) >0){return} 
		
		# Explore the effect of every factor
		#coeff 1 is intercept;
		# coeff 2 is group1 (starts with 'group')
		# then the other coefficients could be other groups or other factors.
		coefficients <-colnames(fit$coefficients)
		positives <-c()
		names_positives<-c()
		# we will use this later
		edgeR_obj_de <- glmLRT(fit, coef=2)
		positives <-append(positives,sum(p.adjust(edgeR_obj_de$table$PValue, method="BH")< FDR))
		names_positives<- append(names_positives,(paste('Coefficient group vs others:',sep="")))
		
		for (i in 3:length(coefficients)){
			edgeR_obj_de.t <- glmLRT(fit, coef=i)
			positives <-append(positives,sum(p.adjust(edgeR_obj_de.t$table$PValue, method="BH")< FDR))
			names_positives<- append(names_positives,(paste("Coeffs ",i,' vs others:',sep="")))
		}
		if (length(coefficients) > 3){
			for (i in 3:length(coefficients)){
				for (k in 4:length(coefficients)){
					if (k == i){next}
					edgeR_obj_de.t <- glmLRT(fit, coef=i:k)
					positives <-append(positives,sum(p.adjust(edgeR_obj_de.t$table$PValue, method="BH")< FDR))
					names_positives <- append(names_positives,(paste("Coeffs ",i,':',k,' vs others:',sep="")))
				}
			}
		}
		rm(edgeR_obj_de.t)
		positives <-as.data.frame(positives)
		rownames(positives)<-names_positives
		write.table(positives,row.names=T,quote=F,sep="\t",file=paste(baseout,".factors.positives",sep=""))
		write(c("factors:"),file=paste(baseout,".factors",sep=""))
		write.table(edgeR_obj$factors_to_do,col.names=F,row.names=T,quote=F,sep="\t",file=paste(baseout,".factors",sep=""),append=T)
		write(c("\ncoefficients:"),file=paste(baseout,".factors",sep=""),append=T)
		write.table(coefficients,col.names=F,row.names=T,quote=F,sep="\t",file=paste(baseout,".factors",sep=""),append=T)
		write(c("\nTests:"),file=paste(baseout,".factors",sep=""),append=T)
		write.table(col.names=F,quote=F,row.names=F,length(edgeR_obj_de$table$PValue),file=paste(baseout,".factors",sep=""),append=T)
		write(c("\nPositives past ",FDR,':'),file=paste(baseout,".factors",sep=""),append=T)
		write.table(positives,col.names=F,row.names=T,quote=F,sep="\t",file=paste(baseout,".factors",sep=""),append=T)
		
		# pvalues
		pdf(file=paste(baseout,'.qc_graphs2.pdf',sep=''))
		logFC.box <-boxplot(edgeR_obj_de$table$logFC,main="log2 fold-change boxplot")
		logFC.dev <- abs(0-logFC.box$stats[3])
		fdrresults<-fdrtool(edgeR_obj_de$table$PValue,statistic='pvalue',plot=T)
		edgeR_obj_de$table$QValue <- fdrresults$qval
		edgeR_obj_de$table$lfdr <- fdrresults$lfdr
		edgeR_obj_de$logFC.dev <- logFC.dev
		qval.box <- boxplot(edgeR_obj_de$table$QValue,main="Q-values (bayesian, tail of area) boxplot")
		dev.off()
		
		# give us the results for the group factor
		edgeR_obj_de <- edgeR_DE_postanalysis(aliases,edgeR_obj, edgeR_obj_de, baseout, dispersion, FDR,kclusters)
	}
	save.image(file=paste(baseout,".edgeR.Rdata",sep=""));
#	return (edgeR_obj_de)
}


edgeR_DE_analysis_no_reps_housekeeping = function (edgeR_obj, baseout='tmp', dispersion='auto', FDR=0.05,kclusters = 10) {
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


edgeR_DE_postanalysis_create_housekeeping = function (aliases,edgeR_obj_de,baseout='tmp'){
	housekeeping<-as.vector(rownames(edgeR_obj_de$table[edgeR_obj_de$table$logFC > -edgeR_obj_de$logFC.dev & edgeR_obj_de$table$logFC < edgeR_obj_de$logFC.dev ,]))
	unlink(paste(baseout,'.housekeeping.txt',sep=''))
	lapply(housekeeping, write, paste(baseout,'.housekeeping.txt',sep=''), append=T, ncolumns=1,sep ="\n")
	housekeeping_alias <- convert_from_uid(aliases,housekeeping)
	unlink(paste(baseout,'.housekeeping.txt.genes',sep=''))
	lapply(housekeeping_alias, write, paste(baseout,'.housekeeping.txt.genes',sep=''), append=T, ncolumns=1,sep ="\n")
	
	# TODO ensure expressed; add minimal dispersion as a requirement too (near common?)
	return(housekeeping)
}

edgeR_DE_postanalysis = function (aliases,edgeR_obj, edgeR_obj_de, baseout='tmp', dispersion='auto', FDR=0.05,kclusters = 10) {
	library(gdata,quietly=T,warn.conflicts=F,verbose=F); 
	library(gplots,quietly=T,warn.conflicts=F,verbose=F)
	library(cluster,quietly=T,warn.conflicts=F,verbose=F)
	library(rjson,quietly=T,warn.conflicts=F,verbose=F)
	library(ape,quietly=T,warn.conflicts=F,verbose=F)
	
	# potentially housekeeping genes
	housekeeping<-edgeR_DE_postanalysis_create_housekeeping(aliases,edgeR_obj_de,baseout)
	
	
	topTagCount = sum(edgeR_obj_de$table$QValue < FDR) 
	toptgw = topTags(edgeR_obj_de, n = topTagCount);
	detags_toptgw = toptgw$table[toptgw$table$FDR < FDR,];
	detags_names = rownames(detags_toptgw)
	#detags_alias <- convert_from_uid(aliases,detags_names) 
	
	
	postscript(file=paste(baseout,'.ps',sep=''))
	plotSmear(edgeR_obj, de.tags=detags_names)
	abline(h=c(-1,1), col="blue")
	dev.off()
	
	# table with entries of interest
	unlink(paste(baseout,'.results.txt',sep=''))
	write.table(edgeR_obj_de$table, quote=F, file=paste(baseout,'.all.txt',sep=''), sep="\t",col.names=NA)
	write.table(detags_toptgw, quote=F, file=paste(baseout,'.results.txt',sep=''), sep="\t",col.names=NA)
	cat ("Comparing groups for significant genes...\n")
	logCPM_edgeR <- predFC(edgeR_obj, prior.count.total=2*ncol(edgeR_obj))
	ordered.data <-matrix(ncol=length(colnames(logCPM_edgeR)),nrow=length(detags_names))
	for (i in 1:length(detags_names)){
		ordered.data[i,] <- logCPM_edgeR[detags_names[i],]
	}
	rownames(ordered.data) <- detags_names
	colnames(ordered.data)<-colnames(logCPM_edgeR)
	sample_names <- colnames(ordered.data)
	
	# gene plots for each library
	pdf(file=paste(baseout,'.per_gene_plots.de.pdf',sep=''))
	par(mfrow=c(3, 2))
	#remember [row,column]
	for (i in 1:length(ordered.data[,1])) {
		gene_name <- rownames(ordered.data)[i]
		if (gene_name %in% detags_names){
			gene_alias <- aliases[gene_name,]
			d <- as.numeric(ordered.data[i,])
			ymax <- max(d) 
			#if (is.infinite(min(d))){ymin=0}
			ymin <- min(d) 
			val <- formatC(edgeR_obj_de$table[gene_name,'QValue'],format='e',digits=2)
			plot(d, type="l", ylim=c(ymin,ymax), main=paste(gene_alias,paste('FDR:',val,sep=' ') ,sep="\n"), col="blue", xaxt="n", xlab="", ylab="CPM")
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
			gene_alias <- aliases[gene_name,]
			d <- as.numeric(ordered.data[i,])
			ymax <- max(d) 
			#if (is.infinite(min(d))){ymin=0}
			ymin <- min(d) 
			val <- formatC(edgeR_obj_de$table[gene_name,'QValue'],format='e',digits=2)
			plot(d, type="l", ylim=c(ymin,ymax), main=paste(gene_alias,paste('FDR:',val,sep=' ') ,sep="\n"), col="blue", xaxt="n", xlab="", ylab="CPM")
			axis(side=1, at=1:length(d), labels=sample_names, las=0)
		}
	}  
	dev.off()
	rm(d,i,ymin,ymax)
	# add real gene names:
	clustering.data <- ordered.data
	rownames(clustering.data)<-convert_from_uid(aliases,rownames(clustering.data))
	#clustering.data <- scale_with_DESeq(ordered.data) #deseq normalize getVarianceStabilizedData; suspiciously similar to logged.data; check with Goldfeld–Quandt test gqtest
	
	# cluster data
	cat ("Clustering...\n")
	#morgan: 8-9gb; 3.5h for arabidopsis
	#library(hopach,quietly=T,warn.conflicts=F,verbose=F)
	#	library(ctc,quietly=T,warn.conflicts=F,verbose=F)
	#hc_gene_dist <- distancematrix(clustering.data, d = "cosangle")
	#hc_gene_hopach <- hopach(clustering.data, dmat = hc_gene_dist)
	#hc_samples3full<-hopach(t(clustering.data), d = "euclid")	
	#tree_hc_gene_hopach<-hopach2tree(clustering.data,
	#	file=paste(baseout,'.hopach.tree',sep=''),hopach.genes=hc_gene_hopach,dist.genes =hc_gene_dist,
	#	gene.names = rownames(clustering.data),hopach.arrays = hc_samples3full
	#)
	#tree_hc_gene_r<-xcluster2r(file=paste(baseout,'.hopach.tree.gtr',sep='')
	
	
	#morgan: 2.3 gb; 4h  
	myheatcol <- redgreen(75)
	hc_genes <- agnes(clustering.data, metric="euclidean") # cluster genes
	hc_samples <- hclust(as.dist(1-cor(clustering.data, method="spearman")), method="complete") # cluster conditions; 
	gene_partition_assignments <- cutree(as.hclust(hc_genes), k=kclusters)
	partition_colors <- rainbow(length(unique(gene_partition_assignments)), start=0.4, end=0.95)
	gene_colors <- partition_colors[gene_partition_assignments]
	
	# heatmap
	cat ("Producing info for heatmap and other graphs...\n")
	postscript(file=paste(baseout,'.heatmap.ps',sep=''), horizontal=FALSE, width=8, height=18, paper="special")
	figure.heatmap1<-heatmap.2(clustering.data, dendrogram="both", Rowv=as.dendrogram(hc_genes), Colv=as.dendrogram(hc_samples), col=myheatcol, RowSideColors=gene_colors,
			scale="none", density.info="none", trace="none", key=TRUE, keysize=1.2, cexCol=2.5, margins=c(15,15), lhei=c(0.4,2), lwid=c(2.5,4))
	dev.off();
	
	## data for JS interactive map
	# detags_toptgw and clustering.data are synced
	fdr.data<-as.data.frame(detags_toptgw$FDR)
	rownames(fdr.data)<-rownames(clustering.data)
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
	
	return (edgeR_obj_de)
}


edgeR_gene_plots_all = function(genesaliases_file,baseout='tmp'){
	require(ggplot2)
	# variables
	aliases <- as.matrix(read.table(file=genesaliases_file,sep="\t",header=F,col.names=c('uid','alias')))
	data <- read.table(file=baseout, header=T, com="", sep="\t",row.names=1)
	gene_names <- convert_from_uid(aliases,rownames(data))
	rownames(data)<-gene_names
	
	sample_names <- colnames(data)
	ordered.data <- data[with(data, order(rownames(data))), ]
	
	# transform data
	#clustering.data <- as.matrix(log2(1+ordered.data)) 
	
	# gene plots
	stacked_data <-stack(ordered.data)
	main_graph <- ggplot(stacked_data, aes(x = ind, y = values))
	main_graph <- main_graph + coord_flip()
	main_graph <- main_graph + stat_boxplot(geom ='errorbar',na.rm=T,alpha=0.3)
	main_graph <- main_graph + geom_boxplot(outlier.size = 1,na.rm=T,color='grey',outlier.colour='grey')
	main_graph <- main_graph + scale_y_log10(name=quote('log10( TMM-normalized FPKM from eff. counts)'))
	main_graph <- main_graph + scale_x_discrete(name=quote('Library'))
	main_graph <- main_graph + theme(axis.title.y = element_text(size=17, colour = rgb(0,0,0)))
	main_graph <- main_graph + theme(axis.title.x = element_text(size=10, colour=rgb(0,0,0)))
	main_graph <- main_graph + theme(axis.text.x = element_text(size=12, colour=rgb(0,0,0)))
	main_graph <- main_graph + theme(axis.text.y = element_text(size=12, colour = rgb(0,0,0),vjust=0.5))
	main_graph <- main_graph + theme(plot.title = element_text(lineheight=0.8, face=quote(bold)))

	pdf(file=paste(baseout,'.per_gene_plots.pdf',sep=''))
	if (length(sample_names) < 20){
		par(mfrow=c(1, 2),oma=c(5,0.1,1,0.1))
	}else{
		par(mfrow=c(1, 1),oma=c(5,0.1,1,0.1))
	}
	for (i in 1:length(ordered.data[,1])) {
		d<-data.frame(val=as.numeric(ordered.data[i,]),cat=as.factor(colnames(ordered.data[i,])))
		#ymin <- min(d$cat);
		#ymax <- max(d$cat);
		print ( main_graph + ggtitle(row.names(ordered.data[i,])) + geom_point(data=d,aes(x=cat,y=val),colour = "red",size=3) )
		#normal plot; it's a line plot so counter-intuitive
		#d <- ordered.data[i,]
		#plot(as.numeric(d), type="l", ylim=c(ymin,ymax), main=row.names(ordered.data[i,]),col="blue", xaxt="n", xlab="", ylab="TMM-norm. eff. counts")
		#axis(side=1, at=1:length(d), labels=sample_names, las=2,cex.axis=0.5)
	}  
	dev.off()
}

edgeR_differential_expression = function(genesaliases_file,baseout='tmp',kclusters = 10){
	library(gdata,quietly=T,warn.conflicts=F,verbose=F)
	library(gplots,quietly=T,warn.conflicts=F,verbose=F)
	library(cluster,quietly=T,warn.conflicts=F,verbose=F)
	library(Biobase,quietly=T,warn.conflicts=F,verbose=F)
	library(rjson,quietly=T,warn.conflicts=F,verbose=F)
	library(ape,quietly=T,warn.conflicts=F,verbose=F)
	
	# variables
	aliases <- as.matrix(read.table(file=genesaliases_file,sep="\t",header=F,col.names=c('uid','alias')))
	## last column is now FDR column
	data <- read.table(file=baseout, header=T, com="", sep="\t",row.names=1)
	ordered.data <- data[with(data, order(data[,length(data[1,])])), ]
	#save FDR column
	fdr.data <- ordered.data$FDR
	#remove FDR column
	ordered.data <- data[,1:(length(data[1,])-1)]
	gene_names <- convert_from_uid(aliases,rownames(ordered.data))
	rownames(ordered.data)<-gene_names
	sample_names <- colnames(ordered.data)
	
	
	# gene plots
# temporarely disable
#	pdf(file=paste(baseout,'.per_gene_plots.pdf',sep=''))
#	if (length(sample_names) < 20){
#		par(mfrow=c(1, 2),oma=c(5,0.1,1,0.1))
#	}else{
#		par(mfrow=c(1, 1),oma=c(5,0.1,1,0.1))
#	}
#	#remember [row,column]
#	for (i in 1:length(ordered.data[,1])) {
#		d <- ordered.data[i,]
#		fdr_value <- formatC(fdr.data[i],format='e',digits=2)
#		ymin <- min(d);
#		ymax <- max(d);
#		plot(as.numeric(d), type="l", ylim=c(ymin,ymax), main=paste(gene_names[i],paste('FDR:',fdr_value,sep=' ') ,sep="\n"), col="blue", xaxt="n", xlab="", ylab="TMM-norm. eff. counts")
#		axis(side=1, at=1:length(d), labels=sample_names, las=2,cex.axis=0.5)
#	}  
#	dev.off()
#	rm(d,i,ymin,ymax)
	
	# transform data
	clustering.data <- as.matrix(log2(1+ordered.data)) 
	
	
	# cluster data
	hc_genes <- agnes(clustering.data, diss=FALSE, metric="euclidean") # cluster genes
	hc_samples <- hclust(as.dist(1-cor(clustering.data, method="spearman")), method="complete") # cluster conditions  
	myheatcol <- redgreen(75)
	gene_partition_assignments <- cutree(as.hclust(hc_genes), k=kclusters)
	partition_colors <- rainbow(length(unique(gene_partition_assignments)), start=0.4, end=0.95)
	gene_colors <- partition_colors[gene_partition_assignments]
	
	# heatmap
	postscript(file=paste(baseout,'.heatmap.ps',sep=''), horizontal=FALSE, width=8, height=18, paper="special")
	figure.heatmap1<-heatmap.2(clustering.data, dendrogram="both", Rowv=as.dendrogram(hc_genes), Colv=as.dendrogram(hc_samples), col=myheatcol, RowSideColors=gene_colors, scale="none", density.info="none", trace="none", key=TRUE, keysize=1.2, cexCol=2.5, margins=c(15,15), lhei=c(0.4,2), lwid=c(2.5,4))
	dev.off();
	
	## data for JS interactive map
	fdr.data<-as.data.frame(fdr.data)
	rownames(fdr.data)<-gene_names
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

#NOTES
# qvalues with fdrtool; correlates very well with FDR but more stringent; lfdr essentially bins values to compute fdr and is even more stringent
# $table: data frame containing columns for the log2-fold-change, 'logFC',
# the average log2-counts-per-million, 'logCPM', the two-sided p-value 'PValue'
## so log2 of 0.1375035 and -0.1520031 is 1.10 and 0.90 FC respectively
#
# genes with positive log-fold change are up-regulated in group B compared
# with group A (and vice versa for genes with negative log-fold change).


wrap_text = function(x='', width = 20){
	return (paste(strwrap(width=width,x=x), collapse = "\n"))
} 
