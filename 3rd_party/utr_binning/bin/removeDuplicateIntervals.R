#!/usr/bin/env Rscript

#This script will choose the appropriate gene model when merged reads overlap multiple genes. It also generates gene-level expression stats
message("Running R script to: \n(1) choose the appropriate gene model when merged reads overlap multiple genes\n(2) generate gene-level read counts")
message("")
message("Dependencies: dplyr (written using v1.0.0)")
message("Script written in R v3.6.2")
message("")


# Libraries + Options -----------------------------------------------------
library(dplyr, warn.conflicts = F)
options(dplyr.summarise.inform = FALSE)

# Functions ---------------------------------------------------------------
usage <- function(){
 message("Rscript --vanilla ./removeDuplicateIntervals.R <intersection_of_genes_and_merged_alignments(from bin_utr.sh)> <outputPrefix>") 
}



# Command Line Parser ------------------------------------------------------
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 2) { usage(); stop("Expected 2 arguments. Found ", length(args)) }

filePath = as.character(args[1])
if (!file.exists(filePath)) stop("Can't find file: ", filePath)
outputPrefix = as.character(args[2])



# Reading Data ------------------------------------------------------------
readsIntersectingGenes.df <- read.csv(file=filePath, header = F, sep = "\t")


# Processing Data ---------------------------------------------------------

#Check our strand info to make sure all duplicated intervals have the same strand (if not, something has gone wrong)
if (readsIntersectingGenes.df %>% 
  group_by(V4) %>% 
  summarise(duplicatesHaveDifferentStrands=length(unique(V6))!=1) %>% pull(2) %>% sum() > 0) stop("'Duplicated' interval actually has 2 different strands. Something has gone horribly wrong")


#Group by interval name and choose most appropriate (most 5') gene model (nr=nonredundant)
readsIntersectingGenes_nr.df <- readsIntersectingGenes.df %>% group_by(V4) %>% summarise(.groups="drop_last",
chrom=V1[1], 
start=V2[1],
end=V3[1],
allocated_feature_name=if(V6[1]=="+") V10[which.min(V8)] else V10[which.max(V8)],
reads_mapped=V5[1], 
strand=V6[1])
colnames(readsIntersectingGenes_nr.df)[1] <- "read_names"



#For preparation of BED file:
#Append original interval names to gene model prediction (separated by |)
readsIntersectingGenes_nr.bed.df <- readsIntersectingGenes_nr.df
readsIntersectingGenes_nr.bed.df$allocated_feature_name <- paste(readsIntersectingGenes_nr.bed.df$allocated_feature_name, readsIntersectingGenes_nr.bed.df[[1]], sep="|Reads=")
readsIntersectingGenes_nr.bed.df <- select(readsIntersectingGenes_nr.bed.df, -1)

#Note if this was a proper bed - we would lerp scores between 0 and 1000 or something along this line

#Check if we've actually allocated each merged read-set to a single gene
#nrow(readsIntersectingGenes.df)
#nrow(readsIntersectingGenes_nr.bed.df)
# message("Unique interval names:", length(unique(readsIntersectingGenes_nr.bed.df$name)))
# message("Total interval names:",length(readsIntersectingGenes_nr.bed.df$name))
message("Merged reads have been uniquely mapped to one gene: ", length(readsIntersectingGenes_nr.bed.df$allocated_feature_name[table(readsIntersectingGenes_nr.bed.df$allocated_feature_name) > 1]) == 0)

# Output data -------------------------------------------------------------
#Write feature-oriented expression summary
readsPerGene.df <- readsIntersectingGenes_nr.df %>% group_by(allocated_feature_name) %>% summarise(reads_mapped=sum(reads_mapped))

#Write Output
write.table(x = readsPerGene.df, file = paste0(outputPrefix, ".geneExpressionSummary.tsv"), append = F, quote = F, sep = "\t", row.names = F, col.names = T)
write.table(x = readsIntersectingGenes_nr.df, file = paste0(outputPrefix, ".tsv"), quote = F, row.names = F, col.names = T, append = F, sep = "\t")
write.table(x = readsIntersectingGenes_nr.bed.df, file = paste0(outputPrefix, ".bed"), quote = F, row.names = F, col.names = F, append = F, sep = "\t")

message("")
message("Output: ")
message(paste0(outputPrefix, ".tsv"), " - each merged alignments interval allocated to a single, 'most correct' (i.e. most 5') gene amongst those it overlaps (+ read count info) [TSV w/headers]" )
message( paste0(outputPrefix, ".bed"), " - each merged alignments interval allocated to a single, 'most correct' (i.e. most 5') gene amongst those it overlaps (+ read count info) [BED track (for visualisation)]" )
message(paste0(outputPrefix, ".geneExpressionSummary.tsv"), " - gene-level expression counts")

message("\nComplete")
