#!/bin/bash
echo "Installing dependencies (Perl, R)"
echo you may want to run this as root...
sleep 3
echo Installing Perl dependencies
sleep 1
cpanm Pod::Usage Getopt::Long Digest::MD5 Digest::SHA Statistics::Descriptive Time::Progress Time::localtime File::stat File::Path Text::CSV_XS Bio::SeqIO \
 Cwd JSON JSON::PP FindBin DBI Storable GD::SVG SVG Devel::CheckLib Bio::Graphics Data::Dumper Compress::LZ4
echo Installing R dependencies
sleep 1
R -e 'local({r <- getOption("repos"); r["CRAN"] <- "http://cran.csiro.au/"; options(repos=r)}); install.packages(c("cluster","gplots","rjson","fdrtool","ggplot2","doMC")) ; source("http://bioconductor.org/biocLite.R"); biocLite(c("Biobase","edgeR","ape"))'
echo Completed!

