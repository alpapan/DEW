#!/bin/bash

../bin/utrbin.pl -gff_file test.gff3 \
 -bam_files test_alignment.bam \
 -gene_fasta test.gff3.gene.fasta \
 -verbose
