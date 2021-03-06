#!/usr/bin/env perl

=pod

=head1 TODO

Speed up graph making by using memory (parallelizing doesn't work)

 Version 2
 * Network analysis
 * add DESeq and make Venn diagram


=head1 NAME

 dew - Differential Expression on the Web

=head1 DESCRIPTION

 requires salmon v1.4.0
 
 It accepts a reference gene FASTA and aligns readsets against it using single-end mode of BWA.
 Since alignments take the longest, it sacrifices speed for low memory usage (uses sqlite) so it is good weak for web-server with a few GB of memory or very large datasets (TODO: farm out alignments)
 User can also provide sequences (as FASTA) to be used as a reference from the cmdline
 Makes a PNG image of the coverage of a reference gene(s) for each readset provided
 
 Files are recalculated only if non-existent/older
 if sequence exists in database with alignments, then don't re-align. use md5 on sequence 
 created server to run these things on
 Added kangade
 added fpkm via express - removed it in favour of salmon
 added bowtie2; made it the default; accepts fastq.bz2
 when not in contextual alignment, salmon processes each reference sequence separately (yes, it is incorrect but db results would be incorrect otherwise and are getting the fpkm...)
 added edgeR
 added html5 interactive heatmap and scatterplot
 housekeeping genes are found automatically using pvalue 1 +/-0.05 
 Currently, highest memory use is from samtools sort (2Gb) and the memory used by depth calculations (9gb for fungal dataset). 
 SQLite is fast because it resides in memory and then written out in file; but this prevents parallel runs.

 Sort requires up to 20Gb of RAM but is configurable

 samtools I use is v. 0.1.18

=head1 EXAMPLE CMD

See /demo and try:

dew.pl -in schizosaccharomyces_pombe_972h-_2_genes.fasta.renamed.nr98 -lib lib_alias.txt \
 -1 Sp.plat.1M.left.fq.bz2 Sp.ds.1M.left.fq.bz2 Sp.hs.1M.left.fq.bz2 Sp.log.1M.left.fq.bz2 \
 -2 Sp.plat.1M.right.fq.bz2 Sp.ds.1M.right.fq.bz2 Sp.hs.1M.right.fq.bz2 Sp.log.1M.right.fq.bz2 \
 -threads 6 -contextual -correct_bias -over -uid dew_fungal_demo 2>&1 | tee dew_fungal_demo.dew.log
 
 If there is an error, email us a compressed dew_fungal_demo.dew.log 
 
test times:
 3 fungal datasets all genes;4 threads;direct I/O 
 
 for standalone (new SQL database):
  preliminary < 1'; 'each readset: 11.5' for aln (1.5')+express (9')+depth (<1') each readset; kangade=failed?;stats n graphs=18'; R=failed'; TOTAL=55'
  real    55m43.152s
  user    68m50.390s
  sys     7m18.750s

 for standalone (existing SQL database):

 for context (new SQL database):
  preliminary < 1'; 4.5' each readset: for aln (1.5')+express (2')+depth (<1') each readset; kangade=1.5';stats n graphs=18'; R=<1'; TOTAL=33:15min. 
  real    33m29.720s
  user    46m22.530s
  sys     1m16.300s
  DB: size = dew_webserver.db 992K
  
 for context (existing SQL database):
  real    31m51.256s
  user    43m29.010s
  sys     1m14.580s
  DB: size = dew_webserver.db 3.4M : why?

=head1 USAGE

 Legend (:s => string; :i => integer; :f => float; {1,} => one or more); shortened names are valid if they are unique (i.e. -1 instead of -1read) 

 Mandatory:

            -infile :s              => Reference file of genes
            -sequence :s            => Instead of a file, provide a single sequence (in FASTA format with \n as line separator);
            -format :s              => The reads can be in BAM or FASTQ format. FASTQ can be .bz2 if bowtie2 is used
            -1read|readset1|r :s{1,}=> Sets of files (one per library). Tested with Phred33 FASTQ format
            -2read|readset2: s{1,}  => If provided, do paired end alignment. Sets of paired right files (synced to readset1). Optional.

 File paths if not in $PATH:

            -samtools_exec :s       => Executable to samtools if not in your path
            -bwa_exec :s            => Executable to BWA if not in your path
            -bowtie2_exec :s        => Executable to Bowtie2 if not in your path
            -bamtools_exec :s       => Executable to bamtools if not in your path
            -kangade_exec :s        => Executable to biokanga if not in your path

 Optional:

            -output or -uid    :s   => A uid for naming output files. Optional, otherwise generate
            -threads :i             => Number of CPUs to use for alignment. BWA has no advantage over 4 threads
            -library_name_file :s   => An tag value tab delimited file for giving a friendly alias for each readset library. Needs a header line to describe columns ("file" and "name" in that order). Only include -1read files.
	    -sample_names    :s{,}  => A list of names to assign to the samples. Only if a -library_name_file is not provided. Must be same order/number as -readsets
            -need_all_readsets      => All sets of reads must have alignments against the gene in order for it to be processed. Otherwise, 1+ is sufficient. 
            -over                   => Allow overwriting of any files with the same name
            -nographs               => Do not produce any graphs. Graphs can take a very long time when there are many readsets (e.g. 30+ libraries and 30k+ genes).
            -gene_coverage_graph    => Only do enough work to get the gene depth/coverage graphs and then exit. No expression boxplots
            -no_js_graphs           => If producing edgeR graphs, then don't produce javascript based graphs.
            -contextual             => Complete realignment of all genes in order to run a correction of biases properly. Does not read/store data in the cache
            -use_bwa                => Use BWA instead of Bowtie2. Much slower.
            -prepare_only           => Quit after post-processing readsets and writing initial DB
            -seeddb :s              => Initialize database using this database file (e.g. after -prepare_only)
            -kanga                  => Use kanga instead of bowtie2 for alignments. It requires a LOT of memory (ca. 1Gb per million reads) and post-processing paired-end is much slower than bowtie  
            -existing_aln :s{1,}    => Use an existing bam file instead of doing a new alignment (must be read name sorted)
            -coord_existing         => Above existing alignments are co-ordinate sorted rather than read name sorted
            -resume                 => Load existing data from database and do not reprocess existing readsets (good for adding new readsets even with contextual. NB assumes same FASTA input so DO NOT use if you changed the FASTA reference gene file)
            -do_kangade|dokangade   => Use kangade to process pairwise libraries. Experimental
            -db_use_file            => Use file for SQL database rather than system memory (much slower but possible to analyze larger datasets)
            -dispersion             => For edgeR: if we have replicates, dispersion can be set to auto. otherwise set it to a float such as 0.1 (def)
            -fdr_cutof              => Cut off FDR for DE as float (0.001 def.)
            -cpm_cutoff             => Cut off counts-per-million for parsing alignments to edgeR (def 2)
            -library_cutoff         => Number of libraries that need to have cpm_cutoff to keep alignment (def 1)
            -binary_min_coverage :f => Minimum gene (length) coverage (0.00 ~ 1.00) for counting a gene as a expressed or not (def. 0.3)
            -binary_min_reads :i    => Minimum number of aligned reads for counting a gene as a expressed or not (def. 4)
            -never_skip             => Always process graph for gene/readset pair, regardless of cutoffs
            -sort_memory_gb :i      => Max memory (in Gb) to use for sorting (def 10). Make sure this is less than the available/requested RAM
            -remove_redund          => Remove redundant sequences (identical) in -infile
            -only_alignments        => Stop after all alignments are completed. Good for large data/alignments and HPC environments. Use without -contextual (and use with -nographs). 
            -cleanup                => Delete alignments after successful completion
            -no_pairwise            => Do not do pairwise comparisons (kangade and edgeR). Otherwise, this can be VERY memory intense for genomewide for many (20+) libraries (160Gb)
            -no_check               => When re-starting, do not check database if every gene has been stored. Do not use if you're adding new genes or database was incomplete (it will crash later), but use if you're restarting and have lots of genes.
            -verbose                => Print on the screen any system commands that are run. Caution, that will produce a lot of output on the screen they are kept in the .log file anyway).
            -no_pdf                 => Do not convert gene coverage/expression images to multi-page PDF. Otherwise, will print a PDF for every 500 genes per PDF (slow for large genomes & dozens of readsets)
	    -sort_tmp :s	    => Temporary directory to use for sorting files. It will need to be on a fast disk that has sufficient free space (depends on number of -threads)

 Salmon options:

            -isoforms               => Use salmon to correct Illumina sequencing biases and transcript isofrm assignments. Increases runtime. Use -contextual for accuracy 
            -genomewide             => Your input provides all the genes of the genome, i.e. expecting to have all reads in the readset aligning. This influences salmon only. Probably needed for genomewide analyses that have readsets with large amount of non coding sequence (e.g. rDNA). Also stores data in database cache
            -extra_options :s       => Extra options for e.g. salmon, exclude any initial --dashes from the salmon options (eg give as "salmon:minAssignedFrags 5;salmon:noLengthCorrection;salmon:libType SF", including the "quotes" ). I highly recommend you include the salmon:libType
            -readset_separation     => Expected insert size for Salmon (approximately). Defaults to 500.
	    -is_quantseq	    => Tell the program that you've sequenced only the end(s) of the genes


NB: I highly recommend you use either the latest version of Bowtie (2.1.0+) or Kanga for alignments (-kanga). I often had issues with Bowtie2 (2.0.5-) on NFS and High Perfomance Computing... Kanga is faster than Bowtie2 but my post-processing makes the whole process a lot slower
        
=head1 FAQ on errors

1. I get this error:
  $VAR1 = bless( {}, 'DBI::st' );
  Cannot find md5 data for Msex2.02821.2

 You are trying to reuse an existing directory but not using an existing SQLite database. This error happens because the file *.md5 *.checked and *.depth.completed exists in your result directory. Please delete them and restart (rm -f outdir/*md5 outdir/*checked). Do not set -resume

=head1 AUTHORS

 Alexie Papanicolaou

        CSIRO Ecosystem Sciences
        alexie@butterflybase.org

 Many thanks/credit to B. Haas (The Broad) for code bundled in Trinity-RNA-Seq

=head1 DISCLAIMER & LICENSE

Copyright 2013-2014 the Commonwealth Scientific and Industrial Research Organization. 
Copyright 2015-2016 the Western Sydney University
This software is released under the Mozilla Public License v.2.

It is provided "as is" without warranty of any kind.
You can find the terms and conditions at http://www.mozilla.org/MPL/2.0.


=head1 BUGS & LIMITATIONS

Making of BioPerl graphs has a memory leak somewhere...

See TODO, otherwise none other known so far but probably lots

=cut

use strict;
use warnings;
use Carp;
use FindBin qw/$RealBin/;
use lib $RealBin. '/PerlLib';
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use Digest::MD5 qw(md5_hex);

#use Digest::SHA qw/sha1 md5_hex/;
use Statistics::Descriptive;
use Time::Progress;
use Time::localtime;
use File::Basename;
use File::stat;
use File::Path qw(remove_tree);
use File::Copy;
use Text::CSV_XS;
use Fasta_reader;
use Cwd qw(abs_path getcwd);
use JSON::PP;
use Compress::LZ4;

#threads
use Thread_helper;    # threads -1
use threads::shared;

#db
use DBI qw(:sql_types);
use Storable qw(freeze thaw);

#graphics
use SVG;
use GD::SVG;          ## we don't want PNG anymore...
use Bio::Graphics;
use Bio::SeqFeatureI;
use Bio::SeqFeature::Generic;
use Bio::SeqFeature::Lite;
$| = 1;

$ENV{'PATH'} = "$RealBin:$RealBin/3rd_party/bin:$RealBin/util:$RealBin/3rd_party/salmon/bin:" . $ENV{'PATH'};
#################################
my (
     $debug,                   $input_reference_file,
     @readsets,                @housekeeping_ids,
     %housekeeping,            $reference_file_md5sum,
     @reference_sequence_list, %library_aliases,
     $contextual_alignment,    %library_metadata,
     $extra_genes,             %paired_readset_lookup,
     $perform_bias_correction, %fold_changes,
     %skipped_references,      %user_alias,
     %groups_readsets,         $cleanup,
     $md5_aliases_file,        $remove_redund,
     $kangade_exec,            $kangax_exec,
     $existing_bam_is_coordsorted,
     $kanga_exec, %gene_aliases_threads, $do_galaxy_cleanup, @sample_names
);

my $db_hostname = 'localhost';
my $sqlite3_exec =&check_program('sqlite3');
my (
     $ps2pdf_exec,  $inkscape_exec,
     $convert_imagemagick_exec, $pdfcrop_exec, $biokanga_exec,
     $samtools_exec,            $bowtie2_exec, $bwa_exec,
     $salmon_exec,             $bamtools_exec, $sort_exec
  )
  = &check_program_optional(
                              'gs',       'inkscape', 'convert',
                             'pdfcrop', 'biokanga', 'samtools', 'bowtie2',
                             'bwa',     'salmon',  'bedtools', 'sort'
  );

if ($ps2pdf_exec) {
 $ps2pdf_exec .=
   " -sPAPERSIZE=a0 -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -q -sOutputFile=";
}

# TODO If there exist a sizeable number of housekeeping transcripts that should not be DE, the the dispersion could be estimated from them.
my $readset_separation     = 500;
my $fragment_max_length    = 800;
my $edgeR_dispersion       = 0.4;
my $minCPM                 = 2;
my $minLibs                = 1;
my $fdr_pval_cutoff        = 0.001;
my $tree_clusters          = 10;
my $read_format            = 'fastq';
my $threads                = 3;
my $min_housekeeping_genes = 5;
my $process_cutoff         = 1;
my $binary_min_coverage    = 0.3;
my $binary_min_reads       = 4;
my $sort_memory            = 10;
my $salmon_min_bias = 2*10^6; # 2 mil as we often have few reads
my (
     $uid,               $lib_alias_file,     $demand_all_readsets,
     $use_kanga,         @use_existing_bam,   $overwrite_results,
     @readsets2,         $db_file,            $no_checks,
     $db_use_file,       $dbname,             $no_graphs,
     $user_ref_sequence, $do_kangade,         $main_output_dir,
     $use_bwa,           $prepare_input_only, @sample_overlays,
     $initial_db,        $resume,             $gene_coverage_graphs_only,
     %binary_table,      $never_skip,         
     $extra_options,     $genomewide,         $only_alignments,
     $options_single,    $no_js_graphs,       $do_png,
     $no_pdf,            $no_pairwise,        $verbose, $is_quantseq
);
my $cwd = getcwd;
my $given_cmd = $0 . " " . join( " ", @ARGV );
my $sort_tmp = $cwd.'/sort_tmp';

pod2usage $! unless GetOptions(
 'sort_tmp:s'		=> \$sort_tmp,
 'do_galaxy_cleanup' => \$do_galaxy_cleanup, # delete EVERYTHING EXCEPT PDF AND TSV and sqlite
 'nocheck|no_check'      => \$no_checks,            #TMP for Debug; should not be needed
 'infile:s'      => \$input_reference_file,
 'extra_genes:s' => \$extra_genes,
 'sequence:s'    => \$user_ref_sequence,
 'format:s'      => \$read_format,
 '1read|readset1|r:s{1,}' => \@readsets,     # if bowtie2, it can be /bz2
 '2read|readset2:s{1,}'   => \@readsets2,    # if bowtie2, it can be /bz2

 'samtools_exec:s' => \$samtools_exec,
 'bowtie2_exec:s'  => \$bowtie2_exec,
 'bamtools_exec:s' => \$bamtools_exec,       # SE bowtie2 only
 'bwa_exec:s'      => \$bwa_exec,            # modified bwa
 'biokanga_exec:s' => \$biokanga_exec,

 'out|output|uid:s'          => \$uid,
 'outdir:s'                  => \$main_output_dir,
 'threads|cpus:i'            => \$threads,
 'alias|library_name_file:s' => \$lib_alias_file,
 'sample_names:s{,}'        => \@sample_names,
 'prepare_only'              => \$prepare_input_only,
 'need_all_readsets'         => \$demand_all_readsets,
 'over'                      => \$overwrite_results,
 'no_graphs|nographs'        => \$no_graphs,
 'hostname:s'                => \$db_hostname,
 'dbname:s'                  => \$dbname,
 'seeddb:s'                  => \$initial_db,
 'contextual'                => \$contextual_alignment,
 'correct_bias|isoforms'     => \$perform_bias_correction,
 'use_bwa'                   => \$use_bwa,
 'kanga'                     => \$use_kanga,
 'existing_aln:s{1,}'        => \@use_existing_bam,
 'coord_existing'        => \$existing_bam_is_coordsorted,
 'debug:i'               => \$debug,               # should not be used by users
 'verbose'               => \$verbose,
 'resume'                => \$resume,
 'do_kangade|dokangade'  => \$do_kangade,
 'db_use_file'           => \$db_use_file,
 'gene_coverage_graphs_only'      => \$gene_coverage_graphs_only,
 'dispersion:s'          => \$edgeR_dispersion,    #auto for bio.reps
 'fdr_cutoff:f'          => \$fdr_pval_cutoff,
 'cpm_cutoff'            => \$minCPM,
 'cutoff_library'        => \$minLibs,
 'binary_min_coverage:f' => \$binary_min_coverage,
 'binary_min_reads:i'    => \$binary_min_reads,
 'never_skip'            => \$never_skip,
 'memory_gb|sort_memory_gb:i' => \$sort_memory,
 'remove_redund'              => \$remove_redund,
 'readset_separation:i'       => \$readset_separation,
 'extra_options:s'            => \$extra_options,
 'genomewide'                 => \$genomewide,
 'only_alignments'            => \$only_alignments,
 'cleanup'                    => \$cleanup,
 'no_js_graphs'               => \$no_js_graphs,
 'png_graphs'                 => \$do_png,
 'nopairwise|no_pairwise'     => \$no_pairwise,
 'options_single:s'           => \$options_single,
 'nopdf|no_pdf'               => \$no_pdf,
 'salmon_min_bias:i'          => \$salmon_min_bias,
 'is_quantseq'		      => \$is_quantseq,
);

die
  "-binary_min_coverage has to be between 0 and 1 (not $binary_min_coverage)\n"
  unless $binary_min_coverage > 0 && $binary_min_coverage <= 1;
$threads = 2 if $threads < 2;    # no idea what happens with a single thread
my ($bunzip2_exec) = &check_program_optional('pbzip2');
my $bunzip_threads = $threads <= 6 ? $threads : 6;
$bunzip2_exec .= " -p$bunzip_threads " if $bunzip2_exec;
($bunzip2_exec) = &check_program('bzip2') if !$bunzip2_exec;

#MEMORY & threads for sorts
die "Please provide memory as an integer of gigabytes, e.g. -sort_memory_gb 10\n" unless $sort_memory && $sort_memory=~/^\d+$/ && $sort_memory >=1;
$sort_memory = $sort_memory * 1024 * 1024 * 1024;
my $sort_memory_exec = $sort_memory;
my $sort_memory_sub = int($sort_memory/3).'b';
my $sort_version = `sort --version|head -1`;
if ($sort_version=~/^sort \(GNU coreutils\) (\d+)/){
   my $sort_threads = $threads > 5 ? 5 : 5;
   # for pipelines
   my $sort_threads_sub = int($threads / 3);
   $sort_threads_sub = $sort_threads_sub > 5 ? 5 : 5;
   $sort_memory_exec .=  " --parallel=$sort_threads" if ($1 && $1 >= 8);
   $sort_memory_sub .= " --parallel=$sort_threads_sub" if ($1 && $1 >= 8);
}
mkdir($sort_tmp) if !-d $sort_tmp;

# parallelise alignments; each alignment uses 4 threads
my $alignment_threads = 4;
my $alignment_helper_threads = int($threads/$alignment_threads) > 0 ? int($threads/$alignment_threads) : 1;
my $alignment_thread_helper = new Thread_helper($alignment_helper_threads);

my $samtools_threads = $threads > 5 ? 5 : $threads;
my $sam_sort_memory = int($sort_memory / $samtools_threads); 
my ( $extra_salmon, %single_end_readsets, $salmon_options_single );
my $kanga_threads = $threads <= 8 ? $threads : 10;
my $R_threads = $threads <= 8 ? $threads : 10;


$uid = &get_uid_time('dew') unless $uid;
$dbname = $uid . '_transcriptome.sqlite.db' unless $dbname;
my $result_dir =
    $main_output_dir
  ? $cwd . "/" . $main_output_dir . '/'
  : $cwd . "/" . $uid . '_results/';
my $edgeR_dir = $result_dir . '/edgeR/';
my $counts_expression_level_matrix =
  "$result_dir/$uid.counts.expression_levels.matrix";
my $total_timer = new Time::Progress;
&perform_checks_preliminary();
my (
     $dbh,
     $get_from_seqdb,
     $get_hash_from_seqdb,
     $add_seqhash_to_seqdb,
     $init_expression_statistics,
     $add_readset,
     $update_readset,
     $update_readset_size,
     $add_depth_data,
     $check_depth_data,
     $check_depth_readset,
     $delete_depth_data,
     $get_readset,
     $get_expression_statistics,
     $check_expression_statistics,
     $update_expression_statistics,
     $set_housekeeping,
     $set_housekeeping_fold,
     $check_hash_from_seqdb,
     $add_sequence_alias,
     $get_sequence_alias,
     $add_to_kangade_analysis,
     $check_kangade_analysis,
     $delete_kangade_analysis,
     $start_fold_change,
     $check_fold_change,
     $add_kangade_fold_change,
     $add_salmon_fold_change,
     $add_raw_fold_change,
     $get_kangade_fold_change,
     $get_salmon_fold_change,
     $get_raw_fold_change,
     $update_salmon_expression_statistics,
     $update_kangade_expression_statistics,
     $update_rpkm_expression_statistics,
     $get_readset_filename
) = &sqlite_init();
my $file_for_alignment = &prepare_input_data();

# this stores all the expression data for creating the coverage graphs without database acceess
my %expression_coverage_stats_hash;

if ( !-s $counts_expression_level_matrix
     || ( $contextual_alignment && $debug && $debug >= 2 ) )
{
 &starts_alignments($file_for_alignment);
 if ($only_alignments) {
  print "User stop requested after alignments are completed.\n";
  exit(0);
 }
 &perform_stats();
 &process_expression_level();
 &sqlite_backup() unless $db_use_file;

 # from now on, we don't need to write to the database
 close STATS;
 close STATS_RATIO;
 &print_binary_table();
}
else {
 print "Expression data already exists ($counts_expression_level_matrix). Will not reprocess. Acquiring data from database...\n";
 # we need to get all the expression coverage data
 &get_all_expression_data();
 &print_binary_table();
}

if ($only_alignments) {
 print "User stop requested after alignments are completed.\n";
 &process_completion();
}

# get TMM normalized expression. this is relatively fast.
my $check_lines = `wc -l < $counts_expression_level_matrix`;
confess "No expression data available (empty $counts_expression_level_matrix)!\n"
  unless ( -s $counts_expression_level_matrix && $check_lines > 1 );
my ( $effective_expression_level_matrix_TMM_tpm )
  = &perform_TMM_normalization_edgeR($counts_expression_level_matrix);

$check_lines = `wc -l < $effective_expression_level_matrix_TMM_tpm `
  if $effective_expression_level_matrix_TMM_tpm;
confess
"edgeR did not produce any TMM output ($effective_expression_level_matrix_TMM_tpm)!\n"
  unless ( -s $effective_expression_level_matrix_TMM_tpm && $check_lines > 1 );

# this needs the database in order to do pairwise comparisons...
if ( !$no_pairwise ) {
 print
"\nPreparing edgeR differential expression data and graphs, This may take a long time if you have a large number of groups (use -no_pairwise to not do that in the future)\n";
 &perform_edgeR_pairwise();
 &prepare_edgeR_graphs($effective_expression_level_matrix_TMM_tpm);
 my ( $html2d, $html3d ) =
   &prepare_scatter_for_canvas($effective_expression_level_matrix_TMM_tpm);
}
else {
 print "User requested no pairwise edgeR comparisons...\n";
}

## from now on, we have no db access? we can multithread.
&sqlite_destroy();

my $expression_coverage_stats_hashref =
  shared_clone( \%expression_coverage_stats_hash );
undef(%expression_coverage_stats_hash);

&perform_coverage_graphs($expression_coverage_stats_hashref);
undef($expression_coverage_stats_hashref);

if ($gene_coverage_graphs_only) {
 print "User stop requested after gene coverage graphs were completed.\n";
 &process_completion();
}

unless ($no_graphs){
	&process_edgeR_graphs_overview( $effective_expression_level_matrix_TMM_tpm,
        	                        $result_dir . '/gene_expression_tpm/','TPM' );
}

&process_completion();

#########################################################################################################
sub process_cmd {
 my ( $cmd, $dir ) = @_;
 print "CMD: $cmd\n" if $debug || $verbose;
 print LOG "CMD: $cmd\n";
 chdir($dir) if $dir;
 my $ret = system($cmd);
 if ( $ret && $ret != 256 ) {
  confess "Error, cmd died with ret $ret\n";
 }
 chdir($cwd) if $dir;
 return;
}

sub mytime() {
 my @mabbr =
   qw(January February March April May June July August September October November December);
 my @wabbr = qw(Sunday Monday Tuesday Wednesday Thursday Friday Saturday);
 my $sec   = localtime->sec() < 10 ? '0' . localtime->sec() : localtime->sec();
 my $min   = localtime->min() < 10 ? '0' . localtime->min() : localtime->min();
 my $hour =
   localtime->hour() < 10 ? '0' . localtime->hour() : localtime->hour();
 my $wday = $wabbr[ localtime->wday ];
 my $mday = localtime->mday;
 my $mon  = $mabbr[ localtime->mon ];
 my $year = localtime->year() + 1900;
 return "$wday, $mon $mday, $year: $hour:$min:$sec\t";
}

sub sqlite_backup() {
 if ($debug) {
  warn "DEBUG mode: Will not backup database\n";
  return;
 }
 my $keep           = shift;
 my $backup_db_file = shift;
 $backup_db_file = $db_use_file ? $db_file . '.backup' : $db_file
   if !$backup_db_file;
 print "\nCheckpointing database ($backup_db_file). DO NOT halt...\r";
 $dbh->do("VACUUM");
 $dbh->do("PRAGMA shrink_memory");

 # ensure no problems if it crashes while backing up
 $dbh->sqlite_backup_to_file( $backup_db_file . ".tmp" );
 if ( -s $backup_db_file . ".tmp" ) {
  unlink( $backup_db_file . '.old' );
  rename( $backup_db_file, $backup_db_file . '.old' ) if $keep;
  unlink($backup_db_file);
  rename( $backup_db_file . ".tmp", $backup_db_file );
  print
"Checkpointing database ($backup_db_file). Done: SQL database checkpointed!\n";
 }
 else {
  print
"Checkpointing database ($backup_db_file). Error: Checkpointing failed or there were no data to write out!\n";
 }
}

sub sqlite_init() {
 $db_file = $debug ? $cwd . "/$dbname.debug" : $cwd . "/$dbname";
 my $active_db_file = $db_use_file ? $db_file . '.active' : ':memory:';
 if ( !-s $db_file && $initial_db && $initial_db ne $db_file ) {
  die "You asked for a seed database ($initial_db) but it doesn't exist!\n"
    unless -s $initial_db;
  warn "Will use $initial_db as seed database\n";
  copy( $initial_db, $db_file );
 }
 my $db_existed = -s $db_file ? 1 : int(0);
 my $dbh;
 $dbh = DBI->connect(
                      "dbi:SQLite:" . $active_db_file,
                      "", "",
                      {
                        sqlite_see_if_its_a_number => 1,
                        AutoCommit                 => 1,
                        RaiseError                 => 1
                      }
 );
 my $db_version = $dbh->{sqlite_version};
 confess "Cannot create SQLite database" unless $db_version;
 print "\tUsing SQLite DB v$db_version";
 if ($db_use_file) {
  print " (with a file)\n";
 }
 else {
  print " (in memory)\n";
 }
 if ($db_existed) {
  $dbh->sqlite_backup_from_file($db_file);
  my $res = `$sqlite3_exec $db_file "PRAGMA integrity_check"`;
  die "Database integrity check failed for $db_file. Delete and try again\n" if ($res && $res!~/^ok/);
 } 
 else {
  if ($resume || $no_checks){
    print
     "Warning! Database $db_file does not seem to be ok, will recreate and not resume.\n";
    undef($resume);
    undef($no_checks);
  }
  print "\tCreating database...\n";
  $dbh->do("PRAGMA encoding = 'UTF-8'");

#$dbh->do("CREATE TABLE file_alignments(file_md5sum char(32),readset_id integer,bam blob)"    );
#$dbh->do("CREATE UNIQUE INDEX file_alignments_idx1 ON file_alignments(file_md5sum,readset_id)"    );
  $dbh->do("CREATE TABLE sequence_data(seq_md5hash char(32) primary key,seq_length integer,housekeeping integer DEFAULT 0)"  );
  $dbh->do("CREATE TABLE sequence_aliases (seq_md5hash char(32), alias text)");
  $dbh->do("CREATE INDEX sequence_aliases_idx ON sequence_aliases(seq_md5hash)");
  $dbh->do("CREATE TABLE readsets (readset_id INTEGER PRIMARY KEY,readset_file varchar(255),total_reads integer, is_paired varchar, alias varchar(255), readlength_median integer, ctime timestamp)"  );
  $dbh->do("CREATE UNIQUE INDEX readsets_idx1 ON readsets(readset_file)");
  $dbh->do("CREATE TABLE expression_statistics (seq_md5hash char(32), readset_id integer, gene_length_coverage REAL, gene_length_coverage_mean REAL, no_coverage integer, rpkm integer, aligned_reads_per_base REAL, gene_length_coverage_median integer, total_aligned_reads integer, gene_length_coverage_max integer, gene_length_coverage_sd REAL, salmon_eff_counts REAL,  salmon_tpm REAL, kangade_counts integer)"  );
  $dbh->do("CREATE UNIQUE INDEX expression_statistics_idx1 ON expression_statistics(seq_md5hash,readset_id)"  );

  #tmp for uncached
  $dbh->do("CREATE TABLE expression_statistics_tmp (seq_md5hash char(32), readset_id integer, gene_length_coverage REAL, gene_length_coverage_mean REAL, no_coverage integer, rpkm integer, aligned_reads_per_base REAL, gene_length_coverage_median integer, total_aligned_reads integer, gene_length_coverage_max integer, gene_length_coverage_sd REAL, salmon_eff_counts REAL, salmon_tpm REAL , kangade_counts integer)"  );
  $dbh->do("CREATE UNIQUE INDEX expression_statistics_tmp_idx1 ON expression_statistics_tmp(seq_md5hash,readset_id)"  );

  # r-tree?
  $dbh->do("CREATE TABLE depth (seq_md5hash char(32), readset_id integer, data blob)");  $dbh->do("CREATE INDEX depth_idx1 ON depth(seq_md5hash,readset_id)");
  $dbh->do("CREATE TABLE depth_tmp (seq_md5hash char(32), readset_id integer, data blob)"  );
  $dbh->do("CREATE INDEX depth_tmp_idx1 ON depth_tmp(seq_md5hash,readset_id)");
  $dbh->do("CREATE TABLE kangade_analysis (seq_md5hash char(32), readset1_id INTEGER, readset2_id INTEGER,Classification INTEGER,Score INTEGER,DECntsScore INTEGER,PearsonScore INTEGER,CtrlUniqueLoci INTEGER,"
     . "ExprUniqueLoci INTEGER,CtrlExprLociRatio INTEGER,PValueMedian REAL,PValueLow95 REAL,PValueHi95 REAL,TotCtrlCnts INTEGER,TotExprCnts INTEGER,TotCtrlExprCnts INTEGER,ObsFoldChange REAL,FoldMedian REAL,"
     . "FoldLow95 REAL,FoldHi95 REAL,ObsPearson REAL,PearsonMedian REAL,PearsonLow95 REAL,PearsonHi95 REAL)"
  );
  $dbh->do("CREATE UNIQUE INDEX kangade_analysis_idx1 ON kangade_analysis(seq_md5hash,readset1_id,readset2_id)"  );
  $dbh->do("CREATE TABLE fold_changes (seq_md5hash char(32), readset1_id INTEGER, readset2_id INTEGER, raw_rpkm REAL, salmon_effective_counts REAL,salmon_tpm REAL , kangade_observed REAL,housekeeping integer DEFAULT 0)"  );
  $dbh->do("CREATE UNIQUE INDEX fold_changes_idx1 ON fold_changes(seq_md5hash,readset1_id,readset2_id)"  );
 }
 ### PRAGMAS for speed
 $dbh->do("PRAGMA journal_mode = MEMORY");
 $dbh->do("PRAGMA temp_store = 2 ");          # memory
 $dbh->do("PRAGMA cache_size = -448000 ");    # 400mb of RAM for sqlite
 $dbh->do("PRAGMA synchronous = OFF") if $db_use_file;
 $dbh->do("PRAGMA quick_check");
 ### SQL queries
 my $start_fold_change = $dbh->prepare(
"INSERT INTO fold_changes (seq_md5hash,readset1_id,readset2_id) VALUES (?,(SELECT readset_id from readsets where readset_file=?),(SELECT readset_id from readsets where readset_file=?))"
 );
 my $check_fold_change = $dbh->prepare(
"SELECT seq_md5hash,readset1_id,readset2_id FROM fold_changes WHERE seq_md5hash=? AND readset1_id=(SELECT readset_id from readsets WHERE readset_file=?) AND readset2_id=(SELECT readset_id from readsets WHERE readset_file=?)"
 );
 my $add_kangade_fold_change = $dbh->prepare(
"UPDATE fold_changes SET kangade_observed=? WHERE seq_md5hash=? AND readset1_id=(SELECT readset_id from readsets WHERE readset_file=?) AND readset2_id=(SELECT readset_id from readsets WHERE readset_file=?)"
 );
 my $add_salmon_fold_change = $dbh->prepare(
"UPDATE fold_changes SET salmon_effective_counts=?,salmon_tpm=? WHERE seq_md5hash=? AND readset1_id=(SELECT readset_id from readsets where readset_file=?) AND readset2_id=(SELECT readset_id from readsets where readset_file=?)"
 );
 my $add_raw_fold_change = $dbh->prepare(
"UPDATE fold_changes SET raw_rpkm=? WHERE seq_md5hash=? AND readset1_id=(SELECT readset_id from readsets where readset_file=?) AND readset2_id=(SELECT readset_id from readsets where readset_file=?)"
 );
 my $get_kangade_fold_change = $dbh->prepare(
"SELECT kangade_observed FROM fold_changes WHERE seq_md5hash=? AND readset1_id=(SELECT readset_id from readsets where readset_file=?) AND readset2_id=(SELECT readset_id from readsets where readset_file=?)"
 );
 my $get_salmon_fold_change = $dbh->prepare(
"SELECT salmon_tpm, salmon_effective_counts FROM fold_changes WHERE seq_md5hash=? AND readset1_id=(SELECT readset_id from readsets where readset_file=?) AND readset2_id=(SELECT readset_id from readsets where readset_file=?)"
 );
 my $get_raw_fold_change = $dbh->prepare(
"SELECT raw_rpkm FROM fold_changes WHERE seq_md5hash=? AND readset1_id=(SELECT readset_id from readsets where readset_file=?) AND readset2_id=(SELECT readset_id from readsets where readset_file=?)"
 );
 my $add_to_kangade_analysis = $dbh->prepare(
"INSERT INTO kangade_analysis (seq_md5hash, readset1_id, readset2_id,Classification,Score,DECntsScore,PearsonScore,CtrlUniqueLoci,ExprUniqueLoci,CtrlExprLociRatio,PValueMedian,PValueLow95,PValueHi95,TotCtrlCnts,TotExprCnts,TotCtrlExprCnts,ObsFoldChange,FoldMedian,FoldLow95,FoldHi95,ObsPearson,PearsonMedian,PearsonLow95,PearsonHi95) "
    . "VALUES (?,(SELECT readset_id from readsets where readset_file=?),(SELECT readset_id from readsets where readset_file=?),?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)"
 );
 my $check_kangade_analysis = $dbh->prepare(
"SELECT * FROM kangade_analysis WHERE seq_md5hash=? AND readset1_id=(SELECT readset_id from readsets where readset_file=?) AND readset2_id=(SELECT readset_id from readsets where readset_file=?) "
 );
 my $delete_kangade_analysis = $dbh->prepare(
"DELETE FROM kangade_analysis WHERE seq_md5hash=? AND readset1_id=(SELECT readset_id from readsets where readset_file=?) AND readset2_id=(SELECT readset_id from readsets where readset_file=?)"
 );

#my $readset_id_to_filename = $dbh->prepare_cached("SELECT readset_file FROM readsets WHERE readset_id = ? ");
#my $readset_filename_to_id = $dbh->prepare_cached("SELECT readset_id FROM readsets WHERE readset_file = ? ");
#my $add_to_db = $dbh->prepare("INSERT INTO file_alignments (file_md5sum,readset_id,bam) VALUES (?, (SELECT readset_id from readsets where readset_file=?), ?)"  );
#my $delete_from_db = $dbh->prepare("DELETE FROM file_alignments WHERE file_md5sum=? and readset_id=(SELECT readset_id from readsets where readset_file=?)"  );
#my $check_db = $dbh->prepare("SELECT readset_id FROM file_alignments WHERE file_md5sum=? AND readset_id=(SELECT readset_id from readsets where readset_file=?)"  );
#my $get_from_aligndb = $dbh->prepare("SELECT bam FROM file_alignments WHERE file_md5sum=? AND readset_id=(SELECT readset_id from readsets where readset_file=?)"  );
 my $get_from_seqdb = $dbh->prepare_cached(
                    "SELECT seq_length FROM sequence_data WHERE seq_md5hash=?");
 my $get_hash_from_seqdb =
   $dbh->prepare("SELECT seq_md5hash FROM sequence_aliases WHERE alias=?");
 my $check_hash_from_seqdb =
   $dbh->prepare("SELECT seq_md5hash FROM sequence_data WHERE seq_md5hash=?");
 my $add_seqhash_to_seqdb = $dbh->prepare(
             "INSERT INTO sequence_data (seq_md5hash,seq_length) VALUES (?,?)");
 my $depth_table = 'depth';
 $depth_table .= '_tmp' if $contextual_alignment && !$genomewide;
 my $check_depth_data = $dbh->prepare(
"SELECT data from $depth_table WHERE seq_md5hash=? AND readset_id=(SELECT readset_id from readsets where readset_file=?)"
 );
 my $check_depth_readset = $dbh->prepare(
"SELECT count(*) from $depth_table WHERE readset_id=(SELECT readset_id from readsets where readset_file=?)"
 );
 my $delete_depth_data = $dbh->prepare(
"DELETE from $depth_table WHERE seq_md5hash=? AND readset_id=(SELECT readset_id from readsets where readset_file=?) "
 );
 my $add_depth_data = $dbh->prepare(
"INSERT INTO $depth_table (seq_md5hash,readset_id,data) VALUES (?,(SELECT readset_id from readsets where readset_file=?),?)"
 );
 my $expression_statistics_table = 'expression_statistics';
 $expression_statistics_table .= '_tmp'
   if $contextual_alignment && !$genomewide;
 my $init_expression_statistics = $dbh->prepare(
"INSERT INTO $expression_statistics_table (seq_md5hash, readset_id) VALUES (?,(SELECT readset_id from readsets where readset_file=?)) "
 );
 my $update_rpkm_expression_statistics = $dbh->prepare(
"UPDATE $expression_statistics_table set rpkm=?, total_aligned_reads=? WHERE seq_md5hash=? AND readset_id=(SELECT readset_id from readsets where readset_file=?)"
 );
 my $update_expression_statistics = $dbh->prepare(
"UPDATE $expression_statistics_table set gene_length_coverage=?, gene_length_coverage_mean=?, no_coverage=?, aligned_reads_per_base=?, gene_length_coverage_median=?, gene_length_coverage_max=?, gene_length_coverage_sd=? WHERE seq_md5hash=? AND readset_id=(SELECT readset_id from readsets where readset_file=?)"
 );
 my $update_salmon_expression_statistics = $dbh->prepare(
"UPDATE $expression_statistics_table set salmon_eff_counts=?,salmon_tpm=? WHERE seq_md5hash=? AND readset_id=(SELECT readset_id from readsets where readset_file=?)"
 );
 my $update_kangade_expression_statistics = $dbh->prepare(
"UPDATE $expression_statistics_table set kangade_counts=? WHERE seq_md5hash=? AND readset_id=(SELECT readset_id from readsets where readset_file=?)"
 );
 my $get_expression_statistics = $dbh->prepare(
"SELECT * FROM $expression_statistics_table WHERE seq_md5hash=? AND readset_id=(SELECT readset_id from readsets where readset_file=?)"
 );
 my $check_expression_statistics = $dbh->prepare(
"SELECT seq_md5hash FROM $expression_statistics_table WHERE seq_md5hash=? AND readset_id=(SELECT readset_id from readsets where readset_file=?)"
 );
 my $set_housekeeping =
   $dbh->prepare("UPDATE sequence_data set housekeeping=? WHERE seq_md5hash=?");
 my $set_housekeeping_fold = $dbh->prepare(
"UPDATE fold_changes set housekeeping=? WHERE seq_md5hash=? AND readset1_id=? AND readset2_id=? "
 );
 my $get_readset = $dbh->prepare("SELECT * from readsets WHERE readset_file=?");
 my $get_readset_filename =
   $dbh->prepare("SELECT readset_file from readsets WHERE alias=?");
 my $add_readset = $dbh->prepare(
"INSERT INTO readsets (readset_file,is_paired,alias,total_reads,readlength_median,ctime) VALUES (?,?,?,?,?,?)"
 );
 my $update_readset = $dbh->prepare( "UPDATE readsets SET is_paired=?, alias=? WHERE readset_file=?");
 my $update_readset_size = $dbh->prepare( "UPDATE readsets SET total_reads=? WHERE readset_file=?");
 my $add_sequence_alias = $dbh->prepare(
               "INSERT INTO sequence_aliases (seq_md5hash,alias) VALUES (?,?)");
 my $get_sequence_alias =
   $dbh->prepare("SELECT alias from sequence_aliases WHERE seq_md5hash=?");
 return (
          $dbh,
          $get_from_seqdb,
          $get_hash_from_seqdb,
          $add_seqhash_to_seqdb,
          $init_expression_statistics,
          $add_readset,
          $update_readset,
          $update_readset_size,
          $add_depth_data,
          $check_depth_data,
          $check_depth_readset,
          $delete_depth_data,
          $get_readset,
          $get_expression_statistics,
          $check_expression_statistics,
          $update_expression_statistics,
          $set_housekeeping,
          $set_housekeeping_fold,
          $check_hash_from_seqdb,
          $add_sequence_alias,
          $get_sequence_alias,
          $add_to_kangade_analysis,
          $check_kangade_analysis,
          $delete_kangade_analysis,
          $start_fold_change,
          $check_fold_change,
          $add_kangade_fold_change,
          $add_salmon_fold_change,
          $add_raw_fold_change,
          $get_kangade_fold_change,
          $get_salmon_fold_change,
          $get_raw_fold_change,
          $update_salmon_expression_statistics,
          $update_kangade_expression_statistics,
          $update_rpkm_expression_statistics,
          $get_readset_filename
 );
}

sub sqlite_set_get_as_housekeeping() {
 my $seq_md5hash = shift;
 my $set         = shift;
 my $readset1    = shift;
 my $readset2    = shift;
 if ( $set && $set == 2 ) {
  $set_housekeeping->execute( 1, $seq_md5hash );
 }
 if ( $readset1 && $readset2 ) {
  $set_housekeeping_fold->execute( 1, $seq_md5hash, $readset1, $readset2 );
 }
}

sub sqlite_destroy() {
 return unless $dbh;
 print &mytime() . "Closing SQL connections and backing-up\n";
 $get_readset_filename->finish();
 $get_from_seqdb->finish();
 $get_hash_from_seqdb->finish();
 $add_seqhash_to_seqdb->finish();
 $check_hash_from_seqdb->finish();
 $check_depth_data->finish();
 $check_depth_readset->finish();
 $delete_depth_data->finish();
 $add_depth_data->finish();
 $set_housekeeping->finish();
 $set_housekeeping_fold->finish();
 $init_expression_statistics->finish();
 $get_expression_statistics->finish();
 $check_expression_statistics->finish();
 $get_readset->finish();
 $update_expression_statistics->finish();
 $update_salmon_expression_statistics->finish();
 $update_kangade_expression_statistics->finish();
 $add_readset->finish();
 $update_readset->finish();
 $update_readset_size->finish();
 $get_sequence_alias->finish();
 $add_sequence_alias->finish();
 $add_to_kangade_analysis->finish();
 $check_kangade_analysis->finish();
 $delete_kangade_analysis->finish();
 $start_fold_change->finish();
 $check_fold_change->finish();
 $add_kangade_fold_change->finish();
 $add_salmon_fold_change->finish();
 $add_raw_fold_change->finish();
 $get_kangade_fold_change->finish();
 $get_salmon_fold_change->finish();
 $update_rpkm_expression_statistics->finish();
 $get_raw_fold_change->finish();

 if ( $contextual_alignment && !$genomewide && !$debug ) {
  print "Contextual non-genomewide alignment was requested, I will delete the tmp directories from the database\n";
  #empty temporary tables and re-create their schema
  $dbh->do("DROP TABLE expression_statistics_tmp");
  $dbh->do(
"CREATE TABLE expression_statistics_tmp (seq_md5hash char(32), readset_id integer, gene_length_coverage REAL, gene_length_coverage_mean REAL, no_coverage integer, rpkm integer, aligned_reads_per_base REAL, gene_length_coverage_median integer, total_aligned_reads integer, gene_length_coverage_max integer, gene_length_coverage_sd REAL, salmon_eff_counts REAL, salmon_tpm REAL , kangade_counts integer)"
  );
  $dbh->do(
"CREATE UNIQUE INDEX expression_statistics_tmp_idx1 ON expression_statistics_tmp(seq_md5hash,readset_id)"
  );
  $dbh->do("DROP TABLE depth_tmp");
  $dbh->do(
"CREATE TABLE depth_tmp (seq_md5hash char(32), readset_id integer, data blob)"
  );
  $dbh->do("CREATE INDEX depth_tmp_idx1 ON depth_tmp(seq_md5hash,readset_id)");
  $dbh->do("VACUUM");
 }

 #&sqlite_backup() unless $db_use_file;
 my $backup_db_file = $db_use_file ? $db_file . '.backup' : $db_file;
 print "\nCheckpointing database ($backup_db_file). DO NOT halt...\r";

 # ensure no problems if it crashes while backing up
 $dbh->sqlite_backup_to_file( $backup_db_file . ".tmp" );
 if ( -s $backup_db_file . ".tmp" ) {
  unlink( $backup_db_file . '.old' );
  rename( $backup_db_file, $backup_db_file . '.old' );
  unlink($backup_db_file);
  rename( $backup_db_file . ".tmp", $backup_db_file );
  print
"Checkpointing database ($backup_db_file). Done: SQL database checkpointed!\n";
 }
 $dbh->disconnect();
 undef($dbh);
 unlink( $db_file . '.active' ) if $db_use_file;
}

sub not_used1() {
####
## how to store the entire alignment in the database.... not used.
 #sub sqlite_check_align_old($) {
 #  my $readset = shift;
 #  $check_db->execute( $reference_file_md5sum, $readset );
 #  my $result = $check_db->fetchrow_arrayref();
 #  $result = $result ? int(1) : int(0);
 #  print "Readset $readset was already in DB. Not re-aligning.\n" if $result;
 #  return $result;
 #}
 #sub sqlite_add_align_old($) {
 #  my $readset    = shift;
 #  my $bamfile    = shift;
 #  my $bamcontent = `cat $bamfile`;
 #  $add_to_db->bind_param( 1, $reference_file_md5sum );
 #  $add_to_db->bind_param( 2, $readset );
 #  $add_to_db->bind_param( 3, $bamcontent, SQL_BLOB );
 #  $add_to_db->execute();
 #  $check_db->execute( $reference_file_md5sum, $readset );
 #  my $result = $check_db->fetchrow_arrayref();
 #  $result = $result ? int(1) : int(0);
 #  print "Readset $readset added to DB\n" if $result;
 #  return $result;
 #}
 #sub sqlite_get_align_old($) {
 #  my $readset = shift;
 #  my $bam     = shift;
 #  $get_from_aligndb->execute( $reference_file_md5sum, $readset );
 #  my $row = $get_from_aligndb->fetchrow_arrayref();
 #  my $result = $row ? int(1) : int(0);
 #  open( BAM, ">$bam" );
 #  print BAM $row->[0];
 #  close BAM;
 #  &process_cmd("$samtools_exec index $bam")
 #    unless -s $bam . '.bai' && ( -s $bam . '.bai' ) > 200;
 #  return $result;
 #}
}

sub sqlite_get_seq_length($) {
 my $md5sum = shift;
 $get_from_seqdb->execute($md5sum);
 my $row = $get_from_seqdb->fetchrow_arrayref();
 return $row->[0];
}

sub sqlite_check_kangade_analysis() {
 my ( $md5_sum, $read1, $read2 ) = @_;
 $check_kangade_analysis->execute( $md5_sum, $read1, $read2 );
 my $res = $check_kangade_analysis->fetchrow_hashref();
 return $res;
}

sub sqlite_get_md5($) {
 my $id = shift;
 $get_hash_from_seqdb->execute($id);
 my $row = $get_hash_from_seqdb->fetchrow_arrayref();
 unless ($row) {
  confess "Cannot find md5 data for $id\n";
 }
 return $row->[0];
}

sub sqlite_add_seq_md5($$$) {
 my $id         = shift;
 my $md5sum     = shift;
 my $seq_length = shift;
 $check_hash_from_seqdb->execute($md5sum);
 my $check = $check_hash_from_seqdb->fetchrow_arrayref();
 $add_seqhash_to_seqdb->execute( $md5sum, $seq_length ) if !$check;
 $get_sequence_alias->execute($md5sum);
 my %existing_ids;

 while ( my $result = $get_sequence_alias->fetchrow_arrayref() ) {
  $existing_ids{ $result->[0] } = 1;
 }
 $add_sequence_alias->execute( $md5sum, $id ) unless $existing_ids{$id};
 $gene_aliases_threads{$md5sum}{$id}=1;
}

sub sqlite_get_seq_aliases() {
 my $md5sum = shift;
 my %existing_ids;
 $get_sequence_alias->execute($md5sum);
 while ( my $result = $get_sequence_alias->fetchrow_arrayref() ) {
  $existing_ids{ $result->[0] } = 1;
 }
 return keys %existing_ids;
}

sub sqlite_get_readset_metadata($) {
 my $readset_filename = shift;
 $get_readset->execute($readset_filename);
 my $result = $get_readset->fetchrow_hashref();
 return $result;
}

sub sqlite_update_readset_size(){
  my $readset_filename  = shift;
  my $new_size = shift;
  $update_readset_size->execute( $new_size ,$readset_filename);
}

sub sqlite_add_readset_metadata($$$) {
 my $readset_filename  = shift;
 my $lib_alias         = shift;
 my $total_reads       = shift;
 my $readlength_median = shift;
 my $is_paired =
     $paired_readset_lookup{$readset_filename}
   ? $paired_readset_lookup{$readset_filename}
   : '0';
 my $result = &sqlite_get_readset_metadata($readset_filename);
 if ( !$result ) {
  $add_readset->execute( $readset_filename, $is_paired,         $lib_alias,
                         $total_reads,      $readlength_median, localtime() );
 }
 elsif ( !$result->{'alias'}
         || ( $result->{'alias'} && $result->{'alias'} ne $lib_alias ) )
 {
  $update_readset->execute( $is_paired, $lib_alias, $readset_filename );
 }

}

sub sqlite_check_expression_statistics($$) {
 my ( $seq_md5hash, $readset ) = @_;
 $check_expression_statistics->execute( $seq_md5hash, $readset );
 my $result = $check_expression_statistics->fetchrow_arrayref();
 return $result && $result->[0] ? 1 : undef;
}

sub sqlite_init_expression_statistics($) {
 my ( $seq_md5hash, $readset ) = @_;
 confess "No readset ID for initializing sqlite\n" unless $readset;
 unless (
      defined( &sqlite_check_expression_statistics( $seq_md5hash, $readset ) ) )
 {
  $init_expression_statistics->execute( $seq_md5hash, $readset );
 }
}

sub sqlite_add_rpkm_expression_statistics($) {
 my ( $seq_md5hash, $readset, $rpkm, $total_reads_hit ) = @_;
 $update_rpkm_expression_statistics->execute( $rpkm, $total_reads_hit,
                                              $seq_md5hash, $readset );
 my $check = &sqlite_get_expression_statistics( $seq_md5hash, $readset );
 if ( !$check || !defined( $check->{'rpkm'} ) ) {
  die
"Could not add RPKM expression statistics for: $seq_md5hash,$readset,$rpkm,$total_reads_hit\n";
 }
}

sub sqlite_add_main_expression_statistics($) {
 my (
      $seq_md5hash, $readset,     $gene_length_coverage_mean, $no_coverage,
      $aligned_reads_per_base,  $gene_length_coverage_median, $gene_length_coverage_max,  $gene_length_coverage_sd, $gene_length_coverage
 ) = @_;
 $update_expression_statistics->execute(
                           $gene_length_coverage, $gene_length_coverage_mean, $no_coverage, $aligned_reads_per_base,  $gene_length_coverage_median,
                           $gene_length_coverage_max,  $gene_length_coverage_sd,      $seq_md5hash, $readset );
}

sub sqlite_add_salmon_expression_statistics() {
 my ( $seq_md5hash, $readset, $eff_counts, $tpm ) = @_;
 my $r =
   $update_salmon_expression_statistics->execute( $eff_counts, $tpm,
                                                   $seq_md5hash, $readset );
 my $check = &sqlite_get_expression_statistics( $seq_md5hash, $readset );
 if ( !$check || !defined( $check->{'salmon_tpm'} ) ) {

   # warn Dumper $check if $debug;
  die "Could not add salmon expression statistics for: $seq_md5hash,$readset,$eff_counts,$tpm\n";
 }
}

sub sqlite_add_kangade_expression_statistics() {
 my ( $seq_md5hash, $readset, $kangade_counts ) = @_;
 $update_kangade_expression_statistics->execute( int($kangade_counts),
                                                 $seq_md5hash, $readset );
 my $check = &sqlite_get_expression_statistics( $seq_md5hash, $readset );
 if ( !$check || !defined( $check->{'kangade_counts'} ) ) {

  #warn Dumper $check if $debug;
  die
"Could not add kangade expression statistics for: $seq_md5hash,$readset,$kangade_counts\n";
 }
}

sub sqlite_get_expression_statistics($$) {

 my ( $seq_md5hash, $readset ) = @_;
 $get_expression_statistics->execute( $seq_md5hash, $readset );
 return $get_expression_statistics->fetchrow_hashref();
}

sub sqlite_check_depth_readset($) {

 # this checks if a readset has been processed
 my $readset_filename = shift;
 $check_depth_readset->execute($readset_filename);
 my $ref = $check_depth_readset->fetchrow_arrayref();
 return $ref->[0] if $ref;
}

sub sqlite_add_depth_data($$$) {
 my $seq_md5hash                = shift;
 my $readset_filename           = shift;
 my $hash_ref                   = shift;
 my $freeze                     = freeze($hash_ref);
 my $base_depth_data_serialized = compress( \$freeze )
   || confess "compression failed\n";

 $check_depth_data->execute( $seq_md5hash, $readset_filename );
 my $result = $check_depth_data->fetchrow_arrayref();

 $delete_depth_data->execute( $seq_md5hash, $readset_filename ) if $result;
 $add_depth_data->bind_param( 1, $seq_md5hash );
 $add_depth_data->bind_param( 2, $readset_filename );
 $add_depth_data->bind_param( 3, $base_depth_data_serialized, SQL_BLOB );
 $add_depth_data->execute();
 $check_depth_data->execute( $seq_md5hash, $readset_filename );
 undef($result);
 $result = $check_depth_data->fetchrow_arrayref();
 confess "Could not add depth data for $seq_md5hash, $readset_filename " unless $result;

 #  warn "Adding depth data for $seq_md5hash vs $readset_filename\n" if $debug;
}

sub sqlite_get_depth_data($$) {
 my $seq_md5hash         = shift;
 my $readset_filename    = shift;
 my $return_uncompressed = shift;
 $check_depth_data->execute( $seq_md5hash, $readset_filename );
 my $result = $check_depth_data->fetchrow_arrayref();

 if ( $result->[0] ) {
  return $result->[0] if $return_uncompressed;
  my $ref = decompress( \$result->[0] ) || confess "Decompression failed\n";
  return thaw($ref);
 }
 return;
}

sub sqlite_start_fold_change() {
 my ( $seq_md5hash, $readset1, $readset2 ) = @_;
 my $check = $check_fold_change->execute( $seq_md5hash, $readset1, $readset2 );
 if ( !$check ) {
  $start_fold_change->execute( $seq_md5hash, $readset1, $readset2 );
 }
 $check = $check_fold_change->execute( $seq_md5hash, $readset1, $readset2 );
 confess "Failed to populate fold change database\n" unless $check;
}

sub perform_checks_preliminary() {
 print &mytime() . "Performing preliminary checks\n";

 pod2usage "SQLite 3 not found\n" unless $sqlite3_exec;
 pod2usage "BWA not found\n"
   if $use_bwa && ( !$bwa_exec || !-s $bwa_exec || !-x $bwa_exec );
 pod2usage "Samtools not found\n"
   unless $samtools_exec && -s $samtools_exec && -x $samtools_exec;

 my $samtools_version = `$samtools_exec 2>&1 |head -n 3|tail -n 1`;
# die "Because the Samtools folks keep changing the option lists, now Samtools 0.1.19 is required. See 3rd_party/ and make sure you ran make\n" unless $samtools_version && $samtools_version=~/0\.1\.19/;

 $read_format = lc($read_format);
 pod2usage "Read format can only be BAM or FASTQ\n"
   unless $read_format eq 'bam' || $read_format eq 'fastq';
 pod2usage "No input\n"
   unless (    ( $input_reference_file && -s $input_reference_file )
            || ( $user_ref_sequence && length($user_ref_sequence) > 100 ) );
 pod2usage "Insufficient readsets (need at least 2!)\n"
   unless @readsets && scalar(@readsets) > 1;

 # biokanga:
 if ( $biokanga_exec && -x $biokanga_exec ) {
  $kangade_exec = $biokanga_exec . ' rnade';
  $kangax_exec  = $biokanga_exec . ' index';
  $kanga_exec   = $biokanga_exec . ' align';
 }
 else {
  undef($do_kangade);
 }

 if ( scalar(@readsets) < 2 ) {
  print
"Only one readset provided. Differential expression will not be run and will exit after the gene graphs are drawn.\n";
  $no_pairwise = 1;
 }

 undef($do_kangade) if ($no_pairwise);

 $never_skip = 1 if $genomewide;

 if (@readsets2) {
  if ( scalar(@readsets2) < scalar(@readsets) ) {
   warn
"Number of files for second readset does not equal to number of 1st readset. Extra files will be treated as single end readsets.\n";
  }
  elsif ( scalar(@readsets2) > scalar(@readsets) ) {
   die
"Number of files for second readset is larger than number for 1st readset. Single end datasets have to be provided with -1read not -2read. Stopping...\n";
  }
 }
 print "Using " . scalar(@readsets) . " readsets\n";
 for ( my $i = 0 ; $i < @readsets ; $i++ ) {
   next unless $readsets[$i];
  if ( !-s $readsets[$i] && -s $readsets[$i] . '.bam' ) {
   $readsets[$i] .= '.bam';
  }
  elsif ( !-s $readsets[$i] && -s $readsets[$i] . '.bz2' ) {
   $readsets[$i] .= '.bz2';
  }
  pod2usage "File " . $readsets[$i] . " not found\n"
    unless $readsets[$i] && -s $readsets[$i];
  $paired_readset_lookup{ $readsets[$i] } = 1;
  if ( @use_existing_bam && !$readsets[$i] ) {
   confess "Too many bam files provided. No readset number "
     . ( $i + 1 )
     . " provided\n";
  }
  elsif (@use_existing_bam) {
   $use_existing_bam[$i] = abs_path($use_existing_bam[$i]) if $use_existing_bam[$i];
   confess "Cannot find user-provided BAM file ("
     . $use_existing_bam[$i]
     . ") number "
     . ( $i + 1 ) . " for "
     . $readsets[$i] . "\n"
     unless $use_existing_bam[$i] && -s $use_existing_bam[$i];
     #convert to full path
  }
 }

 if (@readsets2) {
  for ( my $i = 0 ; $i < @readsets2 ; $i++ ) {
   if ( !-s $readsets2[$i] && -s $readsets2[$i] . '.bam' ) {
    $readsets2[$i] .= '.bam';
    $read_format = '.bam';
   }
   elsif ( !-s $readsets2[$i] && -s $readsets2[$i] . '.bz2' ) {
    $readsets2[$i] .= '.bz2';
   }
   if ( $readsets2[$i] && -s $readsets2[$i] ) {
    confess "Sorry, paired end Bowtie does not work with BAM files\n"
      if ( $read_format eq 'bam' && !$use_bwa );
    $paired_readset_lookup{ $readsets[$i] } = $readsets2[$i];
   }
  }
 }

 if ($extra_genes) {
  pod2usage "Cannot find $extra_genes file\n" unless -s $extra_genes;
 }
 if ( -d $result_dir && !$overwrite_results ) {
  pod2usage
"Result dir $result_dir already exists. If you wish to overwrite give the option -over otherwise delete the directory\n";
 }
 unless ($debug) {
  print "Checking for R installations\n";
  system(
'R --slave --no-restore --no-save -e \'source("$RealBin/R/dew_funcs.R");dew_install();\'  2> R.err'
  );
  my @error_check = `grep -i 'there is no package called' R.err`;
  if (@error_check) {
   undef(@error_check);
   system(
'R --slave --no-restore --no-save -e \'source("$RealBin/R/dew_funcs.R");dew_install();\'  2> R.err'
   );
   @error_check = `grep -i 'there is no package called' R.err`;
   confess "Some R packages are not installed:\n" . join( '', @error_check )
     if @error_check;
  }
  unlink("R.err");
 }
 mkdir($result_dir) unless -d $result_dir;
 open( LOG, ">>$result_dir/$uid.log" ) || die($!);
 print LOG "#Command: " . $given_cmd . "\n";
 print LOG "#Started: " . &mytime() . "\n";

 #"express:r-stranded;express:max-read-len 250"

 my $is_lib_provided; # to add ' --libType a ';

 if ($extra_options) {
  my @opt1 = split( ';', $extra_options );
  foreach my $op1 (@opt1) {
   my @opt2 = split( ':', $op1 );
   next unless $opt2[1];
   if ($opt2[1] =~/libType/){$is_lib_provided++;}
   $extra_salmon .= ' --' . $opt2[1] if $opt2[0] && $opt2[0] eq 'salmon';
  }
  print "All data: Applying extra salmon options: $extra_salmon\n";
  print LOG "All data: Applying extra salmon options: $extra_salmon\n";
 }
 if ($options_single) {
  my @opt1 = split( ';', $options_single );
  foreach my $op1 (@opt1) {
   my @opt2 = split( ':', $op1 );
   if ($opt2[1] =~/libType/){$is_lib_provided++;}
   $salmon_options_single .= ' --' . $opt2[1] if $opt2[0] eq 'salmon';
  }
  print "SE data: Applying extra salmon options: $salmon_options_single\n";
  print LOG
    "SE data: Applying extra salmon options: $salmon_options_single\n";
 }

 $extra_salmon .= ' --libType a ' if (!$is_lib_provided);

 mkdir($edgeR_dir) unless -d $edgeR_dir;
 mkdir( $edgeR_dir . '/js_plots' ) unless -d $edgeR_dir . '/js_plots';
 mkdir( $result_dir . "gene_coverage_plots" )
   unless -d $result_dir . "gene_coverage_plots";

 # maybe one day:
 #  if (@housekeeping_ids) {
 #    print LOG "Housekeeping genes requested:\n";
 #    if ( -s $housekeeping_ids[0] ) {
 #      foreach my $f (@housekeeping_ids) {
 #        open( H, $f ) || die($!);
 #        while ( my $ln = <H> ) {
 #          next if $ln =~ /^\s*$/;
 #          chomp($ln);
 #          $housekeeping{'ids'}{$ln} = 1;
 #          print LOG "$ln\n";
 #        }
 #        close(H);
 #      }
 #    } else {
 #      foreach my $id (@housekeeping_ids) {
 #        $housekeeping{'ids'}{$id} = 1 if $id;
 #        print LOG "$id\n";
 #      }
 #    }
 #    print LOG "\n";
 #  }
}

sub prepare_library_alias() {

 # these are currently not stored in the database
 my ( %library_metadata_headers, $print );
 if ( $lib_alias_file && -s $lib_alias_file ) {
  print "Parsing alias file $lib_alias_file...\n";
  open( IN, $lib_alias_file ) || die($!);
  my $header = <IN>;

  # 'group' is used for differential expression with R, otherwise no replicates
  my ($group_exists);
  chomp($header);
  my @headers = split( "\t", $header );
  my $original_header_number = scalar(@headers);
  die "Library alias ($lib_alias_file) must have 'file' and 'name' as the first two columns."
    . " The other columns are free to be any kind of metadata."
    . " Also the special metadata column 'group' is used for differential expression."
    unless $headers[0] eq 'file' && $headers[1] eq 'name';
  foreach my $h (@headers) {
   $group_exists++ if $h eq 'group';
  }

  # if no group, add it to the end.
  if ( !$group_exists ) {
   push( @headers, 'group' );
  }
  $print = join( "\t", @headers ) . "\n";
  while ( my $ln = <IN> ) {
   next if $ln =~ /^\s*$/;
   chomp($ln);
   my @data = split( "\t", $ln );
   next unless $data[0];
   die "Peculiar lib_alias file: number of headers ($original_header_number) does not equal no. of data (".scalar(@data).") at line:\n$ln\n"
    unless scalar(@data) == $original_header_number;
   
   if ( !-s $data[0] && $paired_readset_lookup{ $data[0] } ) {
    confess "Cannot find file " . $data[0] . "\n";
   }
   elsif ( !$paired_readset_lookup{ $data[0] } ) {
    warn "Lib alias entry " . $data[0] . " not in readset request. Skipping\n";
    next;
   }
   elsif ( !$data[1] ) {
    warn "Lib alias entry "
      . $data[0]
      . " has no name entry. Will use basename of filename\n";
    $data[1] = fileparse( $data[0] );

   }

   $data[1] =~ s/\W+/_/g;

   confess "Library name ("
     . $data[1]
     . ") cannot start with a digit (an R issue)\n"
     if $data[1] =~ /^\d+/;
   confess "Library name ("
     . $data[1]
     . ") cannot have the ^ hat/carrot symbol\n"
     if $data[1] =~ /^\d+/;

   $library_aliases{ $data[0] } = $data[1];
   if ( $data[2] ) {
    for ( my $i = 2 ; $i < @data ; $i++ ) {
     $library_metadata{ $data[1] }{ $headers[$i] } = $data[$i];
     $library_metadata_headers{ $headers[$i] }++;
    }
   }

   if ( !$library_metadata{ $data[1] }{'group'} ) {
    $library_metadata_headers{'group'}++;
    $library_metadata{ $data[1] }{'group'} = $data[1];
    $print .= join( "\t", @data ) . "\t" . $data[1] . "\n";
   }
   else {
    $print .= join( "\t", @data ) . "\n";
   }

   confess "Readset "
     . $data[1]
     . " has already been linked with group "
     . $library_metadata{ $data[1] }{'group'} . "\n"
     if $groups_readsets{ $library_metadata{ $data[1] }{'group'} }{ $data[1] };

   # for edgeR counts
   $groups_readsets{ $library_metadata{ $data[1] }{'group'} }{ $data[1] } =
     $edgeR_dir . $data[1] . '.dat';
  }
  close IN;
 }
 else {
  print "Creating a new LIBRARY alias file\n";
  # there is no lib_alias. make one
  $print = "file\tname\tgroup\n";
 }

 #some checks if there is a readset without data. specifically the group
  for (my $r=0;$r<scalar(@readsets);$r++){
   my $readset=$readsets[$r];

   if ($library_aliases{ $readset }){
	$sample_names[$r] = $library_aliases{ $readset };
   }elsif(!$sample_names[$r]){
	$sample_names[$r] = fileparse($readset);
   }
   # don't add if already a group
   next if !$readset; next if ($library_aliases{$readset} && $library_metadata{$sample_names[$r]}{'group'});

   $library_aliases{$readset} = $sample_names[$r];
   $library_metadata{$sample_names[$r]}{'group'} = $sample_names[$r];
   $groups_readsets{$library_metadata{$sample_names[$r]}{'group'}}{$sample_names[$r]} =$edgeR_dir . $sample_names[$r] . '.dat';;

   $print .=
       $readset . "\t"
     . $library_aliases{$readset} . "\t"
     . $library_metadata{$sample_names[$r]}{'group'} . "\n";
 }
 if ($print) {
  open( OUT, ">$result_dir/lib_alias.txt" );
  print OUT $print;
  close OUT;
 }
 @sample_overlays = sort keys %library_metadata_headers;
}

sub prepare_input_data() {
 my $do_backup = 0;

 &prepare_library_alias();

 print &mytime() . "Preparing readset metadata tables...\n";
 for ( my $i = 0 ; $i < @readsets ; $i++ ) {
  if ( $readsets2[$i] ) {
   my ( $readset_name, $library_size, $do_backup ) =
     &perform_readset_metadata( $readsets[$i], $readsets2[$i] );
  }
  else {
   my ( $readset_name, $library_size, $do_backup ) =
     &perform_readset_metadata( $readsets[$i] );
  }
 }
 &sqlite_backup( 1, $db_file . '.counted_readsets.seed' ) unless $initial_db || -s $db_file . '.counted_readsets.seed';

 &sqlite_backup() if !$db_use_file && $do_backup;
 print "\n" . &mytime() . "Preparing references...\n";

 # prepare alignment file.
 my $file_to_align = $result_dir . "$uid.toalign";

 # Once: create input file by validating user input
 if ( !-s $file_to_align ) {
  open( FASTA, ">$file_to_align" ) || die($!);

  # create input file if user provided sequence on cmd line
  my $reference_file = abs_path($input_reference_file);
  if ($user_ref_sequence) {
   $user_ref_sequence = ~s/^(>.+\\n)//;
   my $user_id = $1 ? $1 : ">query\n";
   $user_ref_sequence =~ s/\\n/\n/g;
   $user_ref_sequence = &seq_cleanup($user_ref_sequence);
   confess "Sequence from user contains non-ATCGN characters!\n"
     if $user_ref_sequence =~ /[^ATCGN]/;
   print FASTA $user_id . $user_ref_sequence . "\n";
  }
  else {

# prepare a cleaned up fasta, optionally adding extra genes for housekeeping checks etc
   my $reference_file_obj = new Fasta_reader($reference_file);
   while ( my $seq_obj = $reference_file_obj->next() ) {
    my $id = $seq_obj->get_accession() ? $seq_obj->get_accession() : 'query';
    my $sequence = $seq_obj->get_sequence();
    $sequence = &seq_cleanup($sequence);
    print FASTA ">$id\n$sequence\n";
   }
   undef($reference_file_obj);
  }
  confess "Nothing to align!" unless $reference_file && -s $reference_file;
  if ($extra_genes) {
   my $file_obj = new Fasta_reader($extra_genes);
   while ( my $seq_obj = $file_obj->next() ) {
    my $id = $seq_obj->get_accession() ? $seq_obj->get_accession() : 'extra';
    my $sequence = $seq_obj->get_sequence();
    $sequence = &seq_cleanup($sequence);
    print FASTA ">$id\n$sequence\n";
   }
   undef($file_obj);
  }
  close FASTA;
 }

 # Every time: process input files
 confess "Cannot find input FASTA: $file_to_align\n" unless -s $file_to_align;
 if ($remove_redund) {
  print "Checking for and removing redundant sequences...\n";
  $file_to_align = &remove_redundant_sequences($file_to_align);
 }
 $md5_aliases_file = $file_to_align . '.md5';
 if ( -s "$file_to_align.checked" && -s $md5_aliases_file && $resume ) {
  print "Parsing existing checksums using $md5_aliases_file\n";
  open( IN, $md5_aliases_file );
  while ( my $ln = <IN> ) {
   chomp($ln);
   my @data = split( "\t", $ln );
   next unless $data[1];
   my $seq_md5hash = $data[0];
   push( @reference_sequence_list, $data[1] );
   $user_alias{$seq_md5hash} = { 'id' => $data[1], 'length' => $data[2] };
   $check_hash_from_seqdb->execute($seq_md5hash);
   my $check = $check_hash_from_seqdb->fetchrow_arrayref();
   die
"Could not find the md5sum $seq_md5hash in the database, did you do something manually?! Please delete the $file_to_align.checked file and restart\n"
     if !$check;
  }
  close IN;
  print "Found " . scalar( keys %user_alias ) . " unique genes\n";
  unless ($no_checks) {
   print &mytime() . "Checking database tables for each gene...\n";
   for ( my $i = 0 ; $i < @readsets ; $i++ ) {
    foreach my $seq_md5hash ( keys %user_alias ) {
     &sqlite_init_expression_statistics( $seq_md5hash, $readsets[$i] );
     for ( my $k = $i + 1 ; $k < scalar(@readsets) ; $k++ ) {
      next if $readsets[$i] eq $readsets[$k];
      &sqlite_start_fold_change( $seq_md5hash, $readsets[$i], $readsets[$k] );
     }
    }
   }
  }

  # we assume we don't need to run the sql commands because they already exist.
 }
 else {
  print "Creating checksums for each sequence\n";

  # first time or every time .checked file is deleted.
  $do_backup = 1;
  open( MD5SUMS, ">$file_to_align.md5" );
  open( BED, ">$file_to_align.bed" ) || die($!);
  my $file_obj = new Fasta_reader($file_to_align);
  while ( my $seq_obj = $file_obj->next() ) {
   my $id = $seq_obj->get_accession();
   push( @reference_sequence_list, $id );
   my $sequence = $seq_obj->get_sequence();
   $sequence = &seq_cleanup($sequence);
   confess "No sequence for $id!\n" if ( !$sequence );
   confess "Sequence for $id contains non-ATCGN characters!\n"
     if $sequence =~ /[^ATCGN]/;
   my $seq_md5hash = md5_hex($sequence);
   $id = 'Query' unless $id;
   my $seq_length = length($sequence);
   print BED "$id\t1\t$seq_length\t$id\t0\t+\n";
   print MD5SUMS "$seq_md5hash\t$id\t$seq_length\n";
   &sqlite_add_seq_md5( $id, $seq_md5hash, $seq_length );
   confess "Some of your data has identical sequences!\n $id vs " 
     . $user_alias{$seq_md5hash}{'id'}
     . "\nUse the option -remove_redund to remove them before processing.\n"
     if $user_alias{$seq_md5hash}{'id'};
   $user_alias{$seq_md5hash} = { 'id' => $id, 'length' => $seq_length };
  }
  close MD5SUMS;
  close BED;
  open( CHECK, ">$file_to_align.checked" );
  print CHECK "Done\n";
  close CHECK;

  print &mytime() . "Preparing database tables for each gene...\n";
  for ( my $i = 0 ; $i < @readsets ; $i++ ) {
   foreach my $seq_md5hash ( keys %user_alias ) {
    &sqlite_init_expression_statistics( $seq_md5hash, $readsets[$i] );
    for ( my $k = $i + 1 ; $k < scalar(@readsets) ; $k++ ) {
     next if $readsets[$i] eq $readsets[$k];
     &sqlite_start_fold_change( $seq_md5hash, $readsets[$i], $readsets[$k] );
    }
   }
  }
  &sqlite_backup(1) if !$db_use_file && $do_backup;
 }

 if ($prepare_input_only) {
  warn "User requested to stop after preparing files\n";
  sqlite_destroy();
  exit;
 }

 # prepare file for alignments
 if ($use_kanga) {
  unless ( -s $file_to_align . '.kangax' ) {
   print "\t\t\tBuilding reference file for kanga...\r";
   &process_cmd("$kangax_exec -i $file_to_align -o $file_to_align.kangax -r $file_to_align -t $file_to_align  2> /dev/null >/dev/null"
   );
   print " Done!\n";
  }
 }
 elsif ($use_bwa) {
  unless ( -s $file_to_align . '.bwt' ) {
   print "\t\t\tBuilding reference file for bwa...\r";
   &process_cmd("$bwa_exec index $file_to_align 2> /dev/null >/dev/null");
   print " Done!\n";
  }
 }
 elsif (@use_existing_bam) {
  print "Will use user-provided BAM (read-name sorted) files\n" if !$existing_bam_is_coordsorted;
  print "Will use user-provided BAM (coordinate sorted) files\n" if $existing_bam_is_coordsorted;
 }
 else {
  unless ( -s "$file_to_align.1.bt2" ) {
   my $build_exec = $bowtie2_exec . '-build';
   print "\t\t\tBuilding reference file for bowtie2...\r";
   &process_cmd("$build_exec --offrate 1 $file_to_align $file_to_align >/dev/null 2>> $result_dir/$uid.log"
   ) unless -s "$file_to_align.rev.1.bt2";
   print " Done!\n";
  }
 }
 return $file_to_align;
}

sub starts_alignments() {
 my $file_to_align = shift;
 my $baseout       = $file_to_align;
 $baseout =~ s/.toalign$//;
 my (%already_done_alignments);
 print "\n" 
   . &mytime
   . "Starting alignments against up to "
   . scalar(@readsets)
   . " readsets\n";
 print "\tChecking database for completed alignments...\n";
 my $todo;
 ## This is very slow as it checks every readset?
 for ( my $i = 0 ; $i < (@readsets) ; $i++ ) {
  my $readset = $readsets[$i];
  if ($resume) {
   my $count = &sqlite_check_depth_readset($readset);
   if ($count) {
    $already_done_alignments{$readset} = $count;
    next;
   }
  }
  $todo++;
  my $readset_basename = fileparse($readset);
  my $alnbase          = fileparse($baseout) . '_vs_' . $readset_basename;
  my $alnout_text      = $alnbase . '.bam';
  print "\t* Paired readsets "
    . $readsets[$i] . " and "
    . $readsets2[$i]
    . " : $alnout_text\n"
    if $readsets2[$i];
  if ( !$readsets2[$i] ) {
   print "\t* Unpaired readset " . $readsets[$i] . " : $alnout_text\n";
   $single_end_readsets{ $readsets[$i] } = 1;
  }
 }

 print "\tNo alignments need to be done. Processing all existing...\n"
   if !$todo;    # SLOW
 print "\t$todo alignments need to be done...\n" if $todo;

 # testing
 return if ( !$todo && $no_checks );
 my ( %alignment_sam_files, %alignment_bam_files, $aligned_ids_hashref, %hash );

 if ($use_bwa){
	&align_bwa_index($file_to_align);
 }elsif ($use_kanga) {
	&align_kanga_index($file_to_align);
 }else{
	&align_bowtie2_index($file_to_align);
 }

 if ($contextual_alignment) {
  for ( my $i = 0 ; $i < (@readsets) ; $i++ ) {
   my $readset          = $readsets[$i];
   my $readset_basename = fileparse($readset);
   my $alnbase          = $baseout . '_vs_' . $readset_basename;
   my $alignment_bam_file = $alnbase . '.bam';
   my $alignment_sam_file = $alnbase . '.sam';

   if (!$already_done_alignments{$readset} && @use_existing_bam) {    #name or coord sorted
      my $thread = threads->create('prepare_alignment_from_existing',$file_to_align, $i,  $alnbase, $alignment_bam_file, $alignment_sam_file  );
      $alignment_thread_helper->add_thread($thread);
   }
   else {
    my $thread = threads->create('prepare_alignment',$file_to_align, $i , $alnbase, $alignment_bam_file, $alignment_sam_file );
    $alignment_thread_helper->add_thread($thread);
   }
  }

  # All alignments are done but they have not been processed yet.
  $alignment_thread_helper->wait_for_all_threads_to_complete();

  # the following cannot be threaded readily due to DB operations
  for ( my $i = 0 ; $i < (@readsets) ; $i++ ) {

   my $readset          = $readsets[$i];
   next unless $readset;
   my $readset2 = $readsets2[$i] if $readsets2[$i];
   my $readset_basename = fileparse($readset);
   my $alnbase          = $baseout . '_vs_' . $readset_basename;
   my $alignment_bam_file = $alnbase . '.bam';
   my $alignment_sam_file = $alnbase . '.sam';

   my ($rpkm_hashref, $eff_counts_hashref);

   if ( $already_done_alignments{$readset} ) {
    my $salmon_results = $alignment_bam_file . ".salmon.results";
    die "Supposedly existing alignment $alignment_bam_file does not exist. Any chance you're using an old database?\n" unless -s $alignment_bam_file;
     # can't thread this as it does DB operations
    &process_salmon_bias( $salmon_results, $readset );
   }else{
    if ($perform_bias_correction){
	my $namesort = -s $alignment_sam_file . '.namesorted' ? $alignment_sam_file . '.namesorted' : $alignment_bam_file . '.namesorted';
       ($alignment_bam_file, $rpkm_hashref, $eff_counts_hashref ) = &perform_correct_bias( $alignment_bam_file, $file_to_align, $readset, $namesort );
    }
    if ($alignment_bam_file && -s $alignment_bam_file){
      &process_alignments( $alignment_bam_file, $readset, $readset2 ) if $readset2;
      &process_alignments( $alignment_bam_file, $readset ) if !$readset2;
    }else{
	warn  "Alignment for $readset did not complete. Will skip.\n";
	print LOG "Alignment for $readset did not complete. Will skip.\n";
	delete($library_metadata{$readsets[$i]});
	delete($readsets[$i]);
	delete($readsets2[$i]) if $readsets2[$i];
    }
   }
 } 

  for ( my $i = 0 ; $i < (@readsets) ; $i++ ) {
   next unless $readsets[$i];

   print "    Processing Depth of Coverage "
     . ( $i + 1 ) . " of "
     . scalar(@readsets)
     . " readsets                                         \r";

   my $readset          = $readsets[$i];
   my $readset2 = $readsets2[$i] if $readsets2[$i];
   my $readset_basename = fileparse($readset);
   my $alnbase          = $baseout . '_vs_' . $readset_basename;
   my $alignment_bam_file = $alnbase . '.bam';
   my $alignment_sam_file = $alnbase . '.sam';

   next if $only_alignments || ($resume && -s $alignment_bam_file . ".depth.completed");

   $aligned_ids_hashref = &process_depth_of_coverage( $file_to_align, $readset, $alignment_bam_file );
   $alignment_sam_files{$readset} = $alignment_sam_file;
   $alignment_bam_files{$readset} = $alignment_bam_file;

   # so if $aligned_ids_hashref from process_depth above
   # is not equal to the references added, then re-enter depth data
   my $do_backup;
   if ( keys %{$aligned_ids_hashref} != scalar(@reference_sequence_list) ) {
    foreach my $id (@reference_sequence_list) {
     next if $aligned_ids_hashref->{$id};
     my $seq_md5hash = &sqlite_get_md5($id);

     if ( &sqlite_get_depth_data( $seq_md5hash, $readset, 1 ) ) {
      &sqlite_add_depth_data( $seq_md5hash, $readset, \%hash )
        if !$resume && !$debug;
      $do_backup++ if !$resume && !$debug;
     }
     else {
      &sqlite_add_depth_data( $seq_md5hash, $readset, \%hash );
      $do_backup++;
     }

    }
   }
   &sqlite_backup() if !$db_use_file && $do_backup;
  }
  print "\n";
  &perform_kangade( $file_to_align, \%alignment_sam_files,
                    \%alignment_bam_files )
    unless !$do_kangade;
 }
 else {

#not contextual
# due to kangade only one aln file so we have re-align if it has not been tested for a readset
  my $new_file_to_align = $file_to_align . '_unaligned';
  open( OUTSEQ, ">$new_file_to_align" );
  my (%unaligned);
  for ( my $i = 0 ; $i < (@readsets) ; $i++ ) {
   next unless $readsets[$i];

   my $file_obj = new Fasta_reader($file_to_align);
   while ( my $seq_obj = $file_obj->next() ) {
    my $id          = $seq_obj->get_accession();
    my $seq_md5hash = &sqlite_get_md5($id);
    my $aln_exists  = &sqlite_get_depth_data( $seq_md5hash, $readsets[$i], 1 );
    if ( !$aln_exists ) {
     print OUTSEQ ">$id\n" . &wrap_text( $seq_obj->get_sequence() ) . "\n"
       if !$unaligned{$seq_md5hash};
     $unaligned{$seq_md5hash} = 1;
    }
    else {
     warn "Sequence $id has already been aligned against "
       . $readsets[$i]
       . " ($seq_md5hash)\n"
       if $debug;
    }
   }
  }
  close OUTSEQ;
  link( $file_to_align . '.bed', $new_file_to_align . '.bed' );
  if ( -s $new_file_to_align ) {
   for ( my $i = 0 ; $i < (@readsets) ; $i++ ) {

   my $readset          = $readsets[$i];
   next unless $readset;
   my $readset2 = $readsets2[$i] if $readsets2[$i];
   my $readset_basename = fileparse($readset);
   my $alnbase          = $baseout . '_vs_' . $readset_basename;
   my $alignment_bam_file = $alnbase . '.bam';
   my $alignment_sam_file = $alnbase . '.sam';
   my ($rpkm_hashref, $eff_counts_hashref);

    #  print '.' x ( $i + 1 ) . "\r";
    print "\tProcessing "
      . ( $i + 1 )
      . " out of "
      . scalar(@readsets)
      . " readsets                \r";
      &prepare_alignment( $new_file_to_align, $i , $alnbase, $alignment_bam_file, $alignment_sam_file );

   ($alignment_bam_file, $rpkm_hashref, $eff_counts_hashref ) = &perform_correct_bias( $alignment_bam_file, $file_to_align, $readset, $alignment_bam_file . '.namesorted' )
   if ($perform_bias_correction);

   if (!$alignment_bam_file || !-s $alignment_bam_file){
        warn "Alignment for $readset did not complete. Will skip.\n";
        print LOG "Alignment for $readset did not complete. Will skip.\n";
	delete($library_metadata{$readsets[$i]});
        delete($readsets[$i]);
        delete($readsets2[$i]) if $readsets2[$i];
    }

   &process_alignments( $alignment_bam_file, $readset, $readset2 ) if $readset2;
   &process_alignments( $alignment_bam_file, $readset ) if !$readset2;

    next if $only_alignments  || ($resume && -s $alignment_bam_file . ".depth.completed");

    #too slow:
    $aligned_ids_hashref =
      &process_depth_of_coverage( $new_file_to_align, $readsets[$i],
                                  $alignment_bam_file );
    $alignment_sam_files{ $readsets[$i] } = $alignment_sam_file;
    $alignment_bam_files{ $readsets[$i] } = $alignment_bam_file;
    if ( keys %{$aligned_ids_hashref} != scalar(@reference_sequence_list) ) {
     foreach my $id (@reference_sequence_list) {
      next if $aligned_ids_hashref->{$id};
      my $seq_md5hash = &sqlite_get_md5($id);
      unless ( &sqlite_get_depth_data( $seq_md5hash, $readsets[$i], 1 ) ) {
       &sqlite_add_depth_data( $seq_md5hash, $readsets[$i], \%hash );
      }
     }
    }
   }
   print "\n";
  }
  &perform_kangade( $new_file_to_align, \%alignment_sam_files,
                    \%alignment_bam_files )
    unless !$do_kangade;
 }
}

sub prepare_alignment_from_existing() {
 my ( $file_to_align, $i, $baseout, $alignment_bam_file,  $alignment_sam_file ) = @_;
 confess "Given BAM alignment file does not exist: "
   . $use_existing_bam[$i] . "\n"
   unless $use_existing_bam[$i] && -s $use_existing_bam[$i];
 my $readset          = $readsets[$i];
 my $readset2         = $readsets2[$i];
 my $readset_basename = fileparse($readset);
 my $alnbase          = $baseout . '_vs_' . $readset_basename;
 print "\n\n" . &mytime
   . "Aligning: Using user-provided alignment ".$use_existing_bam[$i]." for $readset_basename\n";
 print LOG "\n" . &mytime
   . "Aligning: Using user-provided alignment ".$use_existing_bam[$i]." for $readset_basename\n";

 if (!$existing_bam_is_coordsorted){
 # namesorted bam
 # hard link:
 link( $use_existing_bam[$i], "$alignment_bam_file.namesorted" )
   unless -s "$alignment_bam_file.namesorted";

 # in case hard link cannot occur:
 symlink( $use_existing_bam[$i], "$alignment_bam_file.namesorted" )
   unless -s "$alignment_bam_file.namesorted";
 confess "Cannot link "
   . $use_existing_bam[$i]
   . " as $alignment_bam_file.namesorted\n"
   unless -s "$alignment_bam_file.namesorted";
  &process_cmd("$samtools_exec sort -T $sort_tmp -o $alignment_bam_file -@ $samtools_threads -m $sam_sort_memory "
    . $use_existing_bam[$i]
    . "  2>/dev/null" ) unless -s $alignment_bam_file;
 }else{
 # coord sorted
 # hard link:
 link( $use_existing_bam[$i], "$alignment_bam_file" )
   unless -s "$alignment_bam_file";
 # in case hard link cannot occur:
 symlink( $use_existing_bam[$i], "$alignment_bam_file" )
   unless -s "$alignment_bam_file";
 confess "Cannot link "
   . $use_existing_bam[$i]
   . " as $alignment_bam_file\n"
   unless -s "$alignment_bam_file";
  unless ( -s "$alignment_bam_file.namesorted"){
    # samtools sort -f does'nt work with insuffiecient memory
    &process_cmd(   "$samtools_exec sort -T $sort_tmp -o $alignment_bam_file.namesorted -n -@ $samtools_threads -m $sam_sort_memory "
               . $use_existing_bam[$i] . " 2>/dev/null"
                );
  }
  die "Cannot find $alignment_bam_file.namesorted \n" unless -s "$alignment_bam_file.namesorted";
 }
 unlink($alignment_sam_file) if -s $alignment_sam_file;
   
}

sub prepare_alignment() {
 my ( $file_to_align, $i, $baseout, $alignment_bam_file, $alignment_sam_file ) = @_;
 my $readset  = $readsets[$i];
 my $readset2 = $readsets2[$i];
 &perform_alignments( $file_to_align, $readset, $readset2, $baseout, $alignment_bam_file, $alignment_sam_file );

}

sub perform_alignments() {
 my $file_to_align    = shift;
 my $readset          = shift;
 my $readset_basename = fileparse($readset);
 my $readset2         = shift;
 my $baseout = shift;
 my $bam = shift;
 my $sam = shift;

 confess "File $readset not found!" if (!-s $readset);

 print "\n" . &mytime
   . "Aligning: Performing alignments of $file_to_align vs $readset_basename\n";
 my $readset_time = stat($readset)->mtime;
 die "Cannot get modification time for $readset" unless $readset_time;
 my $readset2_time = $readset2 ? stat($readset2)->mtime : int(0);
 print LOG &mytime() . "Processing $readset as $bam\n";

 # check if alignments need to be recalculated
 if (
      -s $bam
      && (    stat($bam)->mtime > $readset_time
           && stat($bam)->mtime > $readset2_time )
   )
 {
  return ( $bam, $sam );
 }

 if ($use_bwa) {
  unlink("$baseout.sai")
    if ( -s "$baseout.sai"
         && ( stat("$baseout.sai")->mtime < $readset_time ) );
  unlink("$baseout.2.sai")
    if ( -s "$baseout.2.sai"
         && ( stat("$baseout.2.sai")->mtime < $readset2_time ) );
 }
 unlink($bam)
   if (
        -s $bam
        && (    stat($bam)->mtime < $readset_time
             || stat($bam)->mtime < $readset2_time )
   );
 unlink($sam)
   if (
        -s $sam
        && (    stat($sam)->mtime < $readset_time
             || stat($sam)->mtime < $readset2_time )
   );
 unless ( -s $bam && ( -s $bam ) > 1000 ) {
  if ($use_bwa) {
   &align_bwa( $baseout, $file_to_align, $readset, $readset2, $bam, $sam );
  }
  elsif ($use_kanga) {
   &align_kanga( $baseout, $file_to_align, $readset, $readset2, $bam, $sam );
  }
  else {
   &align_bowtie2( $baseout, $file_to_align, $readset, $readset2, $bam, $sam );
  }
  confess "Could not produce SAM file for $readset\n" unless -s $sam;
  &process_cmd("$samtools_exec view -u $sam 2>/dev/null|$samtools_exec sort -T $sort_tmp -o $bam -@ $samtools_threads -m $sam_sort_memory -  2>/dev/null"  ) unless -s $bam;
 }
 confess "Could not convert to BAM file ($bam) for $readset\n" unless -s $bam;

 &process_cmd("$samtools_exec index $bam")
   unless -s $bam . '.bai' && ( -s $bam . '.bai' ) > 200;

 return ( $bam, $sam );
}

sub process_depth_of_coverage($$$) {
 my $file_to_align = shift;
 my $readset       = shift;
 my $bam           = shift;
 return if ( -s $bam . ".depth.completed" && $resume );

 my $readset_metadata = &sqlite_get_readset_metadata($readset);
 my $readset_name     = $readset_metadata->{'alias'};
 print LOG &mytime() . "Processing $readset_name depth from $bam\n";
 print &mytime() . "Processing $readset_name depth from $bam\n";
 my $tmp_depth_file = $bam . ".depth";
 &process_cmd("$samtools_exec depth -s -q 10 -Q 30 $bam > $tmp_depth_file 2>/dev/null")
   unless -s $tmp_depth_file;

 unless ( -s $tmp_depth_file > 100 ) {
  warn "Could not produce depth file for $bam. Skipping\n";
  print LOG "Could not produce depth file for $bam. Skipping\n";
  unlink($tmp_depth_file);
  return;
 }
 my ( %hash, %aligned_ids );

 my $timer_counter = int(0);
 my $timer         = new Time::Progress;
 $timer->attr( min => 0, max => -s $tmp_depth_file );
 open( DEPTH, $tmp_depth_file ) || die($!);
 my $previous_id;
 while ( my $ln = <DEPTH> ) {
  $timer_counter += length($ln);
  print $timer->report( "eta: %E min, %40b %p\r", $timer_counter )
    if ( $timer_counter % 1000000 == 0 );
  chomp($ln);
  my @data = split( "\t", $ln );
  if ( $previous_id && $previous_id ne $data[0] ) {
   my $seq_md5hash = &sqlite_get_md5($previous_id);
   $aligned_ids{$previous_id} = 1;
   $previous_id = $data[0];

# do not add if they already exist (unless it is a contextual alignment without debug, in which case overwrite)
   if ( &sqlite_get_depth_data( $seq_md5hash, $readset, 1 ) ) {
    &sqlite_add_depth_data( $seq_md5hash, $readset, \%hash )
      if ( ( $contextual_alignment && !$debug ) || !$resume );
   }
   else {
    &sqlite_add_depth_data( $seq_md5hash, $readset, \%hash );
   }
   %hash = ();
  }
  $hash{ $data[1] } = $data[2];
  $previous_id = $data[0];
 }

 # process last ID
 my $seq_md5hash = &sqlite_get_md5($previous_id);
 $aligned_ids{$previous_id} = 1;
 if ( &sqlite_get_depth_data( $seq_md5hash, $readset, 1 ) ) {
  &sqlite_add_depth_data( $seq_md5hash, $readset, \%hash )
    if ( ( $contextual_alignment && !$debug ) || !$resume );
 }
 else {
  &sqlite_add_depth_data( $seq_md5hash, $readset, \%hash );
 }
 close DEPTH;
 open( OUT, ">$bam.depth.completed" );
 print OUT &mytime()."\n";
 close OUT;


 unlink($tmp_depth_file) unless $debug;
 return \%aligned_ids;
}

sub perform_kangade() {
 my $file_to_align              = shift;
 my $alignment_samfiles_hashref = shift;
 my $alignment_bamfiles_hashref = shift;
 my $baseout                    = $file_to_align . "_kangade_";
 my $kanga_de_stats_outfile     = "$file_to_align.kangade.stats.tsv";
 if ( !-s $kanga_de_stats_outfile . '.completed' ) {
  print "\n"
    . &mytime()
    . "kangade: Calculating differences across gene length for each readsets pair\n";
  open( STATS_KANGADE, ">$kanga_de_stats_outfile" ) || die($!);
  print STATS_KANGADE
"Checksum\tGene\tReadset1\tReadset2\tClassification\tScore\tDECntsScore\tPearsonScore\tObsPearson\tObsFoldChange\n";

# iterate through each reference sequence and process each pair of readsets. first found those that need to be processed.
  my ( %md5_not_to_process, %files_to_delete );

  #  my (%md5_to_process);
  if ( !$contextual_alignment ) {
   foreach my $id (@reference_sequence_list) {
    my $seq_md5hash = &sqlite_get_md5($id);
    for ( my $i = 0 ; $i < ( scalar(@readsets) - 1 ) ; $i++ ) {
     next unless $readsets[$i];
     my $readset_C = $readsets[$i];
     for ( my $k = $i + 1 ; $k < scalar(@readsets) ; $k++ ) {
      next unless $readsets[$k];
      my $readset_E = $readsets[$k];
      my $exists =
        &sqlite_check_kangade_analysis( $seq_md5hash, $readset_C, $readset_E );
      next if !$exists;
      my $readset_metadata_C = &sqlite_get_readset_metadata($readset_C);
      my $readset_metadata_E = &sqlite_get_readset_metadata($readset_E);
      my $readset_name_C     = $readset_metadata_C->{'alias'};
      my $readset_name_E     = $readset_metadata_E->{'alias'};
      $exists->{'ObsFoldChange'} =
        sprintf( "%.2f", $exists->{'ObsFoldChange'} );

# Classification,Score,DECntsScore,PearsonScore,CtrlUniqueLoci,ExprUniqueLoci,CtrlExprLociRatio,PValueMedian,PValueLow95,PValueHi95,TotCtrlCnts,TotExprCnts,TotCtrlExprCnts,ObsFoldChange,FoldMedian,FoldLow95,FoldHi95,ObsPearson,PearsonMedian,PearsonLow95,PearsonHi95
      print STATS_KANGADE $seq_md5hash . "\t" 
        . $id
        . "\t$readset_name_C\t$readset_name_E\t"
        . $exists->{'Classification'} . "\t"
        . $exists->{'Score'} . "\t"
        . $exists->{'DECntsScore'} . "\t"
        . $exists->{'PearsonScore'} . "\t"
        . $exists->{'ObsPearson'} . "\t"
        . $exists->{'ObsFoldChange'} . "\n";
      $add_kangade_fold_change->execute( $exists->{'ObsFoldChange'},
                                         $seq_md5hash, $readset_C, $readset_E );
      $fold_changes{$seq_md5hash}{'kangade'}{$readset_name_C}{$readset_name_E} =
        $exists->{'ObsFoldChange'};
      &sqlite_add_kangade_expression_statistics( $seq_md5hash, $readset_C,
                                                 $exists->{'TotCtrlCnts'} );
      &sqlite_add_kangade_expression_statistics( $seq_md5hash, $readset_E,
                                                 $exists->{'TotExprCnts'} )
        if $i == scalar(@readsets) - 2;
      $md5_not_to_process{$seq_md5hash} = 1;
     }
    }
   }
  }

# now perform kangade against all data but parse only those that were not stored.
# i know that is it not the most efficient...
  for ( my $i = 0 ; $i < ( scalar(@readsets) - 1 ) ; $i++ ) {
   next unless $readsets[$i];
   my $readset_C          = $readsets[$i];
   my $readset_metadata_C = &sqlite_get_readset_metadata($readset_C);
   my $readset_name_C     = $readset_metadata_C->{'alias'};
   my $readset_reads_C    = $readset_metadata_C->{'total_reads'};
   my $readset_basename_C = fileparse($readset_C);
   print &mytime() . "\tProcessing $readset_name_C\n";
   for ( my $k = $i + 1 ; $k < scalar(@readsets) ; $k++ ) {
    next unless $readsets[$k];
    my $readset_E          = $readsets[$k];
    my $readset_metadata_E = &sqlite_get_readset_metadata($readset_E);
    my $readset_name_E     = $readset_metadata_E->{'alias'};
    my $readset_basename_E = fileparse($readset_E);
    my $out = $baseout . $readset_basename_C . '_vs_' . $readset_basename_E;
    unless ( -s "$out.stats.csv" ) {
     my $readset_reads_E = $readset_metadata_E->{'total_reads'};
     my $library_ratio =
       sprintf( "%.2f", ( $readset_reads_C / $readset_reads_E ) );
     my $alignment_sam_file_C = $alignment_samfiles_hashref->{$readset_C};
     if ( ( !$alignment_sam_file_C || !-s $alignment_sam_file_C )
          && $alignment_bamfiles_hashref->{$readset_C} )
     {
      $alignment_sam_file_C =
        $alignment_bamfiles_hashref->{$readset_C} . '.sam';
      &process_cmd( "$samtools_exec view -h -o $alignment_sam_file_C "
                    . $alignment_bamfiles_hashref->{$readset_C} )
        unless -s $alignment_sam_file_C;
      $files_to_delete{$alignment_sam_file_C} = 1;
     }
     elsif (    !$alignment_sam_file_C
             || !-s $alignment_sam_file_C )
     {
      warn("Cannot find alignment SAM file for $readset_C\n");
      next;
     }
     my $alignment_sam_file_E = $alignment_samfiles_hashref->{$readset_E};
     if ( ( !$alignment_sam_file_E || !-s $alignment_sam_file_E )
          && $alignment_bamfiles_hashref->{$readset_E} )
     {
      $alignment_sam_file_E =
        $alignment_bamfiles_hashref->{$readset_E} . '.sam';
      &process_cmd( "$samtools_exec view -h -o $alignment_sam_file_E "
                    . $alignment_bamfiles_hashref->{$readset_E} )
        unless -s $alignment_sam_file_E;
      $files_to_delete{$alignment_sam_file_E} = 1;
     }
     elsif (    !$alignment_sam_file_E
             || !-s $alignment_sam_file_E )
     {
      warn("Cannot find alignment SAM file for $readset_E\n");
      next;
     }
     &process_cmd(
"$kangade_exec -r 0 -T $kanga_threads -t 3 -g $file_to_align.bed -i $alignment_sam_file_C -I $alignment_sam_file_E -f3 -F $out.log -o $out.stats.csv -O $out.bins.csv >/dev/null 2>/dev/null"
     ) if $kangade_exec;

     # no longer need SAMs
     unless ($debug) {
      foreach my $readset (@readsets) {
       unlink( $alignment_samfiles_hashref->{$readset} )
         if $alignment_samfiles_hashref->{$readset};
      }
      foreach my $file ( keys %files_to_delete ) {
       unlink($file);
      }
     }
    }
    unless ( -s "$out.stats.csv" ) {
     warn "kangade failed for $out\n";
     next;
    }
    my $csv = Text::CSV_XS->new()
      || confess "Cannot use CSV: " . Text::CSV_XS->error_diag();
    open( my $fh, "$out.stats.csv" ) || die($!);

    #my $h = <$fh>;chomp($h);$h =~ s/"//g;my @headers = split( "\t", $h );
    $csv->column_names( $csv->getline($fh) );
    while ( my $row = $csv->getline_hr($fh) ) {
     my $seq_md5hash = &sqlite_get_md5( $row->{'Feat'} );
     next if $md5_not_to_process{$seq_md5hash};
     $row->{'ObsFoldChange'} = sprintf( "%.2f", $row->{'ObsFoldChange'} );
     print STATS_KANGADE $seq_md5hash . "\t"
       . $row->{'Feat'}
       . "\t$readset_C\t$readset_E\t"
       . $row->{'Classification'} . "\t"
       . $row->{'Score'} . "\t"
       . $row->{'DECntsScore'} . "\t"
       . $row->{'PearsonScore'} . "\t"
       . $row->{'ObsPearson'} . "\t"
       . $row->{'ObsFoldChange'} . "\n";
     $add_kangade_fold_change->execute( $row->{'ObsFoldChange'},
                                        $seq_md5hash, $readset_C, $readset_E )
       unless $contextual_alignment;
     $fold_changes{$seq_md5hash}{'kangade'}{$readset_name_C}{$readset_name_E} =
       $row->{'ObsFoldChange'};
     &sqlite_add_kangade_expression_statistics( $seq_md5hash, $readset_C,
                                                $row->{'TotCtrlCnts'} );
     &sqlite_add_kangade_expression_statistics( $seq_md5hash, $readset_E,
                                                $row->{'TotExprCnts'} )
       if $i == scalar(@readsets) - 2;
     my $de_exists =
       &sqlite_check_kangade_analysis( $seq_md5hash, $readset_C, $readset_E );
     next if $de_exists;
     $add_to_kangade_analysis->execute(
                                        $seq_md5hash,
                                        $readset_C,
                                        $readset_E,
                                        $row->{'Classification'},
                                        $row->{'Score'},
                                        $row->{'DECntsScore'},
                                        $row->{'PearsonScore'},
                                        $row->{'CtrlUniqueLoci'},
                                        $row->{'ExprUniqueLoci'},
                                        $row->{'CtrlExprLociRatio'},
                                        $row->{'PValueMedian'},
                                        $row->{'PValueLow95'},
                                        $row->{'PValueHi95'},
                                        $row->{'TotCtrlCnts'},
                                        $row->{'TotExprCnts'},
                                        $row->{'TotCtrlExprCnts'},
                                        $row->{'ObsFoldChange'},
                                        $row->{'FoldMedian'},
                                        $row->{'FoldLow95'},
                                        $row->{'FoldHi95'},
                                        $row->{'ObsPearson'},
                                        $row->{'PearsonMedian'},
                                        $row->{'PearsonLow95'},
                                        $row->{'PearsonHi95'}
     );
    }
    close IN;
   }
  }
  close STATS_KANGADE;
  open( OUT, ">$kanga_de_stats_outfile.completed" );
  print OUT &mytime();
  close OUT;
  &sqlite_backup() unless $db_use_file;
 }
 elsif ($contextual_alignment) {
  print "\n" . &mytime() . "kangade: Parsing existing data...\n";
  open( IN, $kanga_de_stats_outfile );
  while ( my $ln = <IN> ) {
   chomp($ln);
   my @data = split( "\t", $ln );
   next unless $data[9];
   $fold_changes{ $data[0] }{'kangade'}{ $data[2] }{ $data[3] } = $data[9];
  }
  close IN;
 }
}

sub perform_readset_metadata($$) {
 my $readset      = shift;
 my $readset2     = shift;
 my $do_backup    = 0;
 my $library_size = int(0);
 my $readset_name =
     $library_aliases{$readset}
   ? $library_aliases{$readset}
   : fileparse($readset);
 my $is_paired =
   $paired_readset_lookup{$readset} ? $paired_readset_lookup{$readset} : '0';
 $get_readset->execute($readset);
 my $result = $get_readset->fetchrow_hashref();

 if ( $result && $result->{'total_reads'} ) {

  #already exists
  #update alias if needed
  if ( !$result->{'alias'}
       || ( $result->{'alias'} && $result->{'alias'} ne $readset_name ) )
  {
   $update_readset->execute( $is_paired, $readset_name, $readset );
   $get_readset->execute($readset);
   $result    = $get_readset->fetchrow_hashref();
   $do_backup = 1;
  }
  $library_metadata{$readset_name}{'readlength_median'} =
    $result->{'readlength_median'};
  $library_metadata{$readset_name}{'nonnormalized_library_size'} =
    $result->{'total_reads'};

  return ( $result->{'alias'}, $result->{'total_reads'}, $do_backup );
 }
 $do_backup = 1;
 print "Counting reads in PE readsets $readset / $readset2\n" if $readset2;
 print "Counting reads in SE readset $readset\n"              if !$readset2;

 # get library sizes
 if ( $read_format eq 'fastq' ) {
  if ( $readset =~ /.bz2$/ ) {
   $library_size = `$bunzip2_exec -dck $readset| wc -l `;
   chomp($library_size);
   $library_size += `$bunzip2_exec -dck $readset2| wc -l `
     if $readset2;
  }
  else {
   $library_size = `wc -l < $readset`;
   chomp($library_size);
   $library_size += `wc -l < $readset2` if $readset2;
  }
  $library_size /= 4;
  print LOG "Reads in $readset / $readset2 PE files: $library_size\n"
    if $readset2;
  print LOG "Reads in $readset SE file: $library_size\n"
    if !$readset2;
 }
 elsif ( $read_format eq 'bam' ) {
  confess "BAM $readset not found\n" unless -s $readset;
  &process_cmd("$samtools_exec index $readset")
    unless -s $readset . '.bai';
  confess "Readset $readset cannot be indexed" unless -s $readset . '.bai';
  my $d = `$samtools_exec idxstats $readset|grep '^\*'`;
  $d =~ /(\d+)$/;
  $library_size = $1;
  if ($readset2) {
   $d = 0;
   $d = `$samtools_exec idxstats $readset2|grep '^\*'`;
   $d =~ /(\d+)$/;
   $library_size += $1;
  }
  print LOG "Reads in $readset file: $library_size\n"
    if !$readset2;
  print LOG "Reads in $readset / $readset2 PE files: $library_size\n"
    if $readset2;
 }
 confess "Cannot get library size for $readset"
   unless $library_size && $library_size > 0;

 my $readlength_median = &get_readlength_median($readset);

 $library_metadata{$readset_name}{'readlength_median'} = $readlength_median;
 $library_metadata{$readset_name}{'nonnormalized_library_size'} = $library_size;

 &sqlite_add_readset_metadata( $readset,      $readset_name,
                               $library_size, $readlength_median );
 return ( $readset_name, $library_size, $do_backup );
}

sub get_readlength_median() {
 my $file = shift;
 my @array;
 if ( $read_format eq 'fastq' ) {
  my @reads;
  if ( $file =~ /\.bz2/ ) {
   @reads = `bunzip2 -kc $file | head -n 400000 `;
  }
  else {
   @reads = `head -n 400000 $file`;
  }
  for ( my $i = 1 ; $i < @reads ; $i += 4 ) {
   my $seq = $reads[$i] || last;
   chomp($seq);
   push( @array, length($seq) );
  }
 }
 elsif ( $read_format eq 'bam' ) {
  my @reads = `$samtools_exec view $file|head -n 100000 | cut -f 10`;
  foreach my $seq (@reads) {
   chomp($seq);
   push( @array, length($seq) );
  }
 }
 my $readlength_median = sprintf( "%.2f", &median( \@array ) );
 confess "Cannot get read size for $file\n" if !$readlength_median;

 return $readlength_median;
}

sub median() {
 my $array_ref = shift;
 my @sorted = sort { $a <=> $b } @{$array_ref};
 return $sorted[ int( @sorted / 2 ) ];
}

sub align_bowtie2_index(){
 my $file_to_align = shift;
 unless ( -s "$file_to_align.1.bt2" ) {
  print "\t\t\tBuilding reference file for bowtie...\r";
  &process_cmd("$bowtie2_exec -build --offrate 1 $file_to_align $file_to_align >/dev/null 2>> $result_dir/$uid.log"  );
  print " Done!\n";
 }

}

sub align_bowtie2() {
 my ( $baseout, $file_to_align, $readset, $readset2, $bam, $sam ) = @_;
 my $readgroup = '@RG\tID:' . $readset;
 $readgroup =~ s/\.bz2$//;
 $readgroup =~ s/\.fastq$//;
 $readgroup =~ s/\.bam$//;
 if ($readset2) {
  if ( $read_format eq 'fastq' ) {
   my $qformat = &check_fastq_format($readset);
   if ($qformat eq 'fasta'){
   	&process_cmd("$bowtie2_exec --reorder -a -f --threads $alignment_threads --rdg 6,5 --rfg 6,5 --score-min L,-0.6,-0.4 --no-discordant -X $fragment_max_length --no-mixed --no-unal -x $file_to_align -1 $readset -2 $readset2 -S $sam 2> $baseout.log"   ) unless -s $sam;
   }else{
	$qformat = '--'.$qformat;
      # NB --all is -a but --all is a bug
   	&process_cmd("$bowtie2_exec --reorder $qformat -a -q --threads $alignment_threads --rdg 6,5 --rfg 6,5 --score-min L,-0.6,-0.4 --no-discordant -X $fragment_max_length --no-mixed --no-unal -x $file_to_align -1 $readset -2 $readset2 -S $sam 2> $baseout.log"   ) unless -s $sam;
   }
  }
  else {
   confess "Sorry, paired end Bowtie does not work with BAM files\n";

#&process_cmd("$bamtools_exec convert -in $readset -format fastq | $bowtie2_exec --end-to-end --fast -qap $alignment_threads -x $file_to_align -U - -S $sam");
#&process_cmd("$bamtools_exec convert -in $readset2 -format fastq | $bowtie2_exec --end-to-end --fast -qap $alignment_threads -x $file_to_align -U - >> $sam");
  }
 }
 else {

  # single end
  if ( $read_format eq 'fastq' ) {
   my $qformat = &check_fastq_format($readset);
   if ($qformat eq 'fasta'){
	&process_cmd("$bowtie2_exec --reorder -a -f --rdg 6,5 --rfg 6,5 --score-min L,-0.6,-0.4 --no-unal --threads $alignment_threads -x $file_to_align -U $readset -S $sam 2> $baseout.log >/dev/null" ) unless -s $sam;
   }else{
	$qformat = '--'.$qformat;
	&process_cmd("$bowtie2_exec --reorder $qformat -a -q --rdg 6,5 --rfg 6,5 --score-min L,-0.6,-0.4 --no-unal --threads $alignment_threads -x $file_to_align -U $readset -S $sam 2> $baseout.log >/dev/null" ) unless -s $sam;
   }
  }
  else {
   &process_cmd("$bamtools_exec convert -in $readset -format fastq | $bowtie2_exec --reorder -a --rdg 6,5 --rfg 6,5 --score-min L,-0.6,-0.4 --no-unal -q --threads $alignment_threads -x $file_to_align -U - -S $sam 2> $baseout.log"
   ) unless -s "$baseout.log";
  }
 }
 confess "Alignment for $file_to_align vs $readset failed\n" unless -s $sam;
 # sometime bowtie2 doesn't actually sort them
 &namesort_sam($sam);
}


sub namesort_sam(){
 my $sam = shift;
 warn "No SAM file or not existing: $sam" unless $sam && -s $sam;
 my $out = "$sam.namesorted";
 return $out if -s $out;

 # test if it is a sam; print headers to which we concatanate results later
 &process_cmd("$samtools_exec view -H $sam > $out 2> /dev/null");
 confess "Can't produce $out. Is $sam a SAM file?\n" unless -s $out;

 # one thing we have to do is remove reads that have an unmapped mate: this is a requirement for salmon etc
 my $cmd = $perform_bias_correction
	 ? "$samtools_exec view -q 30 -F 12 $sam 2> /dev/null| $sort_exec -T $sort_tmp -S $sort_memory_sub -nk4,4|$sort_exec -T $sort_tmp -S $sort_memory_sub -s -k3,3|$sort_exec -T $sort_tmp -S $sort_memory_sub -s -k1,1 >> $out"
	 : "$samtools_exec view -q 30 $sam 2> /dev/null| $sort_exec -T $sort_tmp -S $sort_memory_sub -nk4,4|$sort_exec -T $sort_tmp -S $sort_memory_sub -s -k3,3|$sort_exec -T $sort_tmp -S $sort_memory_sub -s -k1,1 >> $out";
 &process_cmd($cmd);
 unlink($sam);
 #TODO CHECK
 &process_cmd("ln $out $sam");
# link( $result_dir.fileparse($out), $result_dir.fileparse($sam));
 return $sam;
}

sub fix_check_kanga() {

 #	####temp!!! debug
 my $sam = shift;

 #	my $errors;
 #	confess unless $in && -s $in;
 #	open (IN,$in);
 #	open (OUT,">$in.fixed");
 #	while (my $ln=<IN>){
 #       	if ($ln=~/^@/){
 #              	print OUT $ln;next;
 #	        }
 #       	my @data = split("\t",$ln);
 #	        $data[5]=~/(\d+)(M)/;
 #		next unless $2;
 #       	my $match = $1 ? int($1) : int(0);
 #	        if ($match != length($data[9])){
 #       	        $errors++;
 #              	$data[5] = length($data[9]).'M';
 #	        }
 #       	print OUT join("\t",@data);
 #	}
 #	close IN;
 #	close OUT;
 #	unlink($in);
 &namesort_sam($sam);
}

sub align_kanga_index(){
	my $file_to_align = shift;
	unless (-s $file_to_align . '.kangax'){
		print "\t\t\tBuilding reference file for kanga...\r";
	 	&process_cmd("$kangax_exec -i $file_to_align -o $file_to_align.kangax -r $file_to_align -t $file_to_align -T $threads 2> /dev/null >/dev/null" );
		print " Done!                                   \n";
	}

}

sub align_kanga() {
 my ( $baseout, $file_to_align, $readset, $readset2, $bam, $sam ) = @_;
 if ( -s $sam ) {
   &process_cmd("ln $sam $baseout.sam.namesorted");
 # link( $result_dir.fileparse($sam), $result_dir.fileparse($baseout) . ".sam.namesorted" );
  return;
 }
 if ( $read_format eq 'fastq' ) {
  if ( $readset =~ /\.bz2$/ ) {
   confess "Kanga does not support .bz2 FASTQ files (only gz).\n";
  }
  my $qformat =  &check_fastq_format($readset);
  if ($qformat eq 'phred33' ){$qformat = ' -g0 ';}
  elsif($qformat eq 'phred64' ){$qformat = ' -g1 ';}
  elsif($qformat eq 'fasta' ){$qformat = ' -g3 ';}
  if ($readset2) {
   ## Paired END:
#  Current version did not allow multiple mappings.In case it ever does:
#   &process_cmd("$kanga_exec --rptsamseqsthres=5000000 -I $file_to_align.kangax -m 0 $qformat -R 100 -N -r 5 -X -M 5 -U 2 -i $readset -u $readset2 -o $sam -F $sam.log -T $threads -d 50 -D $fragment_max_length >/dev/null"    );
#   &fix_check_kanga("$sam");
   unless ( -s "$sam.1" ) {
    &process_cmd("$kanga_exec -I $file_to_align.kangax -m 0 $qformat -R 100 -N -r 5 -X -M 5 -i $readset -o $sam.1 -F $sam.1.log -T $alignment_threads -w $sam.sqlite -W $readset >/dev/null"
    );
    &fix_check_kanga("$sam.1");
   }
   unless ( -s "$sam.2" ) {
    &process_cmd("$kanga_exec -I $file_to_align.kangax -m 0 $qformat -R 100 -N -r 5 -X -M 5 -i $readset2 -o $sam.2 -F $sam.2.log -T $alignment_threads -w $sam.2.sqlite -W $readset2 >/dev/null"
    );
    &fix_check_kanga("$sam.2");
   }
   &process_cmd("$RealBin/util/merge_left_right_nameSorted_SAMs.pl --left_sam $sam.1 --right_sam $sam.2  -D $fragment_max_length -C 100 > $sam.t"   );
   unlink("$sam.1");
   unlink("$sam.2");
   &process_cmd("$samtools_exec view -F12 -h -T $file_to_align -o $sam $sam.t");
   unlink("$sam.t");
  }
  else {
   ## Single END:
   unless ( -s "$sam" ) {
    &process_cmd("$kanga_exec -I $file_to_align.kangax -m 0 -q 0 -R 100 -N -r 5 -X -M 5 -i $readset -o $sam -F $sam.log -T $alignment_threads -w $sam.sqlite -W $readset >/dev/null"    );
    &fix_check_kanga("$sam");
   }
  }
  confess "Alignment for $file_to_align vs $readset failed\n" unless -s $sam;
  #TODO CHECK
  &process_cmd("ln $sam $baseout.sam.namesorted");
#  link( $result_dir.fileparse($sam), $result_dir.fileparse($baseout) . ".sam.namesorted" );
 }
 else {
  confess "Kanga is currently only supporting FASTQ\n";
 }
}

sub align_bwa_index(){
   my $file_to_align = shift;
   unless (-s $file_to_align . '.bwt'){
      print "\t\t\tBuilding reference file for bwa...\r";
      &process_cmd("$bwa_exec index $file_to_align 2> /dev/null >/dev/null");
      print " Done!\n";
   }
}

sub align_bwa() {
 my ( $baseout, $file_to_align, $readset, $readset2, $bam, $sam ) = @_;
 my $readgroup = '@RG\tID:' . $readset;
 $readgroup =~ s/\.bz2$//;
 $readgroup =~ s/\.fastq$//;
 $readgroup =~ s/\.bam$//;

 if ( $read_format eq 'fastq' ) {
  my $bwa_cmd = "$bwa_exec mem -t $alignment_threads -k 35 -c 100 -E 6,6 -L 3,3 -R '$readgroup' -v 1 -o $sam $file_to_align";
  $bwa_cmd .=" $readset2 " if $readset2;
  $bwa_cmd .=" 2>/dev/null ";
  &process_cmd($bwa_cmd) unless -s $sam;  
#  &process_cmd("$bwa_exec aln -e -1 -E 6 -q 20 -L $qformat -t 2 -f $baseout.sai $file_to_align $readset 2>/dev/null"  ) unless -s "$baseout.sai";
#  if ($readset2) {   &process_cmd("$bwa_exec aln -e -1 -E 6 -q 20 -L $qformat -t 2 -f $baseout.2.sai -q 10 $file_to_align $readset2 2>/dev/null"   ) unless -s "$baseout.2.sai";  }
 }
 elsif ( $read_format eq 'bam' ) {
  if ($readset2) {
    &process_cmd("$bwa_exec aln -e -1 -E 6 -q 20 -L -t $alignment_threads -b -f $baseout.2.sai -q 10 $file_to_align $readset2 2>/dev/null"   ) unless -s "$baseout.2.sai";
    &process_cmd("$bwa_exec sampe -n 20 -N 100 -a 900 -s -r '$readgroup' $file_to_align $baseout.sai $baseout.2.sai $readset $readset2 > $sam 2>/dev/null"  ) unless -s $sam;
  }else{
    &process_cmd("$bwa_exec aln -e -1 -E 6 -q 20 -L -t $alignment_threads -b -f $baseout.sai -q 10 $file_to_align $readset 2>/dev/null"  ) unless -s "$baseout.sai";
    &process_cmd("$bwa_exec samse -n 20 -r '$readgroup' $file_to_align $baseout.sai $readset > $sam 2>/dev/null"  ) unless -s $sam;
  }
 }
 
 confess "Could not produce BWA SAI for $readset\n" unless -s $sam;
 &namesort_sam($sam);
}

sub process_alignments($) {
 my ( $alignment_bam, $readset, $readset2 ) = @_;

 # we assume that if this was done, the data has already gone to the database
 return if $alignment_bam && -s "$alignment_bam.stats";

 print &mytime() . "Post-processing alignment $readset\n";
 print LOG &mytime() . "Post-processing alignment $readset\n";

 &process_cmd("$samtools_exec index $alignment_bam")
   unless -s $alignment_bam . '.bai';
 confess "File $alignment_bam cannot be indexed"
   unless -s $alignment_bam . '.bai';
 &process_cmd("$samtools_exec idxstats $alignment_bam > $alignment_bam.stats")
   unless -s "$alignment_bam.stats";
 confess "File $alignment_bam cannot be indexed"
   unless -s $alignment_bam . '.stats';
 my @sizes              = `grep -v '^\*' $alignment_bam.stats`;
 my $reads_that_aligned = `$samtools_exec view -F260 -c $alignment_bam`; chomp($reads_that_aligned);
 my $readset_metadata = &sqlite_get_readset_metadata($readset);

 #if genomewide then we want to replace the library size with the reads that aligned
 if ($genomewide){
    &sqlite_update_readset_size($readset, $reads_that_aligned);
 }

 foreach my $s (@sizes) {
  next if $s =~ /^\*/;    # unaligned reads
                          # gene_id, length, reads matched, reads,unmatched
  my @data = split( "\t", $s );
  next unless $data[2];

  #    $reads_that_aligned += $data[2];
  my $seq_md5hash = &sqlite_get_md5( $data[0] );
  my $check = &sqlite_get_expression_statistics( $seq_md5hash, $readset );
  next if ( $check && defined( $check->{'rpkm'} ) );

  my $rpkm = sprintf(
                      "%.2f",
                      (
                        $data[2] / (
                              ( $data[1] / 1000 ) *
                                ( $readset_metadata->{'total_reads'} / 1000000 )
                        )
                      )
  );
  &sqlite_add_rpkm_expression_statistics( $seq_md5hash, $readset, $rpkm,
                                          $data[2] );
 }
 my $alignment_proportion = sprintf(
                                     "%.2f%%",
                                     (
                                       $reads_that_aligned /
                                         $readset_metadata->{'total_reads'}
                                       ) * 100
 );
 print "Reads that aligned for $readset: $reads_that_aligned ($alignment_proportion %)\n";
 print LOG "Reads that aligned for $readset: $reads_that_aligned ($alignment_proportion %)\n";
}

sub process_salmon_bias() {
 my $salmon_results = shift;
 my $readset         = shift;

 warn "process_salmon_bias $salmon_results $readset\n" if $debug;

 # final file created when all data are processed and stored.
 my $output = $salmon_results; 
 $output =~s/.salmon.results$/.depth.completed/;
 if ( -s $salmon_results && !-s $output) {
  open( SALMON, $salmon_results ) || die($!);
  my $header = <SALMON>;
  while ( my $ln = <SALMON> ) {
   chomp($ln);
   my @data = split( "\t", $ln );
   next unless $data[4];
   my $seq_md5hash = &sqlite_get_md5( $data[0]);
   # ( $seq_md5hash, $readset, $eff_counts, $tpm )
   &sqlite_add_salmon_expression_statistics( $seq_md5hash, $readset,
                                              sprintf( "%.2f", $data[4] ),
                                              sprintf( "%.2f", $data[3] ),
   );

  }
  close SALMON;
 }
}

sub perform_correct_bias() {
 my ($original_bam, $fasta_file, $readset, $namesorted_sam )    = @_;
 my $original_bam_base = fileparse($original_bam);
 my $namesorted_bam =
#   ( $namesorted_sam =~ /\.bam.namesorted$/ )
#   ? $namesorted_sam :
    $namesorted_sam . '.bam';

 my $salmon_bam_base = $result_dir . $original_bam_base . ".salmon.sampled";
 my $salmon_bam      = $salmon_bam_base . '.bam';
 my $salmon_results  = $result_dir . $original_bam_base . ".salmon.results";
 
 if (-s $salmon_results){
   &process_salmon_bias( $salmon_results, $readset );
   return $original_bam;
 }
 
 print &mytime()   . "salmon: performing Illumina bias and transcript assignment corrections\n";
 warn "Counting output of $original_bam\n" if $debug;
 my $reads_that_aligned = `$samtools_exec view -F260 -c $original_bam` if -s $original_bam;chomp($reads_that_aligned);
 my $genes_that_aligned =  `$samtools_exec view -F260 $original_bam|cut -f 3 |$sort_exec -T $sort_tmp -S $sort_memory_exec -u|wc -l` if -s $original_bam;chomp($genes_that_aligned);
 
 if ( !$reads_that_aligned ) {
  warn "No reads aligned for $readset!\n";
  return;
 }

# if we are to allow salmon to do bias-correction then we need at least 200 genes with enough reads (fragments)
 my $readset_metadata = &sqlite_get_readset_metadata($readset);
 my $readset_name     = $readset_metadata->{'alias'};
 my $salmon_dir      = $result_dir . "$readset_name.bias";

# The library size is determined by the total number of reads in the readset
# this is due to the need to support searches that are not genome wide
# however i've seen pathological genome-wide datasets where only half of the
# reads actually map to the genes... what to do then? implement -genomewide flag
 my $readset_size      = $readset_metadata->{'total_reads'};

# todo get readset distribution from data if genomewide? nah: salmon does it automatically
 my $current_salmon_exec =
   $salmon_exec . " quant --sampleOut --fldMean $readset_separation --threads $threads ";
 $current_salmon_exec .= " $extra_salmon " if $extra_salmon;

 $current_salmon_exec .= " $salmon_options_single "
   if $salmon_options_single && $single_end_readsets{$readset};

 if ($reads_that_aligned < 5000000){
	$current_salmon_exec .= " --mappingCacheMemoryLimit $readset_size --numAuxModelSamples " . int($readset_size - ($readset_size*0.1)) ." ";
 }
 if ($reads_that_aligned < 2000000){
	$current_salmon_exec .= " --numBiasSamples " . int($readset_size - ($readset_size*0.1)) ." ";
 }

 if ($genes_that_aligned < 200 || $reads_that_aligned < $salmon_min_bias ) {
   warn "Read set size is less than $salmon_min_bias reads or 200 genes, salmon will not perform bias correction\n";
   $current_salmon_exec .= ' --noLengthCorrection ' if !$is_quantseq;
  }elsif (!$is_quantseq){
	$current_salmon_exec .= ' --seqBias --posBias --gcBias ';
  }if ($is_quantseq){
 	$current_salmon_exec .= ' --noLengthCorrection --minAssignedFrags 3 ';
  }


 unless ( -s $salmon_results ) {
  mkdir($salmon_dir) unless -d $salmon_dir;
  if ($use_bwa) {
   #&process_cmd("$samtools_exec sort -T $sort_tmp -o $namesorted_bam -@ $samtools_threads -n -m $sam_sort_memory $original_bam  2>/dev/null"   ) unless -s $namesorted_bam;
   if ($contextual_alignment) {
    &process_cmd("$current_salmon_exec --output $salmon_dir --targets $fasta_file --alignments $namesorted_sam > /dev/null 2> $salmon_results.log"
    ) unless -s "$salmon_dir/quant.sf";
    sleep(3);
    confess "Salmon failed to produce output\n" unless -s "$salmon_dir/quant.sf";
    &process_cmd("$samtools_exec sort -T $sort_tmp -o $salmon_bam_base.bam -@ $samtools_threads -m $sam_sort_memory $salmon_dir/postSample.bam 2>/dev/null");
    rename( "$salmon_dir/quant.sf", $salmon_results );
   }
   else {

    # not contextual.
    open( EXPR, ">$salmon_results" ) || die($!);
    print EXPR "target_id\tlength\teff_length\teff_counts\ttpm\n";
    my $fasta_obj      = new Fasta_reader($fasta_file);
    my $timer          = new Time::Progress;
    my $fasta_count    = `grep -c '^>'  < $fasta_file`;
    chomp($fasta_count);
    $timer->attr( min => 0, max => $fasta_count );

    while ( my $seq_obj = $fasta_obj->next() ) {
     my $id          = $seq_obj->get_accession();
     my $seq_md5hash = &sqlite_get_md5($id);
     my $expression_stats =
       &sqlite_get_expression_statistics( $seq_md5hash, $readset );
     next if $expression_stats->{'salmon_tpm'};
     my $seq              = $seq_obj->get_sequence();
     my $len              = length($seq);
     my $tmp_fasta_file   = $salmon_dir . '/' . $seq_md5hash . '.fasta.tmp';
     my $tmp_sam_file     = $salmon_dir . '/' . $seq_md5hash . '.tmp';
     my $tmp_salmon_file = "$tmp_sam_file.dir/quant.sf";
     open( FASTA, ">$tmp_fasta_file" ) || die($!);
     print FASTA ">$id\n$seq\n";
     close FASTA;

     #create sam
     open( SAM, ">$tmp_sam_file" ) || die($!);
     print SAM "\@HD\tVN:1.0\tSO:unsorted\n\@SQ\tSN:$id\tLN:$len\n";
     close SAM;

     #sort for salmon
     &process_cmd("$samtools_exec view -F4 -o $tmp_sam_file $original_bam '$id'  |$sort_exec -T $sort_tmp -S $sort_memory_exec -k1  >>$tmp_sam_file  " );
     if ( -s $tmp_sam_file < 200 ) {
      warn "No alignments for $id. Skipping\n" if $debug;
      print LOG "No alignments for $id. Skipping\n";
      unlink($tmp_fasta_file) unless $debug;
      unlink($tmp_sam_file)   unless $debug;
      next;
     }
     &process_cmd("$current_salmon_exec -o $tmp_sam_file.dir $tmp_fasta_file $tmp_sam_file >/dev/null  2> $salmon_results.log"
     ) unless -s $tmp_salmon_file;
     sleep(3);
     if ( -s $tmp_salmon_file < 200 ) {
      warn "No post-salmon coverage for $id. Skipping\n"
        if $debug;
      print LOG "No post-salmon coverage for $id. Skipping\n";
      unlink($tmp_fasta_file)   unless $debug;
      unlink($tmp_sam_file)     unless $debug;
      unlink($tmp_salmon_file) unless $debug;
      next;
     }
     open( IN, $tmp_salmon_file || die($!) );
     my $h = <IN>;
     while ( my $ln = <IN> ) {
      print EXPR $ln;
     }
     close IN;
     unlink($tmp_salmon_file)        unless $debug;
     unlink($tmp_sam_file)            unless $debug;
     unlink($tmp_fasta_file)          unless $debug;
     remove_tree("$tmp_sam_file.dir") unless $debug;
    }
    undef($fasta_obj);
    close EXPR;
   }
  }
  else {
   if ($contextual_alignment) {
    if ( -s $namesorted_bam ) {
     &process_cmd("$current_salmon_exec --output $salmon_dir --targets $fasta_file --alignments $namesorted_bam >/dev/null  2> $salmon_results.log"
     ) unless -s "$salmon_dir/quant.sf";
     confess "Salmon failed to produce output\n" unless -s "$salmon_dir/quant.sf";
     sleep(3);
     &process_cmd("$samtools_exec sort -T $sort_tmp -o $salmon_bam_base.bam -@ $samtools_threads -m $sam_sort_memory $salmon_dir/postSample.bam 2>/dev/null");
    }
    elsif ( -s $namesorted_sam ) {
	#NB: this needs to include headers
     &process_cmd("$current_salmon_exec --output $salmon_dir --targets $fasta_file --alignments $namesorted_sam >/dev/null  2> $salmon_results.log"
     ) unless -s "$salmon_dir/quant.sf";
     confess "Salmon failed to produce output\n"  unless -s "$salmon_dir/quant.sf";
     sleep(3);
     my $salmon_hits_file = "$salmon_dir/postSample.bam";
     &process_cmd("$samtools_exec sort -T $sort_tmp -o $salmon_bam_base.bam -@ $samtools_threads -m $sam_sort_memory $salmon_dir/postSample.bam 2>/dev/null");
    }
    rename( "$salmon_dir/quant.sf", $salmon_results );
   }
   else {
    my $timer       = new Time::Progress;
    my $fasta_count = `grep -c '^>' < $fasta_file`;
    chomp($fasta_count);
    $timer->attr( min => 0, max => $fasta_count );
    open( EXPR, ">$salmon_results" ) || die($!);
    print EXPR "target_id\tlength\teff_length\teff_counts\ttpm\n";
#    print EXPR "bundle_id\ttarget_id\tlength\teff_length\ttot_counts\tuniq_counts\test_counts\teff_counts\tambig_distr_alpha\tambig_distr_beta\tfpkm\tfpkm_conf_low\tfpkm_conf_high\tsolvable\ttpm\n";
    my $fasta_obj      = new Fasta_reader($fasta_file);

    while ( my $seq_obj = $fasta_obj->next() ) {
     my $id          = $seq_obj->get_accession();
     my $seq_md5hash = &sqlite_get_md5($id);
     my $expression_stats =
       &sqlite_get_expression_statistics( $seq_md5hash, $readset );
     next if $expression_stats->{'salmon_tpm'};
     my $seq              = $seq_obj->get_sequence();
     my $len              = length($seq);
     my $tmp_fasta_file   = $salmon_dir . '/' . $seq_md5hash . '.fasta.tmp';
     my $tmp_sam_file     = $salmon_dir . '/' . $seq_md5hash . '.tmp';
     my $tmp_salmon_file = "$tmp_sam_file.dir/quant.sf";
     open( FASTA, ">$tmp_fasta_file" ) || die($!);
     print FASTA ">$id\n$seq\n";
     close FASTA;

     #create sam
     open( SAM, ">$tmp_sam_file" ) || die($!);
     print SAM "\@HD\tVN:1.0\tSO:unsorted\n\@SQ\tSN:$id\tLN:$len\n";
     close SAM;

     #sort for salmon
     &process_cmd("$samtools_exec view -F4 $original_bam '$id' |$sort_exec -T $sort_tmp -S $sort_memory_exec -k1  >>$tmp_sam_file "     );
     if ( !-s $tmp_sam_file || -s $tmp_sam_file < 200 ) {
      warn "No alignments for $id. Skipping\n" if $debug;
      print LOG "No alignments for $id. Skipping\n";
      unlink($tmp_fasta_file) unless $debug;
      unlink($tmp_sam_file)   unless $debug;
      next;
     }
     &process_cmd("$current_salmon_exec -o $tmp_sam_file.dir --output-align-samp $tmp_fasta_file $tmp_sam_file >/dev/null  2> $salmon_results.log"
     ) unless -s $tmp_salmon_file;
     sleep(3);
     if ( !-s $tmp_salmon_file || -s $tmp_salmon_file < 200 ) {
      warn "No post-salmon coverage for $id. Skipping\n"
        if $debug;
      print LOG "No post-salmon coverage for $id. Skipping\n";
      unlink($tmp_fasta_file)          unless $debug;
      unlink($tmp_sam_file)            unless $debug;
      unlink($tmp_salmon_file)        unless $debug;
      remove_tree("$tmp_sam_file.dir") unless $debug;
      next;
     }
     open( IN, $tmp_salmon_file ) || die($!);
     my $h = <IN>;
     while ( my $ln = <IN> ) {
      print EXPR $ln;
     }
     close IN;
     unlink($tmp_salmon_file)        unless $debug;
     unlink($tmp_sam_file)            unless $debug;
     unlink($tmp_fasta_file)          unless $debug;
     remove_tree("$tmp_sam_file.dir") unless $debug;
    }
    undef($fasta_obj);
    close EXPR;
   }
  }
  unlink($namesorted_sam) unless $debug;
  &remove_tree($salmon_dir) unless $debug;
  if ($contextual_alignment) {

   unlink( $original_bam . '.orig.bai' ) if -s $original_bam . ".bai";
   unlink( $original_bam . '.orig' ) if -s $original_bam;
   rename( $original_bam,          $original_bam . '.orig' );
   rename( $original_bam . ".bai", $original_bam . '.orig.bai' );
   #TODO CHECK
   &process_cmd("ln $salmon_bam $original_bam");
   #link( $result_dir.fileparse($salmon_bam), $result_dir.fileparse($original_bam) );
   &process_cmd( "$samtools_exec index $original_bam" );
  }
 }
 &process_salmon_bias( $salmon_results, $readset );
 return $original_bam;
}

sub perform_stats() {
 if ( !$no_graphs ) {
  print "\n" . &mytime() . "Stats n graphs: Calculating per gene stats";
  print " ($uid.raw_stats.tsv)" if $debug;
  print "\n";
  print LOG &mytime() . "Calculating per gene stats\n";
 }
 else {
  print "\n" . &mytime() . "Calculating per gene stats in $result_dir \n";
  print LOG &mytime() . "Calculating per gene stats in $result_dir \n";
 }
 open( STATS,        ">$result_dir/$uid.raw_stats.tsv" )   || die($!);
 open( STATS_RATIO,  ">$result_dir/$uid.ratio.stats.tsv" ) || die($!);
 print STATS_RATIO   "Checksum\tGene alias\tReadset1\tReadset2\tRaw RPKM\tKANGADE";
 print STATS_RATIO   "\tSalmon-corrected TPM\tsalmon-corrected counts" if $perform_bias_correction;
 print STATS_RATIO   "\n";
 print STATS "Checksum\tReadset\tGene alias\tRaw RPKM\tMean\tstd. dev\tMedian\tUpper limit\tTotal hits\tGene coverage\n";
 my $timer_counter = int(0);
 my $timer         = new Time::Progress;
 $timer->attr( min => 0, max => scalar( keys %user_alias ) );

 foreach my $md5sum ( keys %user_alias ) {
  &process_main_stats($md5sum);
  $timer_counter++;
  print $timer->report( "eta: %E min, %40b %p\r", $timer_counter )
    if ( $timer_counter % 100 == 0 );
 }
}

sub get_all_expression_data() {

 # we need this in order to enable threading for coverage graphs etc.
 foreach my $seq_md5hash ( keys %user_alias ) {
  my $seq_id   = $user_alias{$seq_md5hash}{'id'};
  my $seq_size = $user_alias{$seq_md5hash}{'length'};
  for ( my $i = 0 ; $i < scalar(@readsets) ; $i++ ) {
   next  if !$readsets[$i];
   my $readset_metadata_ref = &sqlite_get_readset_metadata( $readsets[$i] );
   my $readset_name         = $readset_metadata_ref->{'alias'};
   my $stats_ref = &sqlite_get_expression_statistics( $seq_md5hash, $readsets[$i] );
   my $total_aligned_reads =
     $stats_ref->{'total_aligned_reads'} ? $stats_ref->{'total_aligned_reads'} : int(0);
   my $depth_hash_ref =
     &sqlite_get_depth_data( $seq_md5hash, $readsets[$i], 1 );
   my $gene_length_coverage_median = $stats_ref->{'gene_length_coverage_median'};
   my $gene_length_coverage_mean = $stats_ref->{'gene_length_coverage_mean'};
   my $gene_length_coverage = $stats_ref->{'gene_length_coverage'};
   if ( !$depth_hash_ref ) {
    warn "$readset_name: No depth data for $seq_id vs $readset_name. Skipping\n" if $debug;
    print LOG "$readset_name: No depth data for $seq_id vs $readset_name. Skipping\n" if $debug;
    next;
   }
   $expression_coverage_stats_hash{$seq_md5hash}{ $readsets[$i] } = {
                                              'total_reads' => $total_aligned_reads,
                                              'gene_length_coverage' => $gene_length_coverage,
                                              'gene_length_coverage_mean' => $gene_length_coverage_mean,
                                              'gene_length_coverage_median' => $gene_length_coverage_median,
                                              'depth'       => $depth_hash_ref,
                                              'expression'  => $stats_ref
   };
  }
 }
}

sub freeze_compress() {
 my $ref    = shift;
 my $freeze = freeze($ref);
 return compress( \$freeze );
}

sub thaw_decompress() {
 my $serialized   = shift;
 my $decompressed = decompress( \$serialized );
 return thaw($decompressed);
}

sub perform_coverage_graphs() {
 return if $no_graphs;
 print &mytime()
   . "Printing coverage graphs for each gene in $result_dir/gene_coverage_plots/\n";
 print LOG &mytime()
   . "Printing coverage graphs for each gene in $result_dir/gene_coverage_plots/\n";
 my $timer_counter = int(0);
 my $timer         = new Time::Progress;
 $timer->attr( min => 0, max => scalar( keys %user_alias ) );
 my $thread_helper = new Thread_helper($R_threads);
 my @ids_to_do;

 foreach my $seq_md5hash ( keys %user_alias ) {
  my $imagefile = "$result_dir/gene_coverage_plots/$seq_md5hash" . "_gene_coverage.svg";
  next if -s $imagefile;
  $timer_counter++;
  if ( $timer_counter % 100 == 0 && @ids_to_do ) {
   print $timer->report( "eta: %E min, %40b %p\r", $timer_counter );
   my $thread = threads->create( 'make_coverage_graph_wrapper', @ids_to_do );
   $thread_helper->add_thread($thread);
   @ids_to_do = ();
   sleep(5) if ( $timer_counter % 1000 == 0 );
  }
  push( @ids_to_do, $seq_md5hash );
 }

 if (@ids_to_do) {
  print $timer->report( "eta: %E min, %40b %p\r", $timer_counter );
  my $thread = threads->create( 'make_coverage_graph_wrapper', @ids_to_do );
  $thread_helper->add_thread($thread);
  @ids_to_do = ();
 }
 $thread_helper->wait_for_all_threads_to_complete();
 print &mytime() . "Done.                                           \n";
 print LOG &mytime() . "Done.\n";
 &rename_graph_files_md52gene("$result_dir/gene_coverage_plots/");
}

sub make_coverage_graph_wrapper() {
 my @ids = @_;
 return unless @ids;
 foreach my $seq_md5hash (@ids) {
  &make_coverage_graph($seq_md5hash);
 }
}

sub make_coverage_graph($$$) {
 my $seq_md5hash     = shift;
 my $seq_id          = $user_alias{$seq_md5hash}{'id'};
 my $seq_size        = $user_alias{$seq_md5hash}{'length'};
 my $memory_hash_ref = $expression_coverage_stats_hashref->{$seq_md5hash};

 # the testws Bio::Graphics had a bug and had to
 # comment out the integer addition at
 # Bio/Graphics/Glyph/xyplot.pm ca line 315
 # my $positive = $self->pos_color;# + 1073741824;
 # my $negative = $self->neg_color;# + 1073741824;
 # GD is up to date (2.53).
 #SVG is up to date (2.59).
 #GD::SVG is up to date (0.33).
 #Bio::Graphics is up to date (2.37).

 my $graph_file =
   "$result_dir/gene_coverage_plots/$seq_md5hash" . "_gene_coverage.svg";
 return if -s $graph_file;

 # I think there is a memory leak somewhere here....!
 # build a panel with one track per readset
 my $graph;
 if ($do_png) {
  $graph = Bio::Graphics::Panel->new(
                                      -length    => $seq_size,
                                      -width     => 800,
                                      -pad_left  => 40,
                                      -pad_right => 40,
  );
 }
 else {
  $graph = Bio::Graphics::Panel->new(
                                      -length      => $seq_size,
                                      -width       => 800,
                                      -pad_left    => 40,
                                      -pad_right   => 40,
                                      -image_class => 'GD::SVG'
  );
 }
 $graph->add_track(
                    arrow =>
                      Bio::SeqFeature::Generic->new(
                                                     -start        => 1,
                                                     -end          => $seq_size,
                                                     -display_name => $seq_id
                      ),
                    -description => 1,
                    -tick        => 2,
                    -fgcolor     => 'black',
                    -double      => 1,
                    -label       => 1,
                    -bump        => 0
 );

 # find out if the graph is to be printed at all
 my $gene_has_been_printed_before;
 unless ($never_skip) {
  for (my $i=0;$i<scalar(@readsets);$i++) {
   next unless $readsets[$i];
   my $readset = $readsets[$i];
   my $ref = $memory_hash_ref->{$readset};

   next
     if (   !$ref || !$ref->{'total_reads'}
          || $ref->{'total_reads'} < $binary_min_reads
 #         || !$ref->{'gene_length_coverage_mean'}
 #         || $ref->{'gene_length_coverage_mean'} == 0
 #         || ( $process_cutoff && $ref->{'gene_length_coverage_mean'} && $ref->{'gene_length_coverage_mean'} < $process_cutoff )
      );
   $gene_has_been_printed_before++;
  }
 }

 for (my $i=0;$i<scalar(@readsets);$i++) {
  next unless $readsets[$i];
  my $readset = $readsets[$i];
  my $ref          = $memory_hash_ref->{$readset};
  my $readset_name = $library_aliases{$readset};

  next unless $ref && $ref->{'expression'};;
  
  my $expression_statistics_ref = $ref->{'expression'};
  confess "Cannot find any coverage statistics for $seq_id vs $readset !\n"
    unless $expression_statistics_ref;
    
  my $mean =
    $expression_statistics_ref->{'gene_length_coverage_mean'}
    ? int( $expression_statistics_ref->{'gene_length_coverage_mean'} )
    : int(0);
  my $median =
    $expression_statistics_ref->{'gene_length_coverage_median'}
    ? int( $expression_statistics_ref->{'gene_length_coverage_median'} )
    : int(0);
  my $rpkm =
    $expression_statistics_ref->{'rpkm'}
    ? int( $expression_statistics_ref->{'rpkm'} )
    : int(0);
  my $salmon_tpm =
      $expression_statistics_ref->{'salmon_tpm'}
    ? $expression_statistics_ref->{'salmon_tpm'}
    : int(0);
  my $effective_counts =
      $expression_statistics_ref->{'salmon_eff_counts'}
    ? $expression_statistics_ref->{'salmon_eff_counts'}
    : int(0);
  my $gene_length_coverage_sd = $expression_statistics_ref->{'gene_length_coverage_sd'};

  my $pivot =
      $expression_statistics_ref->{'gene_length_coverage_mean'}
    ? $expression_statistics_ref->{'gene_length_coverage_mean'} - ( 2 * $gene_length_coverage_sd )
    : 1
    if $gene_length_coverage_sd;
  $pivot = $expression_statistics_ref->{'gene_length_coverage_mean'}
    if $expression_statistics_ref->{'gene_length_coverage_mean'} && $pivot < 0;
  my $stats_description = " RPKM: $rpkm";
  $stats_description .= " pivot: $pivot" if $pivot;
  $stats_description .= " Eff.TPM: $salmon_tpm Eff.counts: $effective_counts "
    if $expression_statistics_ref->{'salmon_tpm'};
  my $gene_length_coverage_max = $expression_statistics_ref->{'gene_length_coverage_max'};
  

  my $readset_feature =
    Bio::SeqFeature::Generic->new(
                              -display_name => $readset_name,
                              -start        => 1,
                              -end          => $seq_size,
                              -source_tag => $readset_name . $stats_description,
                              -type       => 'expression',
    );

  if (
       $expression_statistics_ref->{'total_aligned_reads'}
       && (    $never_skip
            || $gene_has_been_printed_before
            || $median > 0
            || $expression_statistics_ref->{'total_aligned_reads'} > 50 )
    )
  {

   my $depth_hash_ref = &thaw_decompress( $ref->{'depth'} );

   my $step_size = $seq_size >= 200 ? 10 : 5;

   for ( my $pos = 0 ; $pos < $seq_size ; $pos += $step_size ) {
    last if $pos >= $seq_size;
    my ( $xl, $depth );

    #sliding window for average
    for ( my $x = $pos ; $x < $pos + $step_size ; $x++ ) {
     $depth += $depth_hash_ref->{$x} ? $depth_hash_ref->{$x} : int(0);
     $xl++;
    }
    $depth = ( $depth && $depth > 0 ) ? $depth / $xl : int(0);
    my $t =
      Bio::SeqFeature::Lite->new(
                                  -start => $pos + 1,
                                  -stop  => $pos + 1,
                                  -type  => 'depth',
                                  -score => $depth
      );
    $readset_feature->add_SeqFeature($t);
   }
   $graph->add_track(
        xyplot         => $readset_feature,                               #work!
        -description   => 1,
        -label         => 1,
        -max_score     => $gene_length_coverage_max ? $gene_length_coverage_max + ( $gene_length_coverage_max * 0.1 ) : 1,
        -min_score     => 0,
        -bgcolor       => 'white',
        -fgcolor       => 'black',
        -graph_type    => 'boxes',
        -scale         => 'left',
        -height        => 100,
        -pos_color     => 'blue',
        -neg_color     => 'lightblue',
        -bicolor_pivot => $pivot
   );
  }
  elsif (    $expression_statistics_ref->{'total_aligned_reads'}
          && $expression_statistics_ref->{'total_aligned_reads'} > 500 )
  {

   #warn Dumper $expression_statistics_ref if $debug;
   warn "Something really bad happened: maybe the statistics were not properly processed...\n";
  }
  else {
   warn "No coverage of $seq_id ($seq_md5hash) vs $readset !\n"     if $debug;
   $graph->add_track(
          xyplot       => $readset_feature,                               #work!
          -description => 1,
          -max_score   => $gene_length_coverage_max ? $gene_length_coverage_max + ( $gene_length_coverage_max * 0.1 ) : 1,
          -min_score   => 0,
          -label       => 1,
          -bgcolor     => 'white',
          -fgcolor     => 'black',
          -graph_type  => 'boxes',
          -scale       => 'left',
          -height      => 100,
          -pos_color   => 'blue',
          -neg_color   => 'lightblue',
   );
  }
 }

 open( GRAPH, ">$graph_file" ) || die($!);
 if ($do_png) {
  print GRAPH $graph->png;
 }
 else {
  print GRAPH $graph->svg;
 }
 close GRAPH;
 $graph->finished;
}

sub process_main_stats($) {

 # i don't think there is a memory leak here.
 #having both bp data and rpkm, we can post process
 my $seq_md5hash = shift;
 my $seq_id      = $user_alias{$seq_md5hash}{'id'};

 # my $seq_size         = &sqlite_get_seq_length($seq_md5hash);
 my $seq_size = $user_alias{$seq_md5hash}{'length'};

 my $readsets_covered = int(0);
 $demand_all_readsets = scalar @readsets if $demand_all_readsets;

 for ( my $i = 0 ; $i < scalar(@readsets) ; $i++ ) {
  next if !$readsets[$i];

  my $readset_metadata_ref = &sqlite_get_readset_metadata( $readsets[$i] );
  my $readset_name         = $readset_metadata_ref->{'alias'};

  my $gene_length_coverage_stat             = Statistics::Descriptive::Full->new();
  my $stats_ref = &sqlite_get_expression_statistics( $seq_md5hash, $readsets[$i] );
  my $rpkm = $stats_ref->{'rpkm'} ? $stats_ref->{'rpkm'} : int(0);
  my $total_aligned_reads = $stats_ref->{'total_aligned_reads'} ? $stats_ref->{'total_aligned_reads'} : int(0);
  my ( $no_bp_covered, @gene_length_coverages );
  my $depth_hash_ref_compressed = &sqlite_get_depth_data( $seq_md5hash, $readsets[$i], 1 );

  if ( !$depth_hash_ref_compressed ) {
   warn "No depth data for $seq_id vs $readset_name. Skipping\n"
     if $debug;
   print LOG "No depth data for $seq_id vs $readset_name. Skipping\n" if $debug;
   next;
  }
  my $depth_hash_ref = &thaw_decompress($depth_hash_ref_compressed);

  # parsing depth
  for ( my $pos = 0 ; $pos < $seq_size ; $pos++ ) {
   my $depth =
       $depth_hash_ref->{$pos}
     ? $depth_hash_ref->{$pos}
     : int(0);
   push( @gene_length_coverages, $depth );
   $no_bp_covered++ if $depth == 0;
  }
  $gene_length_coverage_stat->add_data( \@gene_length_coverages );
 
  my $gene_length_coverage_median = $gene_length_coverage_stat->median() ? $gene_length_coverage_stat->median() : int(0);
  my $gene_length_coverage_sd =    $gene_length_coverage_stat->standard_deviation()
    ? sprintf( "%.2f", $gene_length_coverage_stat->standard_deviation() )
    : int(0);
  my $gene_length_coverage_mean =    $gene_length_coverage_stat->mean() ? sprintf( "%.2f", $gene_length_coverage_stat->mean() ) : int(0);
  my $gene_length_coverage_max = $gene_length_coverage_stat->max() ? $gene_length_coverage_stat->max() : int(0);
  my $aligned_reads_per_base =    $total_aligned_reads
    ? sprintf( "%.2f", ( $total_aligned_reads / $seq_size ) )
    : int(0);
  $readsets_covered++;
  my $gene_length_coverage =    sprintf( "%.2f", ( ( $seq_size - $no_bp_covered ) / $seq_size ) );

  &sqlite_add_main_expression_statistics(
                         $seq_md5hash, $readsets[$i], $gene_length_coverage_mean, $no_bp_covered,
                         $aligned_reads_per_base,  $gene_length_coverage_median,   $gene_length_coverage_max,  $gene_length_coverage_sd, $gene_length_coverage );

  print STATS "$seq_md5hash\t$readset_name\t$seq_id\t$rpkm\t$gene_length_coverage_mean\t"
	."$gene_length_coverage_sd\t$gene_length_coverage_median\t$gene_length_coverage_max\t"
	."$total_aligned_reads\t$gene_length_coverage\n";

  if ( !$only_alignments ) {
    $stats_ref =
     &sqlite_get_expression_statistics( $seq_md5hash, $readsets[$i] );

   $expression_coverage_stats_hash{$seq_md5hash}{ $readsets[$i] } = {
                                          'gene_length_coverage' => $gene_length_coverage,
                                          'gene_length_coverage_mean' => $gene_length_coverage_mean,
                                          'gene_length_coverage_median' => $gene_length_coverage_median,
                                          'depth' => $depth_hash_ref_compressed,
                                          'expression' => $stats_ref
   };
  }
 }

 if ( $readsets_covered == 0
      || ( $demand_all_readsets && $readsets_covered < $demand_all_readsets ) )
 {
  warn "Skipping $seq_id as it didn't align to sufficient readsets ($readsets_covered)\n"
    if $debug;
  print LOG "Skipping $seq_id as it didn't align to sufficient readsets ($readsets_covered)\n";
  $skipped_references{$seq_md5hash} = 1;
  return;
 }
}

sub print_binary_table() {
 my $out              = $result_dir . "$uid.expression_levels.binary.tsv";
 my $in               = $result_dir . "$uid.raw_stats.tsv" || die($!);
 my @ordered_readsets = sort keys %library_metadata;
 print "Producing binary table $out\n";
 open( IN, $in );
 my $headers = <IN>;
 while ( my $ln = <IN> ) {
  chomp($ln);
  my (
       $seq_md5hash,     $readset_name, $seq_id,     $rpkm,
       $gene_length_coverage_mean,        $gene_length_coverage_sd,       $gene_length_coverage_median, $gene_length_coverage_max,
       $total_aligned_reads, $gene_length_coverage
  ) = split( "\t", $ln );

  # designate as not expressed:
  if (   !$gene_length_coverage || $gene_length_coverage == 0
   #    || $gene_length_coverage_mean == 0
   #    || ( $process_cutoff && $gene_length_coverage_mean < $process_cutoff )
    )
  {
   print LOG "Gene coverage of $seq_id vs $readset_name was not available or was zero ($gene_length_coverage). Skipping\n";
   warn "Gene coverage of $seq_id vs $readset_name was not available or was zero ($gene_length_coverage). Skipping\n"
     if $debug;
   next;
  }

  $binary_table{$seq_md5hash}{'exists'} = 1  if !$binary_table{$seq_md5hash}{'exists'};
  if ( $total_aligned_reads >= $binary_min_reads ) {
   if ( $gene_length_coverage >= $binary_min_coverage ) {
    $binary_table{$seq_md5hash}{$readset_name} = 1;
   }
  }
 }
 close IN;

 open( OUT, ">$out" );
 my $p = "Gene\t" . join( "\t", @ordered_readsets ) . "\n";
 my @headers = sort keys %{ $library_metadata{ (%library_metadata)[0] } };
 foreach my $metadata_head (@headers) {
  $p .= $metadata_head;
  foreach my $readset (@ordered_readsets) {
   $p .=
     $library_metadata{$readset}{$metadata_head}
     ? "\t" . $library_metadata{$readset}{$metadata_head}
     : "\tnone";
  }
  $p .= "\n";
 }
 print OUT $p;
 foreach my $seq_md5hash ( sort keys %binary_table ) {
  my $seq_id = $user_alias{$seq_md5hash}{'id'};
  my $print  = $seq_id;
  foreach my $readset_name (@ordered_readsets) {
   my $is_expressed = $binary_table{$seq_md5hash}{$readset_name} ? 1 : int(0);
   $print .= "\t" . $is_expressed;
  }
  print OUT $print . "\n";
 }
 close OUT;
}

sub process_expression_level() {

 my $expression_level_tsv = "$result_dir/$uid.expression_levels.stats.tsv";
 my $effective_expression_level_matrix =
   "$result_dir/$uid.expression_levels.matrix.effective";
 open( OUT,           ">$expression_level_tsv" )              || die($!);
 open( COUNTS_MATRIX, ">$counts_expression_level_matrix" )    || die($!);
 open( EFFECT_MATRIX, ">$effective_expression_level_matrix" ) || die($!);
 print OUT "Checksum\tGene alias\tReadset\tRaw_Counts\tRPKM\tsalmon_TPM\tsalmon_eff.counts\tKANGADE_counts\n";
 my $count_matrix_print = "transcript\t";

 for ( my $i = 0 ; $i < scalar(@readsets) ; $i++ ) {
  next  if !$readsets[$i];
  my $readset_metadata = &sqlite_get_readset_metadata( $readsets[$i] );
  $count_matrix_print .= $readset_metadata->{'alias'} . "\t";
 }
 chop($count_matrix_print);
 my $effective_matrix_print = $count_matrix_print;
 print COUNTS_MATRIX $count_matrix_print . "\n";
 print EFFECT_MATRIX $effective_matrix_print . "\n";
 foreach my $seq_id (@reference_sequence_list) {
  my $seq_md5hash = &sqlite_get_md5($seq_id);

  next if !$seq_md5hash || $skipped_references{$seq_md5hash};
  $count_matrix_print     = $seq_md5hash . "\t";
  $effective_matrix_print = $seq_md5hash . "\t";
  for ( my $i = 0 ; $i < scalar(@readsets) ; $i++ ) {
   next  if !$readsets[$i];
   my $readset_metadata_C = &sqlite_get_readset_metadata( $readsets[$i] );
   my $readset_name_C     = $readset_metadata_C->{'alias'};
   my $stats_ref_C =  &sqlite_get_expression_statistics( $seq_md5hash, $readsets[$i] );

   $stats_ref_C->{'total_aligned_reads'}   = int(0) if !$stats_ref_C->{'total_aligned_reads'};
   $stats_ref_C->{'salmon_eff_counts'} = int(0) if !$stats_ref_C->{'salmon_eff_counts'};
   $stats_ref_C->{'salmon_tpm'} = int(0) if !$stats_ref_C->{'salmon_tpm'};
   $stats_ref_C->{'rpkm'} = int(0) if !$stats_ref_C->{'rpkm'};
   $stats_ref_C->{'kangade_counts'} = int(0) if !$stats_ref_C->{'kangade_counts'};

   $count_matrix_print .= $stats_ref_C->{'salmon_eff_counts'} . "\t";
   $effective_matrix_print .=
     sprintf( "%.2f", $stats_ref_C->{'salmon_tpm'} ) . "\t";
   print OUT $seq_md5hash . "\t" 
     . $seq_id . "\t"
     . $readset_name_C . "\t"
     . $stats_ref_C->{'total_aligned_reads'} . "\t"
     . $stats_ref_C->{'rpkm'} . "\t"
     . $stats_ref_C->{'salmon_tpm'} . "\t"
     . $stats_ref_C->{'salmon_eff_counts'} . "\t"
     . $stats_ref_C->{'kangade_counts'} . "\n";

   # now parse each pair of readsets
   unless ($no_pairwise) {

    for ( my $k = $i + 1 ; $k < scalar(@readsets) ; $k++ ) {
     next if $readsets[$i] eq $readsets[$k];
     my $readset_metadata_E = &sqlite_get_readset_metadata( $readsets[$k] );
     my $readset_name_E     = $readset_metadata_E->{'alias'};
     my $stats_ref_E =
       &sqlite_get_expression_statistics( $seq_md5hash, $readsets[$k] );
     $stats_ref_E->{'total_aligned_reads'} = int(0) if !$stats_ref_E->{'total_aligned_reads'};
     $stats_ref_E->{'salmon_eff_counts'} = int(0)   if !$stats_ref_E->{'salmon_eff_counts'};
     $stats_ref_E->{'salmon_tpm'} = int(0)  if !$stats_ref_E->{'salmon_tpm'};
     $stats_ref_E->{'rpkm'} = int(0) if !$stats_ref_E->{'rpkm'};
     $stats_ref_E->{'kangade_counts'} = int(0)  if !$stats_ref_E->{'kangade_counts'};
     my $raw_rpkm_ratio =
       !$stats_ref_E->{'rpkm'} || $stats_ref_E->{'rpkm'} == 0
       ? 'Inf'
       : sprintf( "%.2f", $stats_ref_C->{'rpkm'} / $stats_ref_E->{'rpkm'} );

     $add_raw_fold_change->execute(
                                    $raw_rpkm_ratio, $seq_md5hash,
                                    $readsets[$i],   $readsets[$k]
     ) unless $contextual_alignment;
     $fold_changes{$seq_md5hash}{'raw'}{$readset_name_C}{$readset_name_E} =
       $raw_rpkm_ratio;
     my $kangade_ratio =
       $fold_changes{$seq_md5hash}{'kangade'}{$readset_name_C}{$readset_name_E}
       ? $fold_changes{$seq_md5hash}{'kangade'}{$readset_name_C}
       {$readset_name_E}
       : 'N/A';

     if ($perform_bias_correction) {

      my $salmon_tpm_ratio =
        !$stats_ref_E->{'salmon_tpm'} || $stats_ref_E->{'salmon_tpm'} == 0
        ? 'Inf'
        : sprintf( "%.2f",
                   $stats_ref_C->{'salmon_tpm'} /
                     $stats_ref_E->{'salmon_tpm'} );

      my $eff_counts_ratio =
          !$stats_ref_E->{'salmon_eff_counts'}
        || $stats_ref_E->{'salmon_eff_counts'} == 0
        ? 'Inf'
        : sprintf( "%.2f",
                   $stats_ref_C->{'salmon_eff_counts'} /
                     $stats_ref_E->{'salmon_eff_counts'} );
      $add_salmon_fold_change->execute( $eff_counts_ratio, $salmon_tpm_ratio, $seq_md5hash, $readsets[$i], $readsets[$k] )
        unless $contextual_alignment;
      print STATS_RATIO $seq_md5hash . "\t" 
        . $seq_id . "\t"
        . $readset_name_C . "\t"
        . $readset_name_E . "\t"
        . $raw_rpkm_ratio . "\t"
        . $kangade_ratio . "\t"
        . $salmon_tpm_ratio . "\t"
        . $eff_counts_ratio . "\t" . "\n";
     }
     else {
      print STATS_RATIO $seq_md5hash . "\t" 
        . $seq_id . "\t"
        . $readset_name_C . "\t"
        . $readset_name_E . "\t"
        . $raw_rpkm_ratio . "\t"
        . $kangade_ratio . "\n";
     }
    }

   }

  }
  chop($count_matrix_print);
  chop($effective_matrix_print);
  print COUNTS_MATRIX $count_matrix_print . "\n";
  print EFFECT_MATRIX $effective_matrix_print . "\n";
 }
 close OUT;
 close COUNTS_MATRIX;
 close EFFECT_MATRIX;
}

sub prepare_edgeR_graphs() {
 my $normalized_count_file = shift;

 my ( %edgeR_diff_expressed_genes, %edgeR_housekeeping_genes );
 my @edgeR_results      = glob( $edgeR_dir . "*_vs_*results.txt" );
 my @edgeR_housekeeping = glob( $edgeR_dir . "*_vs_*housekeeping.txt" );
 my $log_change_limit   = 2;
 my $diff_gene_list =
   $edgeR_dir . "diff_expressed.logM2.fdr$fdr_pval_cutoff.list";
 my $unchanged_gene_list = $edgeR_dir . "any.unchanged.p095.list";
 my $comparisons         = int(0);
 open( OUT, ">$diff_gene_list" ) || die($!);

 foreach my $edgeR_result (@edgeR_results) {
  next unless -s $edgeR_result;
  $comparisons++;
  open( IN, $edgeR_result ) || die($!);
  my $header = <IN>;
  while ( my $ln = <IN> ) {
   my @data = split( "\t", $ln );
   next
     unless $data[4]
      || $data[4] > $fdr_pval_cutoff
      || abs( $data[2] ) < 2;
   print OUT $data[0] . "\n"
     unless $edgeR_diff_expressed_genes{ $data[0] };
   $edgeR_diff_expressed_genes{ $data[0] } = $data[4];
  }
  close IN;
 }
 close OUT;
 my $housekeeping_gene_list =
   $edgeR_dir . "all.unchanged$comparisons.p095.list";
 foreach my $res (@edgeR_housekeeping) {
  next unless -s $res;
  $res =~ /(\w+)_vs_(\w+)/;
  my $readset1 = $1;
  my $readset2 = $2;
  open( IN, $res ) || die($!);
  my $header = <IN>;
  while ( my $checksum = <IN> ) {
   chomp($checksum);
   $edgeR_housekeeping_genes{$checksum}++;
   &sqlite_set_get_as_housekeeping( $checksum, 1, $readset1, $readset2 );
  }
  close IN;
 }
 if (%edgeR_housekeeping_genes) {
  open( OUT1,   ">$unchanged_gene_list" )    || die($!);
  open( OUTALL, ">$housekeeping_gene_list" ) || die($!);
  foreach my $checksum ( sort keys %edgeR_housekeeping_genes ) {
   &sqlite_set_get_as_housekeeping( $checksum, 1 )
     if $edgeR_housekeeping_genes{$checksum} == $comparisons;
   print OUT1 $checksum . "\t" . $edgeR_housekeeping_genes{$checksum} . "\n";
   print OUTALL $checksum . "\n"
     if $edgeR_housekeeping_genes{$checksum} == $comparisons;
  }
  close OUT1;
  close OUTALL;
 }
 if ( !%edgeR_diff_expressed_genes ) {
  warn
"No differentially expressed transcripts identified at cuttoffs: FDR p-value:$fdr_pval_cutoff, log2(fold-change):2\n\n";
 }

 &process_edgeR_graphs_differential( $normalized_count_file,
                                     \%edgeR_diff_expressed_genes );

 return ( \%edgeR_diff_expressed_genes, \%edgeR_housekeeping_genes );
}

sub process_housekeeping_stats() {

 # skip $skipped_references{$seq_md5hash};
 # %fold_changes

#ok what do i want here? i want to use one of the matrix files (mine or that after TMM normalization)
# and find which genes are invariable. i think that it would be best to be conservative and do the following:

}

sub process_completion() {

# N.B. Remember: R/FPKM is for looking at differences between libraries and mean/gene_length_coverage_median for looking at differences within a library.
# RPKM is appropriate for Single-End libraries and FPKM for Paired-End libraries
# Community has now rallied behind TPM which is what we're using here
 my $elapsed = $total_timer->report("%L");
 print "\nCompleted in $elapsed!\n";
 print LOG "\nCompleted in $elapsed!\n";
 print LOG "#Ended: " . &mytime() . "\n";
 close LOG;
 &sqlite_destroy();
 ## cleanup
 if ($cleanup) {
  foreach my $readset (@readsets) {
   my $baseout = $result_dir . $uid . '_' . fileparse($readset);
   unlink("$baseout.bam");
   unlink("$baseout.bam.bai");
  }
 }
 print "\n";

 if ($do_galaxy_cleanup){
	my @files = `find $result_dir -type f`;
	chomp(@files);
	foreach my $f (@files){
		next if ($f=~/\.pdf$/ || $f=~/\.tsv$/);
		unlink($f);
	}

 }

 exit(0);
}

sub perform_TMM_normalization_edgeR() {
 my $matrix_file = shift;
 my $matrix_base = fileparse($matrix_file);
 my $TMM_file    = $edgeR_dir . 'TMM_info.txt';

 print "\n" . &mytime
   . "Performing trimmed mean of fold-change normalization";
 print " using $matrix_file" if $debug;
 print "\n";

# AP: accounting for groups (column 2 [1] in target.files and TMM_info.txt files)
 my ($effective_TMM_matrix_file);
 unless ( -s $TMM_file ) {

  # just make some files that point out where the data are...
  # TODO Rewrite this using %groups_readsets
  my ( %ofhs, %readset_R_filenames );
  open( IN, $matrix_file ) || die($!);
  my $readset_header = <IN>;
  chomp($readset_header);
  my @readsets_from_matrix = split( "\t", $readset_header );
  shift @readsets_from_matrix;    # rid gene name
  foreach my $readset (@readsets_from_matrix) {
   my $readset_cf = $readset;
   $readset_cf =~ s/\W/_/g;
   my $count_data_file = $edgeR_dir . "$readset_cf.dat";
   open( my $ofh, ">$count_data_file" )
     || die( "Cannot write to $count_data_file\n" . $! );
   $readset_R_filenames{$readset} = $count_data_file;
   $ofhs{$readset}                = $ofh;
  }

  # create a tab delimited file which is gene tab count
  while ( my $ln = <IN> ) {
   chomp($ln);
   my @x = split( "\t", $ln );
   my $gene_name = shift @x;
   foreach my $readset (@readsets_from_matrix) {
    my $ofh   = $ofhs{$readset};
    my $count = shift @x;
    print $ofh "$gene_name\t$count\n";
   }
  }
  close IN;
  foreach my $ofh ( values %ofhs ) {
   close $ofh;
  }
  undef(@readsets_from_matrix);
  @readsets_from_matrix = keys %readset_R_filenames;
# estimate norm. factors; the fpkm is produced for the user but we are not going to use it for edgeR (we will not transform)
  &run_TMM( \%readset_R_filenames, $TMM_file );
 }
 $effective_TMM_matrix_file =
   &write_normalized_effective_file( $matrix_file, $TMM_file );

 return ( $effective_TMM_matrix_file . '.tpm' );
}

sub perform_edgeR_pairwise() {

# this will now estimate differential expression using the counts for each group
 my @groups = sort keys %groups_readsets;
 print &mytime
   . "Running multithreaded pairwise edgeR comparisons for "
   . scalar( keys %groups_readsets )
   . " groups...\n";
   
 foreach my $g (@groups) {
  print "Group $g readsets: "
    . join( ", ", keys %{ $groups_readsets{$g} } ) . "\n";
 }
 
 my $thread_helper = new Thread_helper($R_threads);
 for ( my $i = 0 ; $i < @groups - 1 ; $i++ ) {
  for ( my $k = $i + 1 ; $k < @groups ; $k++ ) {
   unless (
     -s $edgeR_dir . $groups[$i] . '_vs_' . $groups[$k] . ".edgeR.results.txt" )
   {
    my $thread = threads->create( 'run_edgeR', $groups[$i],$groups[$k] );
    $thread_helper->add_thread($thread);
   }
  }
 }
 $thread_helper->wait_for_all_threads_to_complete();
 print "\n";
}

sub process_edgeR_graphs_overview() {
 my $norm_matrix_file = shift;
 my $plot_dir = shift;
 my $type = shift;
 return if -f "$norm_matrix_file.per_gene_plots.completed";
 mkdir($plot_dir) unless -d $plot_dir;
 print &mytime()
   . "Producing expression graphs for each gene. This can take up to a minute per gene (for 1 CPU and 50 libraries)\n";
 print LOG &mytime()
   . "Producing expression graphs for each gene. This can take up to a minute per gene (for 1 CPU and 50 libraries)\n";
 chdir($edgeR_dir);
 $norm_matrix_file = basename($norm_matrix_file);

 my $do_png_R = $do_png ? 'TRUE' : 'FALSE';
 my $gene_plots_cmd =
"edgeR_gene_plots_all('$md5_aliases_file','$norm_matrix_file','$plot_dir',$do_png_R,$R_threads,'$type')";

 my $R_script = "$norm_matrix_file.R";
 open( OUTR, ">$R_script" ) or confess "Error, cannot write to $R_script. $!";
 print OUTR "
source('$RealBin/R/dew_funcs.R')
$gene_plots_cmd
    ";
 close OUTR;

 if ( !-s "$R_script.log" || -s "$R_script.err" > 50 ) {
  &process_cmd(
"R --no-restore --no-save --slave -f $R_script > $R_script.log 2>$R_script.err"
  );
  system("touch $norm_matrix_file.per_gene_plots.completed");
  &get_figure_legend_effective_plots($plot_dir);
  &rename_graph_files_md52gene($plot_dir);

 }

}

sub process_edgeR_graphs_differential() {

 my $normalized_count_file          = shift;
 my $edgeR_diff_expressed_genes_ref = shift;
 my $norm_matrix_file               = fileparse($normalized_count_file);
 my $plot_dir                       = $result_dir . 'gene_expression_plots/';
 chdir($edgeR_dir);
 my $diff_expr_matrix = "$norm_matrix_file.diff_expressed";
 return if -s $diff_expr_matrix;
 open( OUT, ">$diff_expr_matrix" ) || die($!);
 open( IN,  $norm_matrix_file )    || die($!);
 my $header = <IN>;
 chomp($header);
 print OUT $header . "\tFDR\n";

 while ( my $ln = <IN> ) {
  chomp($ln);
  my @data = split( "\t", $ln );
  if ( $edgeR_diff_expressed_genes_ref->{ $data[0] } ) {
   print OUT $ln . "\t" . $edgeR_diff_expressed_genes_ref->{ $data[0] } . "\n";
  }
 }
 close IN;
 close OUT;

 if ($no_js_graphs) {
  chdir($cwd);
  return;
 }

 # we now have the heatmaps and cluster data. Produce JSON objects for web
 # NEWICK data
 my $html_file =
   &prepare_heatmap_for_canvas(
                                $uid,
                                $diff_expr_matrix,
                                $diff_expr_matrix . '.data.json',
                                $diff_expr_matrix . '.samples.json',
                                $diff_expr_matrix . '.genes.json',
                                $diff_expr_matrix . '.hcsamples.newick',
                                $diff_expr_matrix . '.hcgenes.newick'
   );
 chdir($cwd);
}

sub prepare_heatmap_for_canvas() {
 my ( $title, $diff_expr_matrix, $data_file, $sample_file, $gene_file,
      $hcsamples_newick, $hcgenes_newick )
   = @_;
 return unless $hcgenes_newick && -s $hcgenes_newick;
 my (%fdr);
 open( MATRIX, $diff_expr_matrix );
 my $header = <MATRIX>;
 while ( my $ln = <MATRIX> ) {
  my @data = split( "\t", $ln );
  chomp( $data[-1] );
  $fdr{ $data[0] } = $data[-1];
 }
 close MATRIX;

#print &mytime()     . "canvas: Preparing HTML/JS code for interactive graph of $htmlfile\n";
 my $js_graph_dir = $edgeR_dir . 'js_plots/';
 my $htmlfile     = $js_graph_dir . basename($diff_expr_matrix);
 $htmlfile .= '.heatmap.html';
 my $js_dir     = $js_graph_dir . 'js/';
 my $css_dir    = $js_graph_dir . 'css/';
 my $images_dir = $js_graph_dir . 'images/';
 system("rsync --cvs-exclude -pur $RealBin/web/ $js_graph_dir/");
 my %for_json_array;

# ignore for now $for_json_array{'z'} and $for_json_array{'x'} but they must be defined
 %{ $for_json_array{'z'} } = ();
 %{ $for_json_array{'x'} } = ();
 %{ $for_json_array{'y'} } = ();
 %{ $for_json_array{'t'} } = ();
 $for_json_array{'y'}{'desc'} = 'Scaled log2(fold-change+1)';
 open( IN, $sample_file ) || die($!);
 $for_json_array{'y'}{'smps'} = decode_json(<IN>);    # array
 close IN;
 open( IN, $gene_file ) || die($!);
 $for_json_array{'y'}{'vars'} = decode_json(<IN>);    # array
 close IN;
 open( IN, $hcsamples_newick ) || die($!);
 $for_json_array{'t'}{'smps'} = <IN>;
 chomp( $for_json_array{'t'}{'smps'} );
 $for_json_array{'t'}{'smps'} =~ s/\;$//;

 #negative values not allowed!
 if ( $for_json_array{'t'}{'smps'} =~ /-\d\.\d\d+/ ) {
  warn "NEWICK of samples has negative values!\n";
 }
 close IN;
 open( IN, $hcgenes_newick ) || die($!);
 $for_json_array{'t'}{'vars'} = <IN>;
 chomp( $for_json_array{'t'}{'vars'} );
 $for_json_array{'t'}{'vars'} =~ s/\;$//;

 #negative values not allowed!
 if ( $for_json_array{'t'}{'vars'} =~ /-\d\.\d\d+/ ) {
  warn "NEWICK of genes has negative values!\n";
 }
 close IN;

 # add var data / ancillary data / properties
 open( IN, $data_file ) || die($!);
 my $data = decode_json(<IN>);
 close IN;
 foreach my $readset ( @{ $for_json_array{'y'}{'smps'} } ) {
  if ( $library_metadata{$readset} && %{ $library_metadata{$readset} } ) {
   foreach my $key ( keys %{ $library_metadata{$readset} } ) {
    push( @{ $for_json_array{'x'}{$key} }, $library_metadata{$readset}{$key} );
   }
  }
 }
 for ( my $i = 0 ; $i < @{ $for_json_array{'y'}{'vars'} } ; $i++ ) {
  my $seq_md5sum = $for_json_array{'y'}{'vars'}->[$i];

  #  foreach my $seq_md5sum ( @{ $for_json_array{'y'}{'vars'} } ) {
#  my @aliases = &sqlite_get_seq_aliases($seq_md5sum);
  my @aliases = sort keys %{$gene_aliases_threads{$seq_md5sum}};
  @aliases = qw/none/ unless @aliases;

  #debug trial
  my $alias =
    $user_alias{$seq_md5sum} ? $user_alias{$seq_md5sum}{'id'} : $seq_md5sum;
  $for_json_array{'y'}{'vars'}->[$i] = $alias;

# kind of have to do this because R has no access to the aliases and aliases are important as they are the only way to ensure that
# non-identical sequences are used...
  $for_json_array{'t'}{'vars'} =~ s/\b$seq_md5sum\b/$alias/;
  push( @{ $for_json_array{'z'}{'checksum'} }, $seq_md5sum );
  push( @{ $for_json_array{'z'}{'gene'} },     join( ', ', @aliases ) );
  push( @{ $for_json_array{'z'}{'fdr'} },      $fdr{$seq_md5sum} );
  push( @{ $for_json_array{'y'}{'data'} },     $data->{$seq_md5sum} );

 }
 open( HTML, ">$htmlfile" ) || die( "Cannot write to $htmlfile " . $! );
 print HTML '
  <!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">
  <html>
  <head>
  <title>Heatmap for ' . $title . '</title>
  <meta http-equiv="X-UA-Compatible"  content="chrome=1">
  <meta http-equiv="CACHE-CONTROL"  CONTENT="NO-CACHE">
  <meta http-equiv="Content-Type"  content="text/html; charset=utf-8"/>
  <link rel="stylesheet" href="css/canvasXpress.css" type="text/css"/>
  <!--[if lt IE 9]><script type="text/javascript" src="./js/flashcanvas.js"></script><![endif]-->
  <script type="text/javascript" src="js/canvasXpress.min.js" ></script>

  <script>
    var show_heatmap = function() {
     var cx = new CanvasXpress(
                "canvas",  ' . encode_json( \%for_json_array ) . '  , {
                    "graphType" : "Heatmap",
                    "showVarDendrogram" : true,
                    "showSmpDendrogram" : true,
                    "showSampleNames" : true,
                    "smpTitle" : "Libraries",
                    "showVariableNames" : true,
   ';
 print HTML '
                    "smpOverlays": ' . encode_json( \@sample_overlays ) . ','
   if @sample_overlays;

#   print HTML '
#                   "varOverlays": ["gene"],' if scalar( @{ $for_json_array{'y'}{'smps'} } ) == 2;
 print HTML '
                    "indicatorWidth": 5
                });
   ';
 print HTML '
   cx.kmeansVariables();
   ' if scalar( @{ $for_json_array{'y'}{'smps'} } ) > 2;
 print HTML '             
    }
</script>
</head>
<body onload="show_heatmap();">
        <h1>Heatmap of ' . $title . '</h1>
        <div>
          <canvas
            id="canvas"
            width="613"
            height="500"
          ></canvas>
        </div><br/><br/>
        ' . &return_canvas_instr() . '
</body>
</html>
  ';
 close HTML;
 return $htmlfile;
}

sub prepare_scatter_for_canvas() {
 my $diff_expr_matrix = shift;

 #2D: Y:logFC; X:logCPM; two series: for FDR <= $fdr_pval_cutoff
 #3D: Y:logFC; X:logCPM; Z:FDR
 my @edgeR_results = glob( $edgeR_dir . "*_vs_*.edgeR.all.txt" );

 # preliminary
 my ( %comparisons, %aliases );
 my $js_graph_dir = $edgeR_dir . 'js_plots/';
 my $htmlfile_2d =
   $js_graph_dir . basename($diff_expr_matrix) . '.scatterplots_2D.html';
 my $htmlfile_3d =
   $js_graph_dir . basename($diff_expr_matrix) . '.scatterplots_3D.html';
 my $js_dir     = $js_graph_dir . 'js/';
 my $css_dir    = $js_graph_dir . 'css/';
 my $images_dir = $js_graph_dir . 'images/';
 system("rsync --cvs-exclude -pur $RealBin/web/ $js_graph_dir/");

 my $print_2d = '
  <!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">
  <html>
  <head>
  <title>2D Scatterplots</title>
  <meta http-equiv="X-UA-Compatible"  content="chrome=1">
  <meta http-equiv="CACHE-CONTROL"  CONTENT="NO-CACHE">
  <meta http-equiv="Content-Type"  content="text/html; charset=utf-8"/>
  <link rel="stylesheet" href="css/canvasXpress.css" type="text/css"/>
  <!--[if lt IE 9]><script type="text/javascript" src="./js/flashcanvas.js"></script><![endif]-->
  <script type="text/javascript" src="js/canvasXpress.min.js" ></script>
  <script> var show_Scatter2D = function() {';
 my $print_3d = '
  <!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">
  <html>
  <head>
  <title>3D Scatterplots</title>
  <meta http-equiv="X-UA-Compatible"  content="chrome=1">
  <meta http-equiv="CACHE-CONTROL"  CONTENT="NO-CACHE">
  <meta http-equiv="Content-Type"  content="text/html; charset=utf-8"/>
  <link rel="stylesheet" href="css/canvasXpress.css" type="text/css"/>
  <!--[if lt IE 9]><script type="text/javascript" src="./js/flashcanvas.js"></script><![endif]-->
  <script type="text/javascript" src="js/canvasXpress.min.js" ></script>
  <script> var show_Scatter3D = function() {';
 open( HTML2D, ">$htmlfile_2d" )
   || die( "Cannot write to $htmlfile_2d " . $! );
 open( HTML3D, ">$htmlfile_3d" )
   || die( "Cannot write to $htmlfile_2d " . $! );
 print HTML2D $print_2d;
 print HTML3D $print_3d;

 # not used  %{ $for_json_array{'t'} } = ();
 foreach my $edgeR_result (@edgeR_results) {
  next unless -s $edgeR_result;
  my $res            = fileparse($edgeR_result);
  my $single_2d_file = $js_graph_dir . basename($res) . '.scatterplots_2D.html';
  my $single_3d_file = $js_graph_dir . basename($res) . '.scatterplots_3D.html';
  open( HTML1_2D, ">$single_2d_file" );
  open( HTML1_3D, ">$single_3d_file" );
  print HTML1_2D $print_2d;
  print HTML1_3D $print_3d;
  $res =~ /^(.+)_vs_(.+)\.edgeR/;
  my $sampleA    = $1;
  my $sampleB    = $2;
  my $comparison = $sampleA . ' vs ' . $sampleB;
  $comparisons{$comparison} = 1;
  my $comparison_computer = $comparison;
  $comparison_computer =~ s/[\s\.]+/_/g;
  my %for_json_array;
  %{ $for_json_array{'z'} } = ();    # gene names etc
  %{ $for_json_array{'x'} } = ();    # will probably remain empty
  %{ $for_json_array{'y'} } = ();    #logFC logCPM and qvalue
  $for_json_array{'y'}{'desc'} =
"$comparison: M / Fold change (log) vs Abundance / Counts Per Million (log)";
  @{ $for_json_array{'y'}{'smps'} } =
    ( 'logFC', 'logCPM', 'qvalue', '-log(qvalue)' );
  my $counter = int(0);
  open( IN, $edgeR_result ) || die($!);
  my @headers = split( "\t", <IN> );

  while ( my $ln = <IN> ) {
   chomp($ln);
   my @data = split( "\t", $ln );
   next unless $data[4];
   my $seq_md5sum = $data[0];
   $data[5] = log( $data[4] ) / -1;
   #my @aliases = &sqlite_get_seq_aliases($seq_md5sum);
   my @aliases = sort keys %{$gene_aliases_threads{$seq_md5sum}};
   @aliases = qw/none/ unless @aliases;
   my $alias =
     $user_alias{$seq_md5sum} ? $user_alias{$seq_md5sum}{'id'} : $seq_md5sum;
   push( @{ $for_json_array{'y'}{'vars'} },     $alias );
   push( @{ $for_json_array{'z'}{'checksum'} }, $seq_md5sum );
   push( @{ $for_json_array{'z'}{'aliases'} },  join( ", ", @aliases ) );
   push(
         @{ $for_json_array{'y'}{'data'}[$counter] },
         (
           sprintf( "%.2f", $data[1] ), sprintf( "%.2f", $data[2] ),
           $data[4], $data[5]
         )
   );
   $counter++;
  }
  close IN;
  my $printx_2d = '
  var '
    . $comparison_computer
    . '_2d = new CanvasXpress("'
    . $comparison_computer
    . '_2d", '
    . encode_json( \%for_json_array ) . '
  ,{
                  "graphType": "Scatter2D",
                  "yAxis": ["logFC"],
                  "xAxis": ["logCPM"],
                  "setMaxX": 20,
                  "setMinX": 0,
                  "setMaxY": 20,
                  "setMinY": -20,
                  "yAxisTitle": "logFC (M)",
                  "xAxisTitle": "logCPM (A)",
                  "sizeBy": "-log(qvalue)",
                  "colorBy": "qvalue",
                  "indicatorWidth": 3,
                  "legendPosition": "bottom",
                  "indicatorsPosition": "right"
                }
      );
   ';
  my $printx_3d = ' 
  var '
    . $comparison_computer
    . '_3d = new CanvasXpress("'
    . $comparison_computer
    . '_3d", '
    . encode_json( \%for_json_array ) . '
  ,{
                  "graphType": "Scatter3D",
                  "yAxis": ["logFC"],
                  "xAxis": ["logCPM"],
                  "zAxis":["-log(qvalue)"],
                  "setMaxX": 20,
                  "setMinX": 0,
                  "setMaxY": 20,
                  "setMinY": -20,
                  "yAxisTitle": "logFC (M)",
                  "xAxisTitle": "logCPM (A)",
                  "zAxisTitle": "FDR (-log)",
                  "sizeBy": "-log(qvalue)",
                  "colorBy": "qvalue",
                  "indicatorWidth": 3,
                  "legendPosition": "bottom",
                  "indicatorsPosition": "right",
                  "scatterType": false
                }
      );
   ';
  print HTML2D $printx_2d;
  print HTML3D $printx_3d;
  print HTML1_2D $printx_2d;
  print HTML1_3D $printx_3d;
  print HTML1_2D '}</script></head><body onload="show_Scatter2D();">';
  print HTML1_3D '}</script></head><body onload="show_Scatter3D();">';
  print HTML1_2D '
        <div style="clear: both;">
        <br /><br />
        <h1>2D Scatterplot of ' . $comparison . '</h1>
          <canvas
            id="' . $comparison_computer . '_2d"
            width="1200"
            height="800"
          ></canvas>
        </div>';
  print HTML1_3D '
        <div style="clear: both;">
        <br /><br />
        <h1>3D Scatterplot of ' . $comparison . '</h1>
          <canvas
            id="' . $comparison_computer . '_3d"
            width="1200"
            height="800"
          ></canvas>
        </div>';
  print HTML1_2D '<br/><br/>' . &return_canvas_instr() . '</body></html>';
  print HTML1_3D '<br/><br/>' . &return_canvas_instr() . '</body></html>';
  close HTML1_2D;
  close HTML1_3D;
 }

 # 2d graph
 print HTML2D '}</script></head><body onload="show_Scatter2D();">';
 print HTML3D '}</script></head><body onload="show_Scatter3D();">';
 foreach my $comparison ( sort keys %comparisons ) {
  my $comparison_computer = $comparison;
  $comparison_computer =~ s/[\s\.]+/_/g;
  print HTML2D '
        <div style="clear: both;">
        <br /><br />
        <h1>2D Scatterplot of ' . $comparison . '</h1>
        
          <canvas
            id="' . $comparison_computer . '_2d"
            width="1200"
            height="800"
          ></canvas>
        </div>';
  print HTML3D '
        <div style="clear: both;">
        <br /><br />
        <h1>3D Scatterplot of ' . $comparison . '</h1>
        
          <canvas
            id="' . $comparison_computer . '_3d"
            width="1200"
            height="800"
          ></canvas>
        </div>';
 }
 print HTML2D '<br/><br/>' . &return_canvas_instr() . '</body></html>';
 print HTML3D '<br/><br/>' . &return_canvas_instr() . '</body></html>';
 close HTML2D;
 close HTML3D;
 return ( $htmlfile_2d, $htmlfile_3d );
}

sub return_canvas_instr() {
 return "
  <div style='clear: both;'><br/><br/>
  <h2>Graph instructions</h2>
  <h3>Zooming and Panning</h3>
<p>
Zooming graphs can be done by dragging the mouse over an area while pressing the left mouse button or by clicking close to the axis after the axis resizer appears. You can adjust the maximum and minimum value and / or drag an interval across the range. You can select samples by dragging the mouse over the samples you want to see while you pressing the left mouse button and the 'shift' key. Also, networks and haetmaps can be zoom in and out with the help of the mouse wheel or with the help. Panning can be done also in networks and heatmaps by either dragging the mouse or with the help of the the arrow keys.
</p>

<h3>Selecting data points or nodes</h3>
<p>
You can select data points in the scatter plots or in the networks by simultanously pressing the 'shift' key and dragging the mouse over the data points or nodes. You can also select individual items by simultaneously pressing the 'ctrl' key and clicking with the left mouse button on the individual item. Once selected, you can press simultaneously the 'ctrl' and the 'delete' keys to hide them or the 'ctrl' and the 'insert' keys to show them again. You can reset the selecting by pressing the 'esc' key.
</p>

<h3>Animations</h3>
<p>You can use the arrow keys to rotate the 3D-Scatter plot or cycle over the axes in the 2D-Scatter plots (including the paging keys too). Of course you have to simultaneously press any of those keys and the 'ctrl' key.</p>

<h3>Resizing, reseting and printing the canvas</h3>
<p>You can resize the canvas image using the handle that appears when you mouse over the south and east sides of the canvas.</p>
In order to reset the canvas just press the 'esc' key.
To print the canvas you need to simultaneously click the 'ctrl' and the 'p' keys.</p>
</div>  ";
}

sub prepare_for_dge() {
 my $seq_md5sum        = shift;
 my $dispersion_method = shift;

 # get the table of counts, either via salmon or from survey.

#deseq script from Robles/Alexie
#  my $rscript = '
#    library(DESeq);
#    args <- commandArgs(TRUE);
#    InFile  <- args[1];
#    OuFile  <- args[2];
#    Nreps   <- args[3];
#    RsltsPath  <- args[4];
#    countsTable <- read.delim(InFile, header=TRUE, stringsAsFactors=TRUE );
#    rownames( countsTable ) <- countsTable$X;
#    countsTable <- countsTable[ , -1 ];
#    #Identify ALL conditions:
#    AllConds <- c(rep("ref",Nreps),rep("exp",Nreps));
#    cat("Performing DESeq DE Calculation for n =",length(AllConds)/2,"replicates...\n");
#    cds <- newCountDataSet( countsTable, AllConds );
#    cds <- estimateSizeFactors( cds );
#    cds <- estimateDispersions( cds ); #method is blind without reps...
#    res <- nbinomTest( cds, "ref","exp" );
#    # Write output to a csv file (All isoforms):
#    write.table(res[ order(res$pval), ] ,file=sprintf("%s%s",OuFile,"_DESeq.csv"),sep=",");
#    ';
}

sub run_TMM {
 #from b.haas
 # second column is group and taken from $library_metadata{$name}{'group' }
 # my $readset_name         = $readset_metadata_ref->{'alias'};
 my ( $readset_R_filenames_href, $output ) = @_;
 my $data_file_list = $edgeR_dir . "all_data_files.list";
 open( OUT1, ">$data_file_list" ) || die($!);
 print OUT1 "files\tgroup\tdescription\n";

 foreach my $name ( keys %{$readset_R_filenames_href} ) {
  my $R_filename = $readset_R_filenames_href->{$name};
  print OUT1 join( "\t", $R_filename, $library_metadata{$name}{'group'}, $name )
    . "\n";
 }
 close OUT1;
 my $tmm_norm_script = $edgeR_dir . "tmp_runTMM.R";
 open( OUT2, ">$tmm_norm_script" ) || confess "$!";

 # TODO using minimum CPM/replicates
 print OUT2 "
  source('$RealBin/R/dew_funcs.R')
  TMMresult<-TMM_normalize('$data_file_list',$minCPM,$minLibs,'$output')
    ";
 close OUT2;
 &process_cmd(
"R --no-restore --no-save --slave -f $tmm_norm_script >$tmm_norm_script.log 2>$tmm_norm_script.err",
  $edgeR_dir
 );
 unless ( -s $output ) {
  print LOG "TMM normalization failed.\n";
  confess "TMM normalization failed.\n";
 }
 unlink($tmm_norm_script) unless $debug;
}

sub run_edgeR {

#init. from b.haas
# AP: this is running pair-wise for each dataset
# %groups_readsets has key1 as group and then key2 as readset name and value is data filename
# which is tab-delimited: gene	count
# TODO   m1<-glm(c(c1a,c1b,c1c,c2a,c2b,c2c)~as.factor(c(s1,s1,s1,s2,s2,s2)),family=quasipoisson);
#        m2<-update(m1,~.-as.factor(c(s1,s1,s1,s2,s2,s2)));
#        test<-anova(m1,m2,test="F");
#        local_stat<-test["Pr(>F)"][2,1];
 my ( $group_A, $group_B ) = @_;
 my $base        = $edgeR_dir . $group_A . '_vs_' . $group_B . ".edgeR";
 print &mytime." Group $group_A vs: $group_B as $base\n";
 my $R_script    = $base . "_run.R";
 my $target_file = $base . "_target.files";

 open( OUT1, ">$target_file" )
   or confess "Error, cannot write to $target_file: $!";
 print OUT1 "files\tgroup\tdescription\n";
 foreach my $group ( $group_A, $group_B ) {
  foreach my $readset ( keys %{ $groups_readsets{$group} } ) {
   my $file = $groups_readsets{$group}{$readset}
     || confess "Cannot find file: Group $group; readset $readset\n";
   print OUT1 join( "\t", $file, $group, $readset ) . "\n";
  }
 }
 close OUT1;
 open( OUTR, ">$R_script" ) or confess "Error, cannot write $R_script: $!";
 print OUTR "
  source('$RealBin/R/dew_funcs.R')
  edgeR_DE_explore('$md5_aliases_file','$target_file','$base',$edgeR_dispersion,$fdr_pval_cutoff,$tree_clusters,$minCPM,$minLibs);
  ";
 close OUTR;
 &process_cmd(
"R --no-restore --no-save --slave -f $R_script >$R_script.log 2>$R_script.err",
  $edgeR_dir
 );
 my $html_file =
   &prepare_heatmap_for_canvas(
                                $uid . " of $group_A vs $group_B",
                                $base . '.results.txt',
                                $base . '.data.json',
                                $base . '.samples.json',
                                $base . '.genes.json',
                                $base . '.hcsamples.newick',
                                $base . '.hcgenes.newick'
   );
 unlink($R_script)    unless $debug;
 unlink($target_file) unless $debug;
 print "\n";
}

sub write_normalized_effective_file {

#init. from b.haas
# this will take the matrix file, from e.g. salmon or raw, and convert to FPKM by
# applying the TMM normalization previously done with edgeR.
 my ( $matrix_file, $tmm_info_file ) = @_;
 my $normalized_effective_file =
   $edgeR_dir . fileparse($matrix_file) . ".normalized.effective";
 my ( %eff_lib_sizes, %TMM_hash, %norm_factors, %tpm_readsets );

 # get effective library sizes
 open( IN, $tmm_info_file )
   or confess "Error, cannot open file $tmm_info_file: $!";
 my $header = <IN>;

 while ( my $ln = <IN> ) {
  chomp($ln);
  my @data = split( "\t", $ln );
  next unless $data[2] && $data[5];
  $eff_lib_sizes{ $data[2] } = $data[5];
  $norm_factors{ $data[2] }  = $data[4];    # not used
 }
 close IN; 
 open( OUT2, ">$normalized_effective_file.tpm" ) || die($!);
 open( IN, $matrix_file ) || die($!);
 $header = <IN>;
 print OUT2 $header;
 chomp $header;
  my @readsets_from_matrix = split( "\t", $header ); 
   # first (0) key is just '#transcript'

 # get total transcript abundance first
 while ( my $ln = <IN> ) {
  chomp($ln);
  my @data    = split( "\t", $ln );
  my $gene    = $data[0];
  my $seq_len = $user_alias{$gene}{'length'}
    or confess "Error, no seq length for $gene";
  for ( my $i = 1 ; $i < scalar(@data) ; $i++ ) {
   my $readset_name      = $readsets_from_matrix[$i];
   my $readlength_median = $library_metadata{$readset_name}{'readlength_median'}
     || confess("Cannot find median readlength for $readset_name");
   my $eff_lib_size = $eff_lib_sizes{$readset_name}
     or confess "Error, no eff lib size for $readset_name";
   my $frag_count           = $data[$i];
   my $transcript_abundance = ( $frag_count * $readlength_median ) / $seq_len;
   $tpm_readsets{$readset_name}{'gene_number'}++;
   $tpm_readsets{$readset_name}{'sum_tpm'} += $transcript_abundance;
  }
 }
 close IN;
# print LOG "Average transcript abundances for each readset:\n";
# foreach my $readset ( keys %tpm_readsets ) {
#  $tpm_readsets{$readset}{'average_tpm'} =
#    sprintf( "%.4f",
#             $tpm_readsets{$readset}{'sum_tpm'} /
#               $tpm_readsets{$readset}{'gene_number'} );
#  print LOG "$readset: Gene number:"
#    . $tpm_readsets{$readset}{'gene_number'}
#    . "Average:"
#    . $tpm_readsets{$readset}{'average_tpm'} . " SUM:"
#    . $tpm_readsets{$readset}{'sum_tpm'} . " \n";
# }

 # and reopen
 open( IN, $matrix_file ) || die($!);
 $header = <IN>;
 chomp $header;

 while ( my $ln = <IN> ) {
  chomp($ln);
  my @data = split( "\t", $ln );
  my $gene = $data[0];

  my $seq_len = $user_alias{$gene}{'length'}
    or confess "Error, no seq length for $gene";
  my $print1 = $gene;
  my $print2 = $gene;
  for ( my $i = 1 ; $i < scalar(@data) ; $i++ ) {
   my $readset_name      = $readsets_from_matrix[$i];
   my $readlength_median = $library_metadata{$readset_name}{'readlength_median'}
     || confess("Cannot find median readlength for $readset_name");
   my $eff_lib_size = $eff_lib_sizes{$readset_name}
     or confess "Error, no eff lib size for $readset_name";
   my $norm_factor     = $norm_factors{$readset_name};
   my $frag_count      = $data[$i];
   my $norm_frag_count = sprintf( "%.2f",$norm_factor * $frag_count);

#   my $fpkm = $frag_count / ( $seq_len / 1e3 ) / ( $eff_lib_size / 1e6 );
#   $fpkm = sprintf( "%.2f", $fpkm );

   my $transcript_abundance = ( $frag_count * $readlength_median ) / $seq_len;
   my $tpm                  = ( $transcript_abundance * 1e6 ) / $tpm_readsets{$readset_name}{'sum_tpm'};
   $tpm = sprintf( "%.2f",$tpm);
   
   $print2 .= "\t$tpm";
   $TMM_hash{$gene}{$readset_name}{'norm_frag_count'} = $norm_frag_count;
   $TMM_hash{$gene}{$readset_name}{'TPM'}             = $tpm;
  }
  print OUT2 $print2 . "\n";
 }
 close OUT2;
 unless ( $normalized_effective_file . '.tpm'
          && -s $normalized_effective_file . '.tpm' )
 {
  warn
    "EdgeR failed to produce the output file $normalized_effective_file.tpm'\n";
  return;
 }

 # append FPKM values in the main text outfile and also do ratios
 print "Writing out results as tab delimited files (*stats.tsv)...\n";
 open( STATSIN, "$result_dir/$uid.expression_levels.stats.tsv" )
   || confess($!);
 my $header_S = <STATSIN>;
 my @headers = split( "\t", $header_S );
 if ( scalar(@headers) < 10 ) {

  # has not already been processed
  open( STATSOUT, ">$result_dir/$uid.expression_levels.stats.tsv.t" )
    || die($!);
  chomp($header_S);
  print STATSOUT $header_S . "\tTMM_Normalized.count\tTMM.TPM\n";
  while ( my $ln = <STATSIN> ) {
   chomp($ln);
   my @data = split( "\t", $ln );
   next unless $data[2];
   print STATSOUT $ln . "\t"
     . $TMM_hash{ $data[0] }{ $data[2] }{'norm_frag_count'} . "\t"
     . $TMM_hash{ $data[0] }{ $data[2] }{'TPM'} . "\n";
  }
  close STATSIN;
  close STATSOUT;
  unlink("$result_dir/$uid.expression_levels.stats.tsv");
  rename( "$result_dir/$uid.expression_levels.stats.tsv.t",
          "$result_dir/$uid.expression_levels.stats.tsv" );
 }

 my $R_script = "$result_dir/$uid.expression_levels.stats.tsv.R";
 unless ( -s $R_script ) {
  open( OUTR, ">$R_script" ) or confess "Error, cannot write to $R_script. $!";
  print OUTR "
source('$RealBin/R/dew_funcs.R')
print_statistics_normalized('$result_dir/$uid.expression_levels.stats.tsv')
    ";
  close OUTR;

  if ( !-s "$R_script.log" || -s "$R_script.err" > 50 ) {
   &process_cmd(
"R --no-restore --no-save --slave -f $R_script > $R_script.log 2>$R_script.err"
   );
  }
 }

 open( STATSIN, "$result_dir/$uid.ratio.stats.tsv" ) || die($!);
 my $headerR = <STATSIN>;
 @headers = split( "\t", $headerR );
 if ( scalar(@headers) < 11 ) {
  open( STATSOUT, ">$result_dir/$uid.ratio.stats.tsv.t" ) || die($!);
  chomp($headerR);
  print STATSOUT $headerR . "\tTMM_Normalized.count\tTMM.TPM\n";
  while ( my $ln = <STATSIN> ) {
   chomp($ln);
   my @data = split( "\t", $ln );
   next unless $data[3];

   my $norm_count_ratio =
     (    $TMM_hash{ $data[0] }{ $data[2] }{'norm_frag_count'} > 0
       && $TMM_hash{ $data[0] }{ $data[3] }{'norm_frag_count'} > 0 )
     ? sprintf(
                "%.2f",
                (
                  $TMM_hash{ $data[0] }{ $data[2] }{'norm_frag_count'} /
                    $TMM_hash{ $data[0] }{ $data[3] }{'norm_frag_count'}
                )
     )
     : int(0);

   my $tpm_ratio =
     (    $TMM_hash{ $data[0] }{ $data[2] }{'TPM'} > 0
       && $TMM_hash{ $data[0] }{ $data[3] }{'TPM'} > 0 )
     ? sprintf(
                "%.2f",
                (
                  $TMM_hash{ $data[0] }{ $data[2] }{'TPM'} /
                    $TMM_hash{ $data[0] }{ $data[3] }{'TPM'}
                )
     )
     : int(0);
   print STATSOUT $ln . "\t$norm_count_ratio\t$tpm_ratio\n";
  }
  close STATSIN;
  close STATSOUT;

  unlink("$result_dir/$uid.ratio.stats.tsv");
  rename( "$result_dir/$uid.ratio.stats.tsv.t",
          "$result_dir/$uid.ratio.stats.tsv" );
 }
 return $normalized_effective_file;
}

sub check_fastq_format() {
 my $fq        = shift;
 my $max_seqs  = 100;
 my $max_lines = $max_seqs * 4;
 my ( @lines, $number, $counter );
 if ( $fq =~ /.bz2$/ ) {
  @lines = `$bunzip2_exec -dkc $fq|head -n $max_lines`;
 }
 else {
  @lines = `head -n $max_lines $fq`;
 }
 chomp(@lines);
 for ( my $k = 0 ; $k < @lines ; $k += 4 ) {
  my $ids = $lines[$k];
  if ($ids =~ /^>/){
	return 'fasta';
  }
  confess "$fq: Not in illumina format!\n" unless $ids =~ /^@/;
  $counter++;
  my $seq   = $lines[ $k + 1 ];
  my $idq   = $lines[ $k + 2 ];
  my $qual  = $lines[ $k + 3 ];
  my @quals = split( //, $qual );
  for ( my $i = 0 ; $i <= $#quals ; $i++ ) {
   $number = ord( $quals[$i] );
   if ( $number > 75 ) {
    warn "File $fq is solexa/illumina phred64 format!\n";
    return 'phred64';
   }
  }
  last if $counter >= $max_seqs;
 }
 return 'phred33';
}

sub rename_graph_files_md52gene() {
 my $cwd2  = getcwd;
 my $dir   = shift;
 my $file_type = $dir=~/coverage/ ? 'coverage' : 'expression';
 my @files = glob( $dir . "/*" );
 mkdir( $dir . '/gene_names' ) unless -d $dir . '/gene_names';
 foreach my $file (@files) {
  next unless $file && -s $file;
  next if -d $file;
  my $filename = fileparse($file);
  my $dirname  = dirname($file);
  $filename =~ /^([^_\-\.]+)(.+)$/;
  my $md5 = $1 || next;
  my $xtn = $2 ? $2 : '';
  next unless $user_alias{$md5};
  my $new_name = $user_alias{$md5}{'id'};

  #  $new_name =~ s/[^\w\.\-\_]+/_/g;
  my $new_filename = $dirname . '/gene_names/' . $new_name . $xtn;
  symlink( "../" . $filename, $new_filename );
 }
 @files = ();
 if ( !$no_pdf && $convert_imagemagick_exec && -x $convert_imagemagick_exec ) {
  print &mytime() . "Creating PDF for $dir using imagemagick\n";
  print LOG &mytime() . "Creating PDF for $dir using imagemagick\n";
  chdir( $dir . "/gene_names" );
  my $outdir = "PDF";
  mkdir $outdir unless -d $outdir;
  my ( @file_slice, $count_files );
  my $slice_count = 1;

  my @files = glob("*png");
  push( @files, glob("*pdf") );
  push( @files, glob("*svg") );

  my $thread_helper = new Thread_helper($R_threads);

  foreach my $file ( sort @files ) {
   next unless -s $file;
   $count_files++;
   if ( $count_files && $count_files > 500 && @file_slice ) {
    my $outpdffile = $outdir."/rnaseq_".$file_type."_$slice_count.pdf";
    next if -s $outpdffile;
    my $thread = threads->create(
                                  'process_cmd',
                                  "$convert_imagemagick_exec "
                                    . join( " ", @file_slice )
                                    . " $outpdffile"
    );
    $thread_helper->add_thread($thread);
    $slice_count++;
    @file_slice  = ();
    $count_files = int(0);
   }
   push( @file_slice, "'$file'" );
  }
  my $outpdffile = $outdir."/rnaseq_".$file_type."_$slice_count.pdf";
  if ( @file_slice && !-s $outpdffile ) {

   # last ones
   my $thread = threads->create(
                                 'process_cmd',
                                 "$convert_imagemagick_exec "
                                   . join( " ", @file_slice )
                                   . " $outpdffile"
   );
   $thread_helper->add_thread($thread);
  }

  $thread_helper->wait_for_all_threads_to_complete();
  print &mytime() . "Done.\n";
  chdir($cwd2);
 }
}

sub seq_cleanup ($) {

 #IUPAC:
 # T>G>C>A
 #S: C or G -> G
 #W: A or T -> T
 #R: A or G purine -> G
 #Y: C or T/U pyrimidine -> T
 #M: A or C amino -> C
 #K: G or T/U keto -> T
 #B: C/G/T not A -> T
 #D: A/G/T not C -> T
 #H: A/C/T not G -> T
 #V: A/G/C not T -> G
 #N: A/G/C/T: any -> N
 #X -> N
 my $string = shift;
 return if !$string;
 $string =~ s/\s+//g;
 $string =~ s/\W+//g;
 return if !$string || $string =~ /^\s*$/;
 $string =~ s/S/G/g;
 $string =~ s/R/G/g;
 $string =~ s/M/C/g;
 $string =~ s/V/G/g;
 $string =~ s/X/N/g;
 $string =~ s/U/T/g;
 $string =~ s/Y/T/g;
 $string =~ s/K/T/g;
 $string =~ s/W/T/g;
 $string =~ s/B/T/g;
 $string =~ s/D/T/g;
 $string =~ s/H/T/g;
 return uc($string);
}

sub check_program() {
 my @paths;
 foreach my $prog (@_) {
  my $path = `which $prog`;
  confess "Error, path to required $prog cannot be found\n"
    unless $path =~ /^\//;
  chomp($path);
  $path = readlink($path) if -l $path;
  chomp($path);
  push( @paths, $path );
 }
 return @paths;
}

sub check_program_optional() {
 my @paths;
 foreach my $prog (@_) {
  my $path = `which $prog`;
  warn "Warning: path to optional $prog cannot be found\n"
    unless $path =~ /^\//;
  chomp($path);
  $path = readlink($path) if -l $path;
  chomp($path);
  push( @paths, $path );
 }
 return @paths;
}

sub wrap_text() {
 my $string      = shift;
 my $wrap_length = shift;
 $wrap_length = 120 if !$wrap_length;
 $string =~ s/(.{0,$wrap_length})/$1\n/g;
 $string =~ s/\n{2,}/\n/;
 return $string;
}

sub remove_redundant_sequences() {
 my $file = shift || die;
 my ( %hash, %redundant );
 my $count = int(0);
 open( OUT, ">$file.noredundant.fsa" );
 my $in_obj = new Fasta_reader($file);
 while ( my $object = $in_obj->next() ) {
  my $seq = $object->get_sequence();
  my $id  = $object->get_accession();
  my $md5 = md5_hex($seq);
  if ( !$hash{$md5} ) {
   print OUT ">$id\n" . $seq . "\n";
   $hash{$md5}++;
  }
  else {
   $hash{$md5}++;
   $redundant{$seq} = $hash{$md5};
   $count++;
  }
 }
 close OUT;
 open( OUT, ">$file.redundant.seqs" );
 foreach my $seq ( keys %redundant ) {
  print OUT $redundant{$seq} . "\t$seq\n";
 }
 close OUT;
 print "Excluded $count identical sequences from $file\n";
 return "$file.noredundant.fsa";
}

sub get_uid_time() {
 my $string = shift;
 $string = 'unique' if !$string;
 return md5_hex( time() . $string );
}

sub get_figure_legend_effective_plots() {
 my $file = shift;
 $file =~ s/\/$//;
 $file =~ s/.pdf$//;
 $file .= "_legend.txt";

 my $legend = 'Differential expression as estimated by DEW. For each gene in the assembly an expression value was derived after normalizing effective counts (from salmon) for library size and the Trimmed Mean of M-values (TMM) and converting it to a TPM value. The boxplots in the background of each gene plot show the distribution of expression values for each library: The 1st & 3rd quartile and median of these expression value for the entire library. The bars of each boxplot (whiskers) extend up to 1.5 the length of the box (up to the max/min outlier point) while the grey points denote outliers. Then for each gene, a red dot denotes the actual expression value for each library. If the red dot intercepts the Y axis (i.e. near 0), then the gene deemed as not expressed.';

 open( OUT, ">$file" ) || die;
 print OUT &wrap_text($legend);
 close OUT;

}

