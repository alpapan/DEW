#!/usr/bin/perl -w

=pod

=head1 NAME

 dew

=head1 DESCRIPTION

 It accepts a reference gene FASTA and aligns readsets against it using single-end mode of BWA.
 Since alignments take the longest, it sacrifices speed for low memory usage (uses sqlite) so it is good weak for web-server with a few GB of memory or very large datasets (TODO: farm out alignments)
 User can also provide sequences (as FASTA) to be used as a reference from the cmdline
 Makes a PNG image of the coverage of a reference gene(s) for each readset provided
 
 Files are recalculated only if non-existent/older
 if sequence exists in database with alignments, then don't re-align. use md5 on sequence 
 created server to run these things on
 Added kangade
 added fpkm via express
 added bowtie2; made it the default; accepts fastq.bz2
 when not in contextual alignment, express processes each reference sequence separately (yes, it is incorrect but db results would be incorrect otherwise and are getting the fpkm...)
 added edgeR
 added html5 interactive heatmap and scatterplot
 housekeeping genes are found automatically using pvalue 1 +/-0.05 
 Currently, highest memory use is from samtools sort (2Gb) and the memory used by depth calculations (9gb for fungal dataset). 
  latter only starts exploding when at the graph-making step (i.e. after depth has been calculated and stored) so could be improved by finding the memory leak
 SQLite is fast because it resides in memory and then written out in file; this prevents parallel runs.

 Sort requires up to 20Gb of RAM

 samtools i use is v. 0.1.18

=head1 EXAMPLE CMD

See /demo and try:

time dew.pl -dbname dew_webserver.db -library_name lib_alias.txt  \
 -contextual  -correct_bias  -uid 8f10caff61e2022d056758037df77003 \
 -1read Sp.ds.1M.left.fq \
 -2read Sp.ds.1M.right.fq \
 -1read Sp.hs.1M.left.fq \
 -2read Sp.hs.1M.right.fq \
 -1read Sp.log.1M.left.fq \
 -2read Sp.log.1M.right.fq \
 -i 8f10caff61e2022d056758037df77003.query -o 8f10caff61e2022d056758037df77003.dew.output -thre 4
 
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

            'infile:s'              => Reference file of genes
            'sequence:s'            => Instead of a file, provide a single sequence (in FASTA format with \n as line separator);
            'format:s'              => The reads can be in BAM or FASTQ format. FASTQ can be .bz2 if bowtie2 is used
            '1read|readset1|r:s{1,}'=> Sets of files (one per library). Tested with Phred33 FASTQ format
            '2read|readset2:s{1,}'  => If provided, do paired end alignment. Sets of paired 'right' files (synced to readset1). Optional.
            'samtools_exec:s'       => Executable to samtools if not in your path
            'bwa_exec:s'            => Executable to BWA if not in your path
            'bowtie2_exec:s'        => Executable to Bowtie2 if not in your path
            'bamtools_exec:s'       => Executable to bamtools if not in your path
            'uid:s'                 => A uid for naming output files. Optional, otherwise generate
            'threads:i'             => Number of CPUs to use for alignment. BWA has no advantage over 4 threads
            'library_name_file:s'   => An tag value tab delimited file (filename/alias) for giving a friendly alias for each readset library. Needs a header line to describe columns. Only include -1read files.
            'median_cutoff:i'       => Median number of hits across reference must be above cutoff
            'need_all_readsets'     => All sets of reads must have alignments against the gene in order for it to be processed. Otherwise, 1+ is sufficient. 
            'over'                  => Allow overwriting of any files with the same name
            'nographs'              => Do not produce any graphs. Graphs can take a very long time when there are many readsets (e.g. 30+ libraries and 30k+ genes). Also there is a memory leak somewhere...
            'gene_graphs_only'      =>  The opposite of above; only do enough work to get the gene depth/coverage graphs and then exit
            'contextual'            => Complete realignment of all genes in order to run a correction of biases properly. Does not read/store data in the cache
            'use_bwa'               => Use BWA instead of Bowtie2
            'correct_bias'          => Use eXpress to correct Illumina sequencing biases and transcript isofrm assignments. Increases runtime. Use -contextual for accuracy 
	        'prepare_only'	    => Quit after post-processing readsets and writing initial DB
	        'seeddb:s'              => Initialize database using this database file (e.g. after -prepare_only)
	        'kanga'                 => Experimental:  use kanga
	        'existing_aln:s{1,}'    => Experimental: use existing bam (read name sorted)
	        'resume'		    => Load existing data from database and do not reprocess existing readsets (good for adding new readsets even with contextual. NB assumes same FASTA input so DO NOT use if you changed the FASTA reference gene file)
	        'no_kangade|nokangade'  => Do not use kangade to process pairwise libraries
	        'db_use_file'	    => Use file for SQL database rather than system memory (much slower but possible to analyze larger datasets)
	        'dispersion'		=> For edgeR: if we have replicates, dispersion can be set to auto. otherwise set it to a float such as 0.1 (def)
	        'fdr_cutof'			=> Cut off FDR for DE as float (0.001 def.)

=head1 AUTHORS

 Alexie Papanicolaou

        1 CSIRO Ecosystem Sciences, Canberra, Australia
        alexie@butterflybase.org
        
 Many thanks/credit to B. Haas (The Broad) for code bundled in Trinity-RNA-Seq

=head1 DISCLAIMER & LICENSE

This software is released under the GNU General Public License version 3 (GPLv3).
It is provided "as is" without warranty of any kind.
You can find the terms and conditions at http://www.opensource.org/licenses/gpl-3.0.html.
Please note that incorporating the whole software or parts of its code in proprietary software
is prohibited under the current license.

=head1 BUGS & LIMITATIONS

Making of BioPerl graphs has a memory leak somewhere...

See TODO, otherwise none other known so far but probably lots


=head1 TODO

high: for me 
med:  for bioinformatics alpha tester
low:  for end-user beta tester


& AP; med; think of a way to speed up gene graphs. Eg. fork it out while kangade is running? For 35libs and 30k genes, we're talking about &stats_with_graph taking 1000 min vs 100 min w/o graphs
 
* AP; med; think what to report with the housekeeping genes. 

* AP;low: visualize SNP variation using JBrowse; genes4all?


* AP: low/time consuming; farm out alignments to another computer or farm (without shared filesystem)
* AP: low/time consuming; Gets a list of housekeeping genes based on SD. Plots distribution to get 10% most robust ones. User input also allowed.
      'housekeeping:s{,}'     => A file or a list with IDs from the reference genes that will act as housekeeping (no expression change between readsets)

* AP;low: add canvas graphs for coverage?

* Temi;low: Report Wizard

* AP; low: can use parallel gnu sort but it is only available in newer machines.

=head1 NOT TODO

* AP; tried; Support kanga. Kanga -q1 was tried but issues: Minor (perf): it was slower than bwa. Minor (perf): we had to sort the sam file for both name- and coord-sorting as the
  unmapped reads are placed before mapped ones. Minor (perf): it does not support sam/bam input, like bwa does (neither does bowtie2). 
  Major: kangar are very large so no benefit but does not accept any compressed format (bz2, bam etc);  Major: huge amounts of memory needed

=cut

use strict;
use Pod::Usage;
use Getopt::Long;
use Digest::MD5 qw(md5_hex);
use Statistics::Descriptive;
use Time::Progress;
use Time::localtime;
use File::Basename;
use File::stat;
use File::Path qw(remove_tree);
use File::Copy;
use Text::CSV_XS;
use Bio::SeqIO;
use Cwd;
use JSON;
use IO::Compress::Bzip2 qw(bzip2 $Bzip2Error);
use IO::Uncompress::Bunzip2 qw(bunzip2 $Bunzip2Error);
use FindBin qw/$RealBin/;
use lib $RealBin;

#db
use DBI qw(:sql_types);
use Storable qw(freeze thaw);

#graphics
use Bio::Graphics;
use Bio::SeqFeatureI;
use Bio::SeqFeature::Generic;
use Bio::SeqFeature::Lite;
$| = 1;

#debug
use Data::Dumper;
my $debug = 0;
$ENV{'PATH'} = $ENV{'PATH'} . ':' . $RealBin;
#################################
my (
	$input_reference_file,  @readsets,
	@housekeeping_ids,      %housekeeping,
	$reference_file_md5sum, @reference_sequence_list,
	%library_aliases,       $contextual_alignment,
	%library_metadata,      $extra_genes,
	%readset_lookup,        $perform_bias_correction,
	%fold_changes,          %skipped_references,
	%user_alias,            %groups_readsets
);
my $db_hostname   = 'localhost';
my $samtools_exec = `which samtools`;
chomp($samtools_exec);
my $kangade_exec = `which kangade`;
chomp($kangade_exec);
my $bowtie2_exec = `which bowtie2`;
chomp($bowtie2_exec);
my $bwa_exec = `which bwa`;
chomp($bwa_exec);
my $express_exec = `which express`;
chomp($express_exec);
my $bamtools_exec = `which bedtools`;
chomp($bamtools_exec);
my $ps2pdf_exec = `which gs`;
chomp($ps2pdf_exec);
my $inkscape_exec = `which inkscape`;
chomp($inkscape_exec);
my $pdfcrop_exec = `which pdfcrop`;
chomp($pdfcrop_exec);
my $kanga_exec = `which kanga`;
chomp($kanga_exec);
my $kangax_exec = `which kangax`;
chomp($kangax_exec);

if ($ps2pdf_exec) {
	$ps2pdf_exec .=
	  " -sPAPERSIZE=a0 -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -q -sOutputFile=";
}
my $edgeR_dispersion       = 0.1;
my $fdr_pval_cutoff        = 0.001;
my $tree_clusters          = 10;
my $read_format            = 'fastq';
my $threads                = 3;
my $min_housekeeping_genes = 5;
my $median_cutoff          = 1;
my (
	$uid,               $lib_alias_file,     $demand_all_readsets,
	$use_kanga,         @use_existing_bam,   $overwrite_results,
	@readsets2,         $doing_paired_end,   $db_file,
	$db_use_file,       $dbname,             $no_graphs,
	$user_ref_sequence, $no_kangade,         $main_output_dir,
	$use_bwa,           $prepare_input_only, @sample_overlays,
	$initial_db,        $resume,             $gene_graphs_only
);
GetOptions(
	'infile:s'      => \$input_reference_file,
	'extra_genes:s' => \$extra_genes,
	'sequence:s'    => \$user_ref_sequence,
	'format:s'      => \$read_format,
	'1read|readset1|r:s{1,}' => \@readsets,         # if bowtie2, it can be /bz2
	'2read|readset2:s{1,}'   => \@readsets2,        # if bowtie2, it can be /bz2
	'kangade_exec:s'         => \$kangade_exec,
	'samtools_exec:s'        => \$samtools_exec,
	'bowtie2_exec:s'         => \$bowtie2_exec,
	'bamtools_exec:s'        => \$bamtools_exec,    # SE bowtie2 only
	'bwa_exec:s'             => \$bwa_exec,         # modified bwa
	'uid:s'                  => \$uid,
	'threads:i'              => \$threads,
	'library_name_file:s'    => \$lib_alias_file,
	'prepare_only' => \$prepare_input_only,

	#  'housekeeping:s{,}'      => \@housekeeping_ids,
	'median_cutoff:i'      => \$median_cutoff,
	'need_all_readsets'    => \$demand_all_readsets,
	'over'                 => \$overwrite_results,
	'nographs'             => \$no_graphs,
	'hostname:s'           => \$db_hostname,
	'dbname:s'             => \$dbname,
	'seeddb:s'             => \$initial_db,
	'o:s'                  => \$main_output_dir,
	'contextual'           => \$contextual_alignment,
	'correct_bias'         => \$perform_bias_correction,
	'bwa'                  => \$use_bwa,
	'kanga'                => \$use_kanga,
	'existing_aln:s{1,}'   => \@use_existing_bam,
	'debug:i'              => \$debug,
	'resume'               => \$resume,
	'no_kangade|nokangade' => \$no_kangade,
	'db_use_file'          => \$db_use_file,
	'gene_graphs_only'     => \$gene_graphs_only,
	'dispersion:s'         => \$edgeR_dispersion,          #auto for bio.reps
	'fdr_cutoff:f'         => \$fdr_pval_cutoff
);
my $bunzip2_exec = `which pbzip2`;
chomp($bunzip2_exec);
my $bunzip_threads = $threads <= 6 ? $threads : 6;
$bunzip2_exec .= " -p$bunzip_threads " if $bunzip2_exec;
$bunzip2_exec = `which bzip2` if !$bunzip2_exec;
chomp($bunzip2_exec);
die "Cannot find bzip2\n" unless $bunzip2_exec;
$uid    = md5_hex( time() . 'dew' )         unless $uid;
$dbname = $uid . '_transcriptome.sqlite.db' unless $dbname;
my $cwd = getcwd;
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
	$add_express_fold_change,
	$add_raw_fold_change,
	$get_kangade_fold_change,
	$get_express_fold_change,
	$get_raw_fold_change,
	$update_express_expression_statistics,
	$update_kangade_expression_statistics,
	$update_rpkm_expression_statistics,
	$get_readset_filename
) = &sqlite_init();
my $file_for_alignment = &prepare_input_data();

if ( !-s $counts_expression_level_matrix
	|| ( $contextual_alignment && $debug && $debug >= 2 ) )
{
	&starts_alignments($file_for_alignment);
	&perform_stats();
	&sqlite_backup() unless $db_use_file;
	exit if $gene_graphs_only;
	&process_expression_level();
	&sqlite_backup() unless $db_use_file;
	close STATS;
	close STATS_RATIO;
}
my $check_lines = `wc -l < $counts_expression_level_matrix`;
die "No expression data available"
  unless ( -s $counts_expression_level_matrix && $check_lines > 1 );
my $fpkm_expression_level_matrix_TMM =
  &perform_TMM_normalization_edgeR($counts_expression_level_matrix);
$check_lines = `wc -l < $fpkm_expression_level_matrix_TMM `
  if $fpkm_expression_level_matrix_TMM;
die "No differential expression found"
  unless ( -s $fpkm_expression_level_matrix_TMM && $check_lines > 1 );

&perform_edgeR_pairwise();
&prepare_edgeR_graphs();
my ( $html2d, $html3d ) =
  &prepare_scatter_for_canvas($fpkm_expression_level_matrix_TMM);
&process_completion();
exit();
#########################################################################################################
sub process_cmd {
	my ( $cmd, $dir ) = @_;
	print "CMD: $cmd\n" if $debug;
	chdir($dir) if $dir;
	my $ret = system($cmd);
	if ( $ret && $ret != 256 ) {
		die "Error, cmd died with ret $ret\n";
	}
	chdir($cwd) if $dir;
	return;
}

sub mytime() {
	my $hour =
	  localtime->hour() < 10 ? '0' . localtime->hour() : localtime->hour();
	my $min = localtime->min() < 10 ? '0' . localtime->min() : localtime->min();
	my $sec = localtime->sec() < 10 ? '0' . localtime->sec() : localtime->sec();
	return "$hour:$min:$sec\t";
}

sub sqlite_backup() {
    return if $debug;
	my $keep = shift;
	my $backup_db_file = $db_use_file ? $db_file . '.backup' : $db_file;
	print "\nCheckpointing database. Do not halt...\r";

	# backup file in case it crashes later.
	$dbh->do("VACUUM");

	# ensure no problems if it crashes while backing up
	$dbh->sqlite_backup_to_file( $backup_db_file . ".tmp" );
	if ( -s $backup_db_file . ".tmp" ) {
		unlink( $backup_db_file . '.old' );
		rename( $backup_db_file, $backup_db_file . '.old' ) if $keep;
		unlink($backup_db_file);
		rename( $backup_db_file . ".tmp", $backup_db_file );
		print
		  "Checkpointing database. Do not halt... SQL database checkpointed!\n";
	}
	else {
		print
"Checkpointing database. Do not halt... Checkpointing failed or there were no data to write out!\n";
	}
}

sub sqlite_init() {
	$db_file = $debug ? $cwd . "/$dbname.debug" : $cwd . "/$dbname";
	my $active_db_file = $db_use_file ? $db_file . '.active' : ':memory:';
	if ( $initial_db && $initial_db ne $db_file ) {
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
	die "Cannot create SQLite database" unless $db_version;
	print "\tUsing SQLite DB $db_version\n";
	if ($db_existed) {
		$dbh->sqlite_backup_from_file($db_file);
		### PRAGMAS for speed
		$dbh->do("PRAGMA journal_mode = MEMORY");
		$dbh->do("PRAGMA temp_store = 2 ");          # memory
		$dbh->do("PRAGMA cache_size = -448000 ");    # 400mb of RAM for sqlite
		$dbh->do("PRAGMA quick_check");
		$dbh->do("PRAGMA synchronous = OFF") if $db_use_file;
	}
	else {
		print "\tCreating database...\n";
		$dbh->do("PRAGMA encoding = 'UTF-8'");

#$dbh->do("CREATE TABLE file_alignments(file_md5sum char(32),readset_id integer,bam blob)"    );
#$dbh->do("CREATE UNIQUE INDEX file_alignments_idx1 ON file_alignments(file_md5sum,readset_id)"    );
		$dbh->do(
"CREATE TABLE sequence_data(seq_md5hash char(32) primary key,seq_length integer,housekeeping integer DEFAULT 0)"
		);
		$dbh->do(
			"CREATE TABLE sequence_aliases (seq_md5hash char(32), alias text)");
		$dbh->do(
			"CREATE INDEX sequence_aliases_idx ON sequence_aliases(seq_md5hash)"
		);
		$dbh->do(
"CREATE TABLE readsets (readset_id INTEGER PRIMARY KEY,readset_file varchar(255),total_reads integer, is_paired boolean, alias varchar(255), ctime timestamp)"
		);
		$dbh->do("CREATE UNIQUE INDEX readsets_idx1 ON readsets(readset_file)");
		$dbh->do(
"CREATE TABLE expression_statistics (seq_md5hash char(32), readset_id integer, mean_hits REAL, no_coverage integer, rpkm integer, mean_reads REAL, median_hits integer, total_hits integer, max_hits integer, hit_sd REAL, express_fpkm integer, express_eff_counts REAL, kangade_counts integer)"
		);
		$dbh->do(
"CREATE UNIQUE INDEX expression_statistics_idx1 ON expression_statistics(seq_md5hash,readset_id)"
		);

		#tmp for uncached
		$dbh->do(
"CREATE TABLE expression_statistics_tmp (seq_md5hash char(32), readset_id integer, mean_hits REAL, no_coverage integer, rpkm integer, mean_reads REAL, median_hits integer, total_hits integer, max_hits integer, hit_sd REAL, express_fpkm integer, express_eff_counts REAL, kangade_counts integer)"
		);
		$dbh->do(
"CREATE UNIQUE INDEX expression_statistics_tmp_idx1 ON expression_statistics_tmp(seq_md5hash,readset_id)"
		);

		# r-tree?
		$dbh->do(
"CREATE TABLE depth (seq_md5hash char(32), readset_id integer, data blob)"
		);
		$dbh->do("CREATE INDEX depth_idx1 ON depth(seq_md5hash,readset_id)");
		$dbh->do(
"CREATE TABLE depth_tmp (seq_md5hash char(32), readset_id integer, data blob)"
		);
		$dbh->do(
			"CREATE INDEX depth_tmp_idx1 ON depth_tmp(seq_md5hash,readset_id)");
		$dbh->do(
"CREATE TABLE kangade_analysis (seq_md5hash char(32), readset1_id INTEGER, readset2_id INTEGER,Classification INTEGER,Score INTEGER,DECntsScore INTEGER,PearsonScore INTEGER,CtrlUniqueLoci INTEGER,"
			  . "ExprUniqueLoci INTEGER,CtrlExprLociRatio INTEGER,PValueMedian REAL,PValueLow95 REAL,PValueHi95 REAL,TotCtrlCnts INTEGER,TotExprCnts INTEGER,TotCtrlExprCnts INTEGER,ObsFoldChange REAL,FoldMedian REAL,"
			  . "FoldLow95 REAL,FoldHi95 REAL,ObsPearson REAL,PearsonMedian REAL,PearsonLow95 REAL,PearsonHi95 REAL)"
		);
		$dbh->do(
"CREATE UNIQUE INDEX kangade_analysis_idx1 ON kangade_analysis(seq_md5hash,readset1_id,readset2_id)"
		);
		$dbh->do(
"CREATE TABLE fold_changes (seq_md5hash char(32), readset1_id INTEGER, readset2_id INTEGER, raw_rpkm REAL, express_fpkm REAL, express_effective_counts REAL, kangade_observed REAL,housekeeping integer DEFAULT 0)"
		);
		$dbh->do(
"CREATE UNIQUE INDEX fold_changes_idx1 ON fold_changes(seq_md5hash,readset1_id,readset2_id)"
		);
	}
	### PRAGMAS for speed
	$dbh->do("PRAGMA journal_mode = MEMORY");
	$dbh->do("PRAGMA temp_store = 2 ");          # memory
	$dbh->do("PRAGMA cache_size = -448000 ");    # 400mb of RAM for sqlite
	$dbh->do("PRAGMA synchronous = OFF") if $db_use_file;
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
	my $add_express_fold_change = $dbh->prepare(
"UPDATE fold_changes SET express_fpkm=?,express_effective_counts=? WHERE seq_md5hash=? AND readset1_id=(SELECT readset_id from readsets where readset_file=?) AND readset2_id=(SELECT readset_id from readsets where readset_file=?)"
	);
	my $add_raw_fold_change = $dbh->prepare(
"UPDATE fold_changes SET raw_rpkm=? WHERE seq_md5hash=? AND readset1_id=(SELECT readset_id from readsets where readset_file=?) AND readset2_id=(SELECT readset_id from readsets where readset_file=?)"
	);
	my $get_kangade_fold_change = $dbh->prepare(
"SELECT kangade_observed FROM fold_changes WHERE seq_md5hash=? AND readset1_id=(SELECT readset_id from readsets where readset_file=?) AND readset2_id=(SELECT readset_id from readsets where readset_file=?)"
	);
	my $get_express_fold_change = $dbh->prepare(
"SELECT express_fpkm, express_effective_counts FROM fold_changes WHERE seq_md5hash=? AND readset1_id=(SELECT readset_id from readsets where readset_file=?) AND readset2_id=(SELECT readset_id from readsets where readset_file=?)"
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
	my $check_hash_from_seqdb = $dbh->prepare(
		"SELECT seq_md5hash FROM sequence_data WHERE seq_md5hash=?");
	my $add_seqhash_to_seqdb = $dbh->prepare(
		"INSERT INTO sequence_data (seq_md5hash,seq_length) VALUES (?,?)");
	my $depth_table = 'depth';
	if ($contextual_alignment) { $depth_table .= '_tmp'; }
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
	if ($contextual_alignment) { $expression_statistics_table .= '_tmp'; }
	my $init_expression_statistics = $dbh->prepare(
"INSERT INTO $expression_statistics_table (seq_md5hash, readset_id) VALUES (?,(SELECT readset_id from readsets where readset_file=?)) "
	);
	my $update_rpkm_expression_statistics = $dbh->prepare(
"UPDATE $expression_statistics_table set rpkm=?, total_hits=? WHERE seq_md5hash=? AND readset_id=(SELECT readset_id from readsets where readset_file=?)"
	);
	my $update_expression_statistics = $dbh->prepare(
"UPDATE $expression_statistics_table set mean_hits=?, no_coverage=?, mean_reads=?, median_hits=?, max_hits=?, hit_sd=? WHERE seq_md5hash=? AND readset_id=(SELECT readset_id from readsets where readset_file=?)"
	);
	my $update_express_expression_statistics = $dbh->prepare(
"UPDATE $expression_statistics_table set express_fpkm=?,express_eff_counts=? WHERE seq_md5hash=? AND readset_id=(SELECT readset_id from readsets where readset_file=?)"
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
	my $set_housekeeping = $dbh->prepare(
		"UPDATE sequence_data set housekeeping=? WHERE seq_md5hash=?");
	my $set_housekeeping_fold = $dbh->prepare(
"UPDATE fold_changes set housekeeping=? WHERE seq_md5hash=? AND readset1_id=? AND readset2_id=? "
	);
	my $get_readset = $dbh->prepare(
"SELECT readset_id,is_paired,alias,total_reads,ctime from readsets WHERE readset_file=?"
	);
	my $get_readset_filename =
	  $dbh->prepare("SELECT readset_file from readsets WHERE alias=?");
	my $add_readset = $dbh->prepare(
"INSERT INTO readsets (readset_file,is_paired,alias,total_reads,ctime) VALUES (?,?,?,?,?)"
	);
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
		$add_express_fold_change,
		$add_raw_fold_change,
		$get_kangade_fold_change,
		$get_express_fold_change,
		$get_raw_fold_change,
		$update_express_expression_statistics,
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
		$set_housekeeping_fold->execute( 1, $seq_md5hash, $readset1,
			$readset2 );
	}
}

sub sqlite_destroy() {
	print &mytime . "Closing SQL connections\n";
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
	$update_express_expression_statistics->finish();
	$update_kangade_expression_statistics->finish();
	$add_readset->finish();
	$get_sequence_alias->finish();
	$add_sequence_alias->finish();
	$add_to_kangade_analysis->finish();
	$check_kangade_analysis->finish();
	$delete_kangade_analysis->finish();
	$start_fold_change->finish();
	$check_fold_change->finish();
	$add_kangade_fold_change->finish();
	$add_express_fold_change->finish();
	$add_raw_fold_change->finish();
	$get_kangade_fold_change->finish();
	$get_express_fold_change->finish();
	$update_rpkm_expression_statistics->finish();
	$get_raw_fold_change->finish();

	if ( $contextual_alignment && !$debug ) {

		#empty temporary tables and re-create their schema
		$dbh->do("DROP TABLE expression_statistics_tmp");
		$dbh->do(
"CREATE TABLE expression_statistics_tmp (seq_md5hash char(32), readset_id integer, mean_hits REAL, no_coverage integer, rpkm integer, mean_reads REAL, median_hits integer, total_hits integer, max_hits integer, hit_sd REAL, express_fpkm integer, express_eff_counts REAL, kangade_counts integer)"
		);
		$dbh->do(
"CREATE UNIQUE INDEX expression_statistics_tmp_idx1 ON expression_statistics_tmp(seq_md5hash,readset_id)"
		);
		$dbh->do("DROP TABLE depth_tmp");
		$dbh->do(
"CREATE TABLE depth_tmp (seq_md5hash char(32), readset_id integer, data blob)"
		);
		$dbh->do(
			"CREATE INDEX depth_tmp_idx1 ON depth_tmp(seq_md5hash,readset_id)");
	}
	&sqlite_backup() unless $db_use_file;
	$dbh->disconnect();
	undef($dbh);
	unlink( $db_file . '.active' ) if $db_use_file;
}

=pod

how to store the entire alignment in the database.... not used.

sub sqlite_check_align_old($) {
  my $readset = shift;
  $check_db->execute( $reference_file_md5sum, $readset );
  my $result = $check_db->fetchrow_arrayref();
  $result = $result ? int(1) : int(0);
  print "Readset $readset was already in DB. Not re-aligning.\n" if $result;
  return $result;
}

sub sqlite_add_align_old($) {
  my $readset    = shift;
  my $bamfile    = shift;
  my $bamcontent = `cat $bamfile`;
  $add_to_db->bind_param( 1, $reference_file_md5sum );
  $add_to_db->bind_param( 2, $readset );
  $add_to_db->bind_param( 3, $bamcontent, SQL_BLOB );
  $add_to_db->execute();
  $check_db->execute( $reference_file_md5sum, $readset );
  my $result = $check_db->fetchrow_arrayref();
  $result = $result ? int(1) : int(0);
  print "Readset $readset added to DB\n" if $result;
  return $result;
}

sub sqlite_get_align_old($) {
  my $readset = shift;
  my $bam     = shift;
  $get_from_aligndb->execute( $reference_file_md5sum, $readset );
  my $row = $get_from_aligndb->fetchrow_arrayref();
  my $result = $row ? int(1) : int(0);
  open( BAM, ">$bam" );
  print BAM $row->[0];
  close BAM;
  &process_cmd("$samtools_exec index $bam")
    unless -s $bam . '.bai' && ( -s $bam . '.bai' ) > 200;
  return $result;
}

=cut

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
	die "Cannot find data for $id\n" unless $row;
	return $row->[0];
}

sub sqlite_check_md5($) {
	my $id = shift;
	$get_hash_from_seqdb->execute($id);
	my $row = $get_hash_from_seqdb->fetchrow_arrayref();
	return $row->[0] if $row && $row->[0];
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
}

sub sqlite_get_readset_metadata($) {
	my $readset_filename = shift;
	$get_readset->execute($readset_filename);
	my $result = $get_readset->fetchrow_hashref();
	return $result;
}

sub sqlite_add_readset_metadata($$$) {
	my $readset_filename = shift;
	my $lib_alias        = shift;
	my $total_reads      = shift;
	my $is_paired        = $doing_paired_end ? 1 : int(0);
	my $result           = &sqlite_get_readset_metadata($readset_filename);
	$add_readset->execute( $readset_filename, $is_paired, $lib_alias,
		$total_reads, localtime() )
	  if !$result;
}

sub sqlite_check_expression_statistics($$) {
	my ( $seq_md5hash, $readset ) = @_;
	$check_expression_statistics->execute( $seq_md5hash, $readset );
	my $result = $check_expression_statistics->fetchrow_arrayref();
	return $result && $result->[0] ? 1 : undef;
}

sub sqlite_init_expression_statistics($) {
	my ( $seq_md5hash, $readset ) = @_;
	unless (
		defined(
			&sqlite_check_expression_statistics( $seq_md5hash, $readset )
		)
	  )
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
		warn Dumper($check) if $debug;
		die
"Could not add express expression statistics for: $seq_md5hash,$readset,$rpkm,$total_reads_hit\n";
	}
}

sub sqlite_add_main_expression_statistics($) {
	my (
		$seq_md5hash, $readset,     $mean_hits, $no_coverage,
		$mean_reads,  $median_hits, $max_hits,  $hit_sd
	) = @_;
	$update_expression_statistics->execute(
		$mean_hits, $no_coverage, $mean_reads,  $median_hits,
		$max_hits,  $hit_sd,      $seq_md5hash, $readset
	);
}

sub sqlite_add_express_expression_statistics() {
	my ( $seq_md5hash, $readset, $fpkm, $eff_counts ) = @_;
	my $r =
	  $update_express_expression_statistics->execute( $fpkm, $eff_counts,
		$seq_md5hash, $readset );
	my $check = &sqlite_get_expression_statistics( $seq_md5hash, $readset );
	if ( !$check || !defined( $check->{'express_fpkm'} ) ) {
		warn Dumper $check if $debug;
		die
"Could not add express expression statistics for: $seq_md5hash,$readset,$fpkm,$eff_counts\n";
	}
}

sub sqlite_add_kangade_expression_statistics() {
	my ( $seq_md5hash, $readset, $kangade_counts ) = @_;
	$update_kangade_expression_statistics->execute( int($kangade_counts),
		$seq_md5hash, $readset );
	my $check = &sqlite_get_expression_statistics( $seq_md5hash, $readset );
	if ( !$check || !defined( $check->{'kangade_counts'} ) ) {
		warn Dumper $check if $debug;
		die
"Could not add kangade expression statistics for: $seq_md5hash,$readset,$kangade_counts\n";
	}
}

sub sqlite_get_expression_statistics($$) {

	#my ( $seq_md5hash, $readset, $do_alias ) = @_;
	#if ($do_alias) {
	#  $get_readset_filename->execute($readset);
	#  my $ref = $get_readset_filename->fetchrow_arrayref();
	#  $readset = $ref->[0] if $ref;
	#}
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
	my $seq_md5hash      = shift;
	my $readset_filename = shift;
	my $hash_ref         = shift;
	my $freeze           = freeze($hash_ref);
	my $base_depth_data_serialized;
	bzip2 \$freeze => \$base_depth_data_serialized
	  || die "bzip2 failed: $Bzip2Error\n";
	$check_depth_data->execute( $seq_md5hash, $readset_filename );
	my $result = $check_depth_data->fetchrow_arrayref();
	$delete_depth_data->execute( $seq_md5hash, $readset_filename ) if $result;
	$add_depth_data->bind_param( 1, $seq_md5hash );
	$add_depth_data->bind_param( 2, $readset_filename );
	$add_depth_data->bind_param( 3, $base_depth_data_serialized, SQL_BLOB );
	$add_depth_data->execute();

  #  warn "Adding depth data for $seq_md5hash vs $readset_filename\n" if $debug;
}

sub sqlite_get_depth_data($$) {
	my $seq_md5hash      = shift;
	my $readset_filename = shift;
	$check_depth_data->execute( $seq_md5hash, $readset_filename );
	my $result = $check_depth_data->fetchrow_arrayref();
	if ( $result->[0] ) {
		my $ref;
		bunzip2 \$result->[0] => \$ref || die "bunzip2 failed: $Bunzip2Error\n";
		return thaw($ref);
	}
	return;
}

sub sqlite_start_fold_change() {
	my ( $seq_md5hash, $readset1, $readset2 ) = @_;
	my $check =
	  $check_fold_change->execute( $seq_md5hash, $readset1, $readset2 );
	if ( !$check ) {
		$start_fold_change->execute( $seq_md5hash, $readset1, $readset2 );
	}
	$check = $check_fold_change->execute( $seq_md5hash, $readset1, $readset2 );
	die "Failed to populate fold change database\n" unless $check;
}

sub perform_checks_preliminary() {
	print &mytime . "Performing preliminary checks\n";
	my $sqlite3_exec = `which sqlite3`;
	pod2usage "SQLite 3 not found\n" unless $sqlite3_exec;
	pod2usage "BWA not found\n"
	  unless $bwa_exec && -s $bwa_exec && -x $bwa_exec;
	pod2usage "Samtools not found\n"
	  unless $samtools_exec && -s $samtools_exec && -x $samtools_exec;
	$read_format = lc($read_format);
	pod2usage "Read format can only be BAM or FASTQ\n"
	  unless $read_format eq 'bam' || $read_format eq 'fastq';
	pod2usage "No input\n"
	  unless ( ( $input_reference_file && -s $input_reference_file )
		|| ( $user_ref_sequence && length($user_ref_sequence) > 100 ) );
	pod2usage "Insufficient readsets\n"
	  unless @readsets && scalar(@readsets) > 1;

	if (@readsets2) {
		if ( scalar(@readsets2) != scalar(@readsets) ) {
			pod2usage
"Number of files for second readset must be 0 or equal to number of 1st readset\n";
		}
		$doing_paired_end = 1;
	}
	for ( my $i = 0 ; $i < @readsets ; $i++ ) {
		if ( !-s $readsets[$i] && -s $readsets[$i] . '.bam' ) {
			$readsets[$i] .= '.bam';
		}
		elsif ( !-s $readsets[$i] && -s $readsets[$i] . '.bz2' ) {
			$readsets[$i] .= '.bz2';
		}
		pod2usage "File " . $readsets[$i] . " not found\n"
		  unless $readsets[$i] && -s $readsets[$i];
		$readset_lookup{ $readsets[$i] } = 1;
		if (@use_existing_bam) {
			die "Cannot find user-provided BAM file number " . ( $i + 1 ) . "\n"
			  unless $use_existing_bam[$i] && -s $use_existing_bam[$i];
		}
	}
	if ($doing_paired_end) {
		for ( my $i = 0 ; $i < @readsets2 ; $i++ ) {
			if ( !-s $readsets2[$i] && -s $readsets2[$i] . '.bam' ) {
				$readsets2[$i] .= '.bam';
			}
			elsif ( !-s $readsets2[$i] && -s $readsets2[$i] . '.bz2' ) {
				$readsets2[$i] .= '.bz2';
			}
			pod2usage "File " . $readsets2[$i] . " not found\n"
			  unless $readsets2[$i] && -s $readsets2[$i];
			$readset_lookup{ $readsets[$i] } = $readsets2[$i];
		}
	}
	die "Sorry, paired end Bowtie does not work with BAM files\n"
	  if $doing_paired_end && $read_format eq 'bam' && !$use_bwa;
	if ($extra_genes) {
		pod2usage "Cannot find $extra_genes file\n" unless -s $extra_genes;
	}
	if ( -d $result_dir && !$overwrite_results ) {
		pod2usage
"Result dir $result_dir already exists. If you wish to overwrite give the option -over otherwise delete the directory\n";
	}
	mkdir($result_dir)           unless -d $result_dir;
	mkdir($edgeR_dir)            unless -d $edgeR_dir;
	mkdir($result_dir."graphs") unless -d $result_dir."graphs";
	open( LOG, ">>$result_dir/$uid.log" ) || die($!);
	print LOG "#Started: " . &mytime . "\n";

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
		die
"Library alias ($lib_alias_file) must have 'file' and 'name' as the first two columns."
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
			if ( !-s $data[0] && $readset_lookup{ $data[0] } ) {
				die "Cannot find file " . $data[0] . "\n";
			}
			elsif ( !$readset_lookup{ $data[0] } ) {
				warn "Lib alias entry "
				  . $data[0]
				  . " not in readset request. Skipping\n";
				next;
			}
			elsif ( !$data[1] ) {
				warn "Lib alias entry "
				  . $data[0]
				  . " has no name entry. Will use basename of filename\n";
				$data[1] = basename( $data[0] );
			}
			$data[1] =~ s/\W+/_/g;

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
				$print .= $ln . "\t" . $data[1] . "\n";
			}
			else {
				$print .= $ln . "\n";
			}
			die "Readset "
			  . $data[1]
			  . " has already been linked with group "
			  . $library_metadata{ $data[1] }{'group'} . "\n"
			  if $groups_readsets{ $library_metadata{ $data[1] }{'group'} }
			  { $data[1] };

			# for edgeR counts
			my $computer_friendly_name = $data[1];
			$computer_friendly_name =~ s/\W/_/g;
			$groups_readsets{ $library_metadata{ $data[1] }{'group'} }
			  { $data[1] } = $edgeR_dir . $computer_friendly_name . '.dat';
		}
		close IN;

	}
	else {
		$print = "file\tname\tgroup\n";
		foreach my $readset (@readsets) {
			$library_metadata{$readset}{'group'} = $readset;
			$groups_readsets{$readset}{$readset}++;
			$print .= $readset . "\t" . $readset . "\t" . $readset . "\n";
		}
	}

	open( OUT, ">$result_dir/lib_alias.txt" );
	print OUT $print;
	close OUT;
	@sample_overlays = sort keys %library_metadata_headers;

}

sub prepare_input_data() {
	print &mytime . "Preparing readsets\n";
	&prepare_library_alias();

	for ( my $i = 0 ; $i < @readsets ; $i++ ) {
		if ($doing_paired_end) {
			my ( $readset_name, $library_size ) =
			  &perform_readset_metadata( $readsets[$i], $readsets2[$i] );
		}
		else {
			my ( $readset_name, $library_size ) =
			  &perform_readset_metadata( $readsets[$i] );
		}
	}
	print "\n" . &mytime . "Preparing references...\n";

	# prepare alignment file.
	my $file_to_align = $result_dir . "$uid.toalign";

	#unlink($file_to_align) if -s $file_to_align;
	if ( !-s $file_to_align ) {

		# create input file if user provided sequence on cmd line
		my $reference_file = $input_reference_file;
		if ($user_ref_sequence) {
			$reference_file = $result_dir . "$uid.query";
			$user_ref_sequence =~ s/\\n/\n/g;
			die "Sequence from user contains non-ATCGN characters!\n"
			  if $user_ref_sequence =~ /[^ATCGN]/;
			open( FASTA, ">$reference_file" ) || die($!);
			print FASTA ">query\n" unless ( $user_ref_sequence =~ /^>/ );
			print FASTA uc($user_ref_sequence) . "\n";
			close FASTA;
		}
		die "Nothing to align!" unless $reference_file && -s $reference_file;
		copy( $reference_file, $file_to_align );
		if ($extra_genes) {
			open( FASTA, ">>$file_to_align" ) || die($!);
			open( IN,    $extra_genes )       || die($!);
			while ( my $ln = <IN> ) {
				$ln = uc($ln) unless $ln =~ /^>/;
				print FASTA $ln;
			}
			close FASTA;
		}
		unlink($reference_file) if $user_ref_sequence;
	}
	die "Cannot find input FASTA: $file_to_align\n" unless -s $file_to_align;
	if ( -s "$file_to_align.checked" && -s $db_file ) {
		open( IN, "$file_to_align.md5" );
		while ( my $ln = <IN> ) {
			chomp($ln);
			my @data = split( "\t", $ln );
			next unless $data[1];
			push( @reference_sequence_list, $data[1] );
			$user_alias{ $data[0] } = $data[1];
		}
		close IN;

   # we assume we don't need to run the sql commands because they already exist.
	}
	else {
		open( MD5SUMS, ">$file_to_align.md5" );
		open( BED, ">$file_to_align.bed" ) || die($!);
		my $reference_file_obj =
		  Bio::SeqIO->new( -file => $file_to_align, -format => 'fasta' );
		while ( my $seq_obj = $reference_file_obj->next_seq() ) {
			my $id = $seq_obj->id();
			push( @reference_sequence_list, $id );
			my $sequence = uc( $seq_obj->seq() );
			die "No sequence for $id!\n" if ( !$sequence );
			die "Sequence for $id contains non-ATCGN characters!\n"
			  if $sequence =~ /[^ATCGN]/;
			my $seq_md5hash = md5_hex($sequence);
			$id = 'Query' unless $id;
			my $seq_length = $seq_obj->length();
			print BED "$id\t1\t$seq_length\t$id\t0\t+\n";
			print MD5SUMS "$seq_md5hash\t$id\n";
			&sqlite_add_seq_md5( $id, $seq_md5hash, $seq_length );
			die "Your data has sequences that are identical! $id vs "
			  . $user_alias{$seq_md5hash} . "\n"
			  if $user_alias{$seq_md5hash};
			$user_alias{$seq_md5hash} = $id;

			for ( my $i = 0 ; $i < ( scalar(@readsets) ) ; $i++ ) {
				&sqlite_init_expression_statistics( $seq_md5hash,
					$readsets[$i] );
				for ( my $k = $i + 1 ; $k < scalar(@readsets) ; $k++ ) {
					next if $readsets[$i] eq $readsets[$k];
					&sqlite_start_fold_change( $seq_md5hash, $readsets[$i],
						$readsets[$k] );
				}
			}
		}
		close MD5SUMS;
		close BED;
		open( CHECK, ">$file_to_align.checked" );
		print CHECK "Done\n";
		close CHECK;
	}
	&sqlite_backup(1) unless $db_use_file;

	# prepare file for alignments
	if ($use_kanga) {
		unless ( -s $file_to_align . '.kangax' ) {
			print "\t\t\tBuilding reference file for kanga...\r";
			&process_cmd(
"$kangax_exec -i $file_to_align -o $file_to_align.kangax -r $file_to_align -t $file_to_align  2> /dev/null >/dev/null"
			);
			print " Done!\n";
		}
	}
	elsif ($use_bwa) {
		unless ( -s $file_to_align . '.bwt' ) {
			print "\t\t\tBuilding reference file for bwa...\r";
			&process_cmd(
				"$bwa_exec index -a is $file_to_align 2> /dev/null >/dev/null");
			print " Done!\n";
		}
	}
	elsif (@use_existing_bam) {
		print "Will use user-provided BAM (read-name sorted) files\n";
	}
	else {
		unless ( -s "$file_to_align.1.bt2" ) {
			my $build_exec = $bowtie2_exec . '-build';
			print "\t\t\tBuilding reference file for bowtie2...\r";
			&process_cmd(
"$build_exec --offrate 1 $file_to_align $file_to_align >/dev/null 2>> $result_dir/$uid.log"
			);
			print " Done!\n";
		}
	}
	if ($prepare_input_only) {
		warn "User requested to stop after preparing files\n";
		sqlite_destroy();
		exit;
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
	my $todo;
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
		my $readset_basename = basename($readset);
		my $alnbase          = basename($baseout) . '_vs_' . $readset_basename;
		my $alnout_text      = $alnbase . '.bam';
		print "Paired readsets "
		  . $readsets[$i] . " and "
		  . $readsets2[$i]
		  . " : $alnout_text\n"
		  if $readsets2[$i];
		print "Unpaired readset " . $readsets[$i] . " : $alnout_text\n"
		  if !$readsets2[$i];
		print "\n";
	}
	print "\tNo alignments need to be done. Processing all existing...\n"
	  if !$todo;
	print "\t$todo alignments need to be done...\n" if $todo;
	my ( %alignment_sam_files, %alignment_bam_files, $aligned_ids_hashref,
		%hash );
	if ($contextual_alignment) {
		for ( my $i = 0 ; $i < (@readsets) ; $i++ ) {
			print '.' x ( $i + 1 ) . "\r";
			my $readset          = $readsets[$i];
			my $readset_basename = basename($readset);
			my ( $alignment_bam_file, $alignment_sam_file );
			if ( $already_done_alignments{$readset} ) {
				my $alnbase = $baseout . '_vs_' . $readset_basename;
				$alignment_bam_file = $alnbase . '.bam';
				$alignment_sam_file = $alnbase . '.sam';
				my $express_results = $alignment_bam_file . ".express.results";
				&process_express_bias( $express_results, $readset );
			}
			elsif (@use_existing_bam) {    # name sorted
				( $alignment_bam_file, $alignment_sam_file ) =
				  &prepare_alignment_from_existing( $file_to_align, $i );
			}
			else {
				( $alignment_bam_file, $alignment_sam_file ) =
				  &prepare_alignment( $file_to_align, $i );
			}
			$aligned_ids_hashref =
			  &process_depth_of_coverage( $file_to_align, $readset,
				$alignment_bam_file );
			$alignment_sam_files{$readset} = $alignment_sam_file;
			$alignment_bam_files{$readset} = $alignment_bam_file;

			# so if $aligned_ids_hashref from process_depth above
			# is not equal to the references added, then re-enter depth data
			my $do_backup;
			if (
				keys %{$aligned_ids_hashref} !=
				scalar(@reference_sequence_list) )
			{
				foreach my $id (@reference_sequence_list) {
					next if $aligned_ids_hashref->{$id};
					my $seq_md5hash = &sqlite_get_md5($id);
					unless ( ( $resume || $debug )
						&& &sqlite_get_depth_data( $seq_md5hash, $readset ) )
					{
						&sqlite_add_depth_data( $seq_md5hash, $readset,
							\%hash );
						$do_backup++;
					}
				}
			}
			&sqlite_backup() if !$db_use_file && $do_backup;
		}
		print "\n";
		&perform_kangade( $file_to_align, \%alignment_sam_files,
			\%alignment_bam_files )
		  unless $no_kangade;
	}
	else {

# due to kangade only one aln file so we have re-align if it has not been tested for a readset
		my $new_file_to_align = $file_to_align . '_unaligned';
		my $out_obj =
		  Bio::SeqIO->new( -file => ">$new_file_to_align", -format => 'fasta' );
		my (%unaligned);
		for ( my $i = 0 ; $i < (@readsets) ; $i++ ) {
			my $file_obj =
			  Bio::SeqIO->new( -file => $file_to_align, -format => 'fasta' );
			while ( my $seq_obj = $file_obj->next_seq() ) {
				my $id          = $seq_obj->id();
				my $seq_md5hash = &sqlite_get_md5($id);
				my $aln_exists =
				  &sqlite_get_depth_data( $seq_md5hash, $readsets[$i] );
				if ( !$aln_exists ) {
					$out_obj->write_seq($seq_obj) if !$unaligned{$seq_md5hash};
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
		undef($out_obj);
		link( $file_to_align . '.bed', $new_file_to_align . '.bed' );
		if ( -s $new_file_to_align ) {
			for ( my $i = 0 ; $i < (@readsets) ; $i++ ) {
				print '.' x ( $i + 1 ) . "\r";
				my ( $alignment_bam_file, $alignment_sam_file ) =
				  &prepare_alignment( $new_file_to_align, $readsets[$i],
					$readsets2[$i] );

				#too slow:
				$aligned_ids_hashref =
				  &process_depth_of_coverage( $new_file_to_align, $readsets[$i],
					$alignment_bam_file );
				$alignment_sam_files{ $readsets[$i] } = $alignment_sam_file;
				$alignment_bam_files{ $readsets[$i] } = $alignment_bam_file;
				if (
					keys %{$aligned_ids_hashref} !=
					scalar(@reference_sequence_list) )
				{
					foreach my $id (@reference_sequence_list) {
						next if $aligned_ids_hashref->{$id};
						my $seq_md5hash = &sqlite_get_md5($id);
						unless (
							&sqlite_get_depth_data(
								$seq_md5hash, $readsets[$i]
							)
						  )
						{
							&sqlite_add_depth_data( $seq_md5hash, $readsets[$i],
								\%hash );
						}
					}
				}
			}
			print "\n";
		}
		&perform_kangade( $new_file_to_align, \%alignment_sam_files,
			\%alignment_bam_files )
		  unless $no_kangade;
	}
}

sub prepare_alignment_from_existing() {
	my ( $file_to_align, $i ) = @_;
	my $baseout = $file_to_align;
	$baseout =~ s/.toalign$//;
	die "Given BAM alignment file does not exist: "
	  . $use_existing_bam[$i] . "\n"
	  unless $use_existing_bam[$i] && -s $use_existing_bam[$i];
	my $readset          = $readsets[$i];
	my $readset2         = $readsets2[$i];
	my $readset_basename = basename($readset);
	my $alnbase          = $baseout . '_vs_' . $readset_basename;
	my $alignment_sam_file = $alnbase . '.sam';    # reference sorted
	my $alignment_bam_file = $alnbase . '.bam';    # reference sorted
	print "\n\n" . &mytime
	  . "Aligning: Using user-provided alignment for $readset_basename\n";
	print LOG "\n" . &mytime
	  . "Aligning: Using user-provided alignment for $readset_basename\n";

	# namesorted bam
	# hard link:
	link( $use_existing_bam[$i], "$alignment_bam_file.namesorted" )
	  unless -s "$alignment_bam_file.namesorted";

	# in case hard link cannot occur:
	symlink( $use_existing_bam[$i], "$alignment_bam_file.namesorted" )
	  unless -s "$alignment_bam_file.namesorted";
	die "Cannot link "
	  . $use_existing_bam[$i]
	  . " as $alignment_bam_file.namesorted\n"
	  unless -s "$alignment_bam_file.namesorted";

#         &process_cmd("samtools view -h -o $alignment_sam_file.namesorted ".$use_existing_bam[$i]) unless -s "$alignment_sam_file.namesorted";
	&process_cmd(
		"samtools sort -m 18000000000 " . $use_existing_bam[$i] . " $alnbase" )
	  unless -s $alignment_bam_file;
	&perform_correct_bias( $alignment_bam_file, $file_to_align, $readset,
		$alignment_bam_file . '.namesorted' )
	  if ($perform_bias_correction);
	&process_alignments( $alignment_bam_file, $readset, $readset2 );
	return ( $alignment_bam_file, $alignment_sam_file );
}

sub prepare_alignment() {
	my ( $file_to_align, $i ) = @_;
	my $readset  = $readsets[$i];
	my $readset2 = $readsets2[$i];
	my ( $rpkm_hashref,       $eff_counts_hashref );
	my ( $alignment_bam_file, $alignment_sam_file ) =
	  &perform_alignments( $file_to_align, $readset, $readset2 );
	( $alignment_bam_file, $rpkm_hashref, $eff_counts_hashref ) =
	  &perform_correct_bias( $alignment_bam_file, $file_to_align, $readset,
		$alignment_sam_file . '.namesorted' )
	  if ($perform_bias_correction);
	&process_alignments( $alignment_bam_file, $readset, $readset2 );
	return ( $alignment_bam_file, $alignment_sam_file );
}

sub perform_alignments() {
	my $file_to_align    = shift;
	my $readset          = shift;
	my $readset_basename = basename($readset);
	my $readset2         = shift;
	print "\n" . &mytime
	  . "Aligning: Performing alignments of $file_to_align vs $readset_basename\n";
	my $readset_time  = stat($readset)->mtime;
	my $readset2_time = $readset2 ? stat($readset2)->mtime : int(0);
	my $baseout       = $file_to_align;
	$baseout =~ s/.toalign$//;
	$baseout .= '_vs_' . $readset_basename;
	my $bam = $baseout . '.bam';
	my $sam = $baseout . '.sam';
	print LOG &mytime . "Processing $readset as $bam\n";

	# check if alignments need to be recalculated
	if (
		-s $bam
		&& (   stat($bam)->mtime > $readset_time
			&& stat($bam)->mtime > $readset2_time )
	  )
	{
		return ( $bam, $sam );
	}

	# DB method
	# We don't really want to store the alignments, do we!?
	#my $alignment_exists = &sqlite_check_align($readset);
	#if ($alignment_exists) {
	#  &sqlite_get_align( $readset, $bam );
	#} else {
	#  modtime method
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
		&& (   stat($bam)->mtime < $readset_time
			|| stat($bam)->mtime < $readset2_time )
	  );
	unlink($sam)
	  if (
		-s $sam
		&& (   stat($sam)->mtime < $readset_time
			|| stat($sam)->mtime < $readset2_time )
	  );
	unless ( -s $bam && ( -s $bam ) > 1000 ) {
		if ($use_bwa) {
			&align_bwa( $baseout, $file_to_align, $readset, $readset2, $bam,
				$sam );
		}
		elsif ($use_kanga) {
			&align_kanga( $baseout, $file_to_align, $readset, $readset2, $bam,
				$sam );
		}
		else {
			&align_bowtie2( $baseout, $file_to_align, $readset, $readset2, $bam,
				$sam );
		}
		die "Could not produce SAM file for $readset\n" unless -s $sam;
		&process_cmd(
"$samtools_exec view -S  -u $sam 2>/dev/null|samtools sort -m 18000000000 - $baseout"
		) unless -s $bam;
	}
	die "Could not convert to BAM file for $readset\n" unless -s $bam;
	&process_cmd("$samtools_exec index $bam")
	  unless -s $bam . '.bai' && ( -s $bam . '.bai' ) > 200;
	return ( $bam, $sam );
}

sub process_depth_of_coverage($$$) {
	my $file_to_align = shift;
	my $readset       = shift;
	my $bam           = shift;
	return if ( -s $bam . ".depth.completed" );

	my $readset_metadata = &sqlite_get_readset_metadata($readset);
	my $readset_name     = $readset_metadata->{'alias'};
	print LOG &mytime
	  . "Processing $readset_name depth file $bam from $file_to_align\n";
	print &mytime
	  . "Processing $readset_name depth file $bam from $file_to_align\n";
	my $tmp_depth_file = $bam . ".depth";
	&process_cmd("$samtools_exec depth $bam > $tmp_depth_file 2>/dev/null")
	  unless -s $tmp_depth_file;

	unless ( -s $tmp_depth_file > 100 ) {
		warn "Could not produce depth file for $bam. Skipping\n";
		print LOG "Could not produce depth file for $bam. Skipping\n";
		unlink($tmp_depth_file);
		return;
	}
	my ( %hash, %aligned_ids );

#TODO this is very slow when we have lots and lots of files. concatanate and process once?
	my $timer_counter = int(0);
	my $timer         = new Time::Progress;
	$timer->attr( min => 0, max => -s $tmp_depth_file );
	open( DEPTH, $tmp_depth_file ) || die($!);
	my $previous_id;
	while ( my $ln = <DEPTH> ) {
		$timer_counter += length($ln);
		print $timer->report( "eta: %E min, %40b %p\r", $timer_counter )
		  if ( $timer_counter % 10000 == 0 );
		chomp($ln);
		my @data = split( "\t", $ln );
		if ( $previous_id && $previous_id ne $data[0] ) {
			my $seq_md5hash = &sqlite_get_md5($previous_id);
			$aligned_ids{$previous_id} = 1;
			$previous_id = $data[0];

		 # do not add if they already exist (unless it is a contextual alignment
		 # without debug, in which case overwrite
			unless ( ( !$contextual_alignment || $debug || $resume )
				&& &sqlite_get_depth_data( $seq_md5hash, $readset ) )
			{
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
	unless ( ( !$contextual_alignment || $debug || $resume )
		&& &sqlite_get_depth_data( $seq_md5hash, $readset ) )
	{
		&sqlite_add_depth_data( $seq_md5hash, $readset, \%hash );
	}
	close DEPTH;
	unlink($tmp_depth_file) unless $debug;
	open( OUT, ">$bam.depth.completed" );
	print OUT &mytime();
	close OUT;
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
					my $readset_C = $readsets[$i];
					for ( my $k = $i + 1 ; $k < scalar(@readsets) ; $k++ ) {
						my $readset_E = $readsets[$k];
						my $exists =
						  &sqlite_check_kangade_analysis( $seq_md5hash,
							$readset_C, $readset_E );
						next if !$exists;
						my $readset_metadata_C =
						  &sqlite_get_readset_metadata($readset_C);
						my $readset_metadata_E =
						  &sqlite_get_readset_metadata($readset_E);
						my $readset_name_C = $readset_metadata_C->{'alias'};
						my $readset_name_E = $readset_metadata_E->{'alias'};
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
						$add_kangade_fold_change->execute(
							$exists->{'ObsFoldChange'},
							$seq_md5hash, $readset_C, $readset_E );
						$fold_changes{$seq_md5hash}{'kangade'}{$readset_name_C}
						  {$readset_name_E} = $exists->{'ObsFoldChange'};
						&sqlite_add_kangade_expression_statistics( $seq_md5hash,
							$readset_C, $exists->{'TotCtrlCnts'} );
						&sqlite_add_kangade_expression_statistics( $seq_md5hash,
							$readset_E, $exists->{'TotExprCnts'} )
						  if $i == scalar(@readsets) - 2;
						$md5_not_to_process{$seq_md5hash} = 1;
					}
				}
			}
		}

# now perform kangade against all data but parse only those that were not stored.
# i know that is it not the most efficient...
		for ( my $i = 0 ; $i < ( scalar(@readsets) - 1 ) ; $i++ ) {
			my $readset_C          = $readsets[$i];
			my $readset_metadata_C = &sqlite_get_readset_metadata($readset_C);
			my $readset_name_C     = $readset_metadata_C->{'alias'};
			my $readset_reads_C    = $readset_metadata_C->{'total_reads'};
			my $readset_basename_C = basename($readset_C);
			print &mytime() . "\tProcessing $readset_name_C\n";
			for ( my $k = $i + 1 ; $k < scalar(@readsets) ; $k++ ) {
				my $readset_E = $readsets[$k];
				my $readset_metadata_E =
				  &sqlite_get_readset_metadata($readset_E);
				my $readset_name_E     = $readset_metadata_E->{'alias'};
				my $readset_basename_E = basename($readset_E);
				my $out =
				  $baseout . $readset_basename_C . '_vs_' . $readset_basename_E;
				unless ( -s "$out.stats.csv" ) {
					my $readset_reads_E = $readset_metadata_E->{'total_reads'};
					my $library_ratio =
					  sprintf( "%.2f",
						( $readset_reads_C / $readset_reads_E ) );
					my $kanga_threads = $threads <= 8 ? $threads : 8;
					my $alignment_sam_file_C =
					  $alignment_samfiles_hashref->{$readset_C};
					if ( ( !$alignment_sam_file_C || !-s $alignment_sam_file_C )
						&& $alignment_bamfiles_hashref->{$readset_C} )
					{
						$alignment_sam_file_C =
						  $alignment_bamfiles_hashref->{$readset_C} . '.sam';
						&process_cmd(
							"samtools view -h -o $alignment_sam_file_C "
							  . $alignment_bamfiles_hashref->{$readset_C} )
						  unless -s $alignment_sam_file_C;
						$files_to_delete{$alignment_sam_file_C} = 1;
					}
					elsif (!$alignment_sam_file_C
						|| !-s $alignment_sam_file_C )
					{
						warn("Cannot find alignment SAM file for $readset_C\n");
						next;
					}
					my $alignment_sam_file_E =
					  $alignment_samfiles_hashref->{$readset_E};
					if ( ( !$alignment_sam_file_E || !-s $alignment_sam_file_E )
						&& $alignment_bamfiles_hashref->{$readset_E} )
					{
						$alignment_sam_file_E =
						  $alignment_bamfiles_hashref->{$readset_E} . '.sam';
						&process_cmd(
							"samtools view -h -o $alignment_sam_file_E "
							  . $alignment_bamfiles_hashref->{$readset_E} )
						  unless -s $alignment_sam_file_E;
						$files_to_delete{$alignment_sam_file_E} = 1;
					}
					elsif (!$alignment_sam_file_E
						|| !-s $alignment_sam_file_E )
					{
						warn("Cannot find alignment SAM file for $readset_E\n");
						next;
					}
					&process_cmd(
"$kangade_exec -r 0 -T $kanga_threads -t 3 -g $file_to_align.bed -i $alignment_sam_file_C -I $alignment_sam_file_E -f3 -F $out.log -o $out.stats.csv -O $out.bins.csv >/dev/null 2>/dev/null"
					);

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
				  || die "Cannot use CSV: " . Text::CSV_XS->error_diag();
				open( my $fh, "$out.stats.csv" ) || die($!);

		  #my $h = <$fh>;chomp($h);$h =~ s/"//g;my @headers = split( "\t", $h );
				$csv->column_names( $csv->getline($fh) );
				while ( my $row = $csv->getline_hr($fh) ) {
					my $seq_md5hash = &sqlite_get_md5( $row->{'Feat'} );
					next if $md5_not_to_process{$seq_md5hash};
					$row->{'ObsFoldChange'} =
					  sprintf( "%.2f", $row->{'ObsFoldChange'} );
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
					$fold_changes{$seq_md5hash}{'kangade'}{$readset_name_C}
					  {$readset_name_E} = $row->{'ObsFoldChange'};
					&sqlite_add_kangade_expression_statistics( $seq_md5hash,
						$readset_C, $row->{'TotCtrlCnts'} );
					&sqlite_add_kangade_expression_statistics( $seq_md5hash,
						$readset_E, $row->{'TotExprCnts'} )
					  if $i == scalar(@readsets) - 2;
					my $de_exists =
					  &sqlite_check_kangade_analysis( $seq_md5hash, $readset_C,
						$readset_E );
					next if $de_exists;
					$add_to_kangade_analysis->execute(
						$seq_md5hash,             $readset_C,
						$readset_E,               $row->{'Classification'},
						$row->{'Score'},          $row->{'DECntsScore'},
						$row->{'PearsonScore'},   $row->{'CtrlUniqueLoci'},
						$row->{'ExprUniqueLoci'}, $row->{'CtrlExprLociRatio'},
						$row->{'PValueMedian'},   $row->{'PValueLow95'},
						$row->{'PValueHi95'},     $row->{'TotCtrlCnts'},
						$row->{'TotExprCnts'},    $row->{'TotCtrlExprCnts'},
						$row->{'ObsFoldChange'},  $row->{'FoldMedian'},
						$row->{'FoldLow95'},      $row->{'FoldHi95'},
						$row->{'ObsPearson'},     $row->{'PearsonMedian'},
						$row->{'PearsonLow95'},   $row->{'PearsonHi95'}
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
			$fold_changes{ $data[0] }{'kangade'}{ $data[2] }{ $data[3] } =
			  $data[9];
		}
		close IN;
	}
}

sub perform_readset_metadata($$) {
	my $readset      = shift;
	my $readset2     = shift;
	my $library_size = int(0);
	my $readset_name =
	    $library_aliases{$readset}
	  ? $library_aliases{$readset}
	  : basename($readset);
	$get_readset->execute($readset);
	my $result = $get_readset->fetchrow_hashref();
	if ($result) {
		return ( $result->{'alias'}, $result->{'total_reads'} );
	}
	print "Counting reads in readset $readset / $readset2\n";

	# get library sizes
	if ( $read_format eq 'fastq' ) {
		if ( $readset =~ /.bz2$/ ) {
			$library_size = `$bunzip2_exec -dck $readset| wc -l `;
			chomp($library_size);
			$library_size += `$bunzip2_exec -dck $readset2| wc -l `
			  if $doing_paired_end;
		}
		else {
			$library_size = `wc -l < $readset`;
			chomp($library_size);
			$library_size += `wc -l < $readset2` if $doing_paired_end;
		}
		$library_size /= 4;
		print LOG "Reads in $readset / $readset2 PE files: $library_size\n"
		  if $doing_paired_end;
		print LOG "Reads in $readset file: $library_size\n"
		  if !$doing_paired_end;
	}
	elsif ( $read_format eq 'bam' ) {
		die "BAM $readset not found\n" unless -s $readset;
		&process_cmd("$samtools_exec index $readset")
		  unless -s $readset . '.bai';
		die "Readset $readset cannot be indexed" unless -s $readset . '.bai';
		my $d = `samtools idxstats $readset|grep '^\*'`;
		$d =~ /(\d+)$/;
		$library_size = $1;
		if ($doing_paired_end) {
			$d = 0;
			$d = `samtools idxstats $readset2|grep '^\*'`;
			$d =~ /(\d+)$/;
			$library_size += $1;
		}
		print LOG "Reads in $readset file: $library_size\n"
		  if !$doing_paired_end;
		print LOG "Reads in $readset / $readset2 PE files: $library_size\n"
		  if $doing_paired_end;
	}
	die "Cannot get library size for $readset"
	  unless $library_size && $library_size > 0;
	&sqlite_add_readset_metadata( $readset, $readset_name, $library_size );
	return ( $readset_name, $library_size );
}

sub align_bowtie2() {
	my ( $baseout, $file_to_align, $readset, $readset2, $bam, $sam ) = @_;
	my $build_exec = $bowtie2_exec . '-build';
	unless ( -s "$file_to_align.1.bt2" ) {
		print "\t\t\tBuilding reference file for bowtie...\r";
		&process_cmd(
"$build_exec --offrate 1 $file_to_align $file_to_align >/dev/null 2>> $result_dir/$uid.log"
		);
		print " Done!\n";
	}
	my $readgroup = '@RG\tID:' . $readset;
	$readgroup =~ s/\.bz2$//;
	$readgroup =~ s/\.fastq$//;
	$readgroup =~ s/\.bam$//;
	if ($doing_paired_end) {
		if ( $read_format eq 'fastq' ) {
			my $qformat = '--' . &check_fastq_format($readset);
			&process_cmd(
"$bowtie2_exec $qformat --end-to-end --sensitive -qap $threads -I 100 -X 800 --no-mixed --no-unal -x $file_to_align -1 $readset -2 $readset2 -S $sam 2> $baseout.log"
			) unless -s "$baseout.log";
		}
		else {
			die "Sorry, paired end Bowtie does not work with BAM files\n";

#&process_cmd("$bamtools_exec convert -in $readset -format fastq | $bowtie2_exec --end-to-end --fast -qap $threads -x $file_to_align -U - -S $sam");
#&process_cmd("$bamtools_exec convert -in $readset2 -format fastq | $bowtie2_exec --end-to-end --fast -qap $threads -x $file_to_align -U - >> $sam");
		}
	}
	else {
		if ( $read_format eq 'fastq' ) {
			my $qformat = '--' . &check_fastq_format($readset);
			&process_cmd(
"$bowtie2_exec $qformat --end-to-end --sensitive -qap $threads -x $file_to_align -U $readset -S $sam 2> $baseout.log >/dev/null"
			) unless -s "$baseout.log";
		}
		else {
			&process_cmd(
"$bamtools_exec convert -in $readset -format fastq | $bowtie2_exec --end-to-end --no-unal --sensitive -qap $threads -x $file_to_align -U - -S $sam 2> $baseout.log"
			) unless -s "$baseout.log";
		}
	}
	die "Alignment for $file_to_align vs $readset failed\n" unless -s $sam;
	chdir($result_dir);
	link( basename($sam), basename($baseout) . ".sam.namesorted" );
	chdir("../");
}

sub fix_check_kanga() {

	#	####temp!!! debug
	my $in = shift;

	#	my $errors;
	#	die unless $in && -s $in;
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
	&process_cmd("sort -T /tmp -S 20% -k1,1 -k3,3 -o $in.sorted $in");
	unlink("$in");
	rename( "$in.sorted", $in );
}

sub align_kanga() {
	my ( $baseout, $file_to_align, $readset, $readset2, $bam, $sam ) = @_;
	print "\t\t\tBuilding reference file for kanga...\r";
	&process_cmd(
"$kangax_exec -i $file_to_align -o $file_to_align.kangax -r $file_to_align -t $file_to_align  2> /dev/null >/dev/null"
	) unless -s $file_to_align . '.kangax';
	print " Done!\n";
	if ( $read_format eq 'fastq' ) {
		if ( $readset =~ /\.bz2$/ ) {
			die "Kanga does not support .bz2 FASTQ files.\n";
		}
		my $qformat =
		  ( &check_fastq_format($readset) eq 'phred33' ) ? ' -q0 ' : ' -q1 ';
		if ($doing_paired_end) {
			unless ( -s "$sam.1" ) {
				&process_cmd(
"$kanga_exec -I $file_to_align.kangax -m 0 $qformat -R 100 -r 5 -X -M 5 -i $readset -o $sam.1 -F $sam.1.log -T $threads -d 100 -D 800 >/dev/null"
				);
				&fix_check_kanga("$sam.1");
			}
			unless ( -s "$sam.2" ) {
				&process_cmd(
"$kanga_exec -I $file_to_align.kangax -m 0 $qformat -R 100 -r 5 -X -M 5 -i $readset2 -o $sam.2 -F $sam.2.log -T $threads -d 100 -D 800 >/dev/null"
				);
				&fix_check_kanga("$sam.2");
			}
			####temp!!! debug
			&process_cmd(
"merge_left_right_nameSorted_SAMs.pl --left_sam $sam.1 --right_sam $sam.2  -D 800 -C 100 > $sam.t"
			);
			unlink("$sam.1");
			unlink("$sam.2");
			&process_cmd(
				"samtools view -F12 -h -T $file_to_align -S -o $sam $sam.t");
			unlink("$sam.t");
		}
		else {
			unless ( -s "$sam" ) {
				&process_cmd(
"$kanga_exec -I $file_to_align.kangax -m 0 -q 0 -R 100 -r 5 -X -M 5 -i $readset -o $sam -F $sam.log -T $threads -d 100 -D 800 >/dev/null"
				);
				&fix_check_kanga("$sam");
			}
		}
		die "Alignment for $file_to_align vs $readset failed\n" unless -s $sam;
		chdir($result_dir);
		link( basename($sam), basename($baseout) . ".sam.namesorted" );
		chdir("../");
	}
	else {
		die "Kanga is currently only supporting FASTQ\n";
	}
}

sub align_bwa() {
	my ( $baseout, $file_to_align, $readset, $readset2, $bam, $sam ) = @_;
	print "\t\t\tBuilding reference file for bwa...\r";
	&process_cmd("$bwa_exec index -a is $file_to_align 2> /dev/null >/dev/null")
	  unless -s $file_to_align . '.bwt';
	print " Done!\n";
	if ( $read_format eq 'fastq' ) {
		my $qformat =
		  ( &check_fastq_format($readset) eq 'phred33' ) ? ' ' : '-I';
		&process_cmd(
"$bwa_exec aln $qformat -t $threads -f $baseout.sai -q 10 $file_to_align $readset 2>/dev/null"
		) unless -s "$baseout.sai";
		if ($doing_paired_end) {
			&process_cmd(
"$bwa_exec aln $qformat -t $threads -f $baseout.2.sai -q 10 $file_to_align $readset2 2>/dev/null"
			) unless -s "$baseout.2.sai";
		}
	}
	elsif ( $read_format eq 'bam' ) {
		&process_cmd(
"$bwa_exec aln -t $threads -b -f $baseout.sai -q 10 $file_to_align $readset 2>/dev/null"
		) unless -s "$baseout.sai";
		if ($doing_paired_end) {
			&process_cmd(
"$bwa_exec aln -t $threads -b -f $baseout.2.sai -q 10 $file_to_align $readset2 2>/dev/null"
			) unless -s "$baseout.2.sai";
		}
	}
	die "Could not produce BWA SAI for $readset\n" unless -s "$baseout.sai";
	my $readgroup = '@RG\tID:' . $readset;
	$readgroup =~ s/\.bz2$//;
	$readgroup =~ s/\.fastq$//;
	$readgroup =~ s/\.bam$//;
	if ($doing_paired_end) {
		&process_cmd(
"$bwa_exec sampe -n 20 -N 100 -a 700 -s -r '$readgroup' -o $sam $file_to_align $baseout.sai $baseout.2.sai $readset $readset2 2>/dev/null"
		);
		unlink("$baseout.sai");
		unlink("$baseout.2.sai");
	}
	else {
		&process_cmd(
"$bwa_exec samse -n 20 -s -r '$readgroup' -o $sam $file_to_align $baseout.sai $readset 2>/dev/null"
		);
		unlink("$baseout.sai");
	}
	&process_cmd("sort -S30% -k1,1 -k3,3 -o $sam.namesorted $sam ");
}

sub process_alignments($) {
	my ( $alignment_bam, $readset, $readset2 ) = @_;
	print &mytime . "Post-processing alignment $readset\n";
	&process_cmd("$samtools_exec index $alignment_bam")
	  unless -s $alignment_bam . '.bai';
	die "File $alignment_bam cannot be indexed"
	  unless -s $alignment_bam . '.bai';
	my @sizes              = `samtools idxstats $alignment_bam|grep -v '^\*'`;
	my $reads_that_aligned = `samtools view -F260 -c $alignment_bam`;
	chomp($reads_that_aligned);
	my $readset_metadata = &sqlite_get_readset_metadata($readset);

	foreach my $s (@sizes) {
		next if $s =~ /^\*/;   # unaligned reads
		                       # gene_id, length, reads matched, reads,unmatched
		my @data = split( "\t", $s );
		next unless $data[2];

		#    $reads_that_aligned += $data[2];
		my $seq_md5hash = &sqlite_get_md5( $data[0] );
		my $rpkm        = sprintf(
			"%.0f",
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
	my $alignment_proportion = sprintf( "%.2f%%",
		( $reads_that_aligned / $readset_metadata->{'total_reads'} ) * 100 );
	print LOG "Reads that aligned: $reads_that_aligned ("
	  . $alignment_proportion . ")\n";
}

sub process_express_bias() {
	my $express_results = shift;
	my $readset         = shift;
	if ( -s $express_results ) {
		open( EXPRESS, $express_results ) || die($!);
		my $header = <EXPRESS>;
		while ( my $ln = <EXPRESS> ) {
			chomp($ln);
			my @data = split( "\t", $ln );
			next unless $data[10] && $data[7];
			my $seq_md5hash = &sqlite_get_md5( $data[1] );
			&sqlite_add_express_expression_statistics(
				$seq_md5hash, $readset,
				sprintf( "%.2f", $data[10] ),
				sprintf( "%.2f", $data[7] )
			);
		}
	}
}

sub perform_correct_bias($$$) {
	print &mytime
	  . "eXpress: performing Illumina bias and transcript assignment corrections\n";
	my $original_bam      = shift;
	my $fasta_file        = shift;
	my $readset           = shift;
	my $namesorted_sam    = shift;
	my $readset_metadata  = &sqlite_get_readset_metadata($readset);
	my $readset_name      = $readset_metadata->{'alias'};
	my $original_bam_base = basename($original_bam);
	my $namesorted_bam =
	  ( $namesorted_sam =~ /\.bam.namesorted$/ )
	  ? $namesorted_sam
	  : $namesorted_sam . '.bam';
	my $express_bam_base =
	  $result_dir . $original_bam_base . ".express.sampled";
	my $express_bam     = $express_bam_base . '.bam';
	my $express_results = $result_dir . $original_bam_base . ".express.results";
	my $express_dir     = $result_dir . "$readset_name.bias";

	unless ( -s $express_results ) {
		mkdir($express_dir) unless -d $express_dir;
		if ($use_bwa) {
			&process_cmd(
"$samtools_exec sort -n -m 18000000000 $original_bam $namesorted_sam"
			) unless -s $namesorted_bam;
			if ($contextual_alignment) {
				&process_cmd(
"$express_exec --no-update-check -o $express_dir --output-align-samp -B 2 -p 4 -m 250 $fasta_file $namesorted_bam >/dev/null "
				) unless -s "$express_dir/results.xprs";
				sleep(30);
				&process_cmd(
"$samtools_exec sort -m 18000000000 $express_dir/hits.1.samp.bam $express_bam_base"
				) unless -s "$express_bam_base.bam";
				rename( "$express_dir/results.xprs", $express_results );
				rename( "$express_dir/varcov.xprs",
					$express_results . ".varcov" );
			}
			else {
				open( EXPR, ">$express_results" ) || die($!);
				print EXPR
"bundle_id\ttarget_id\tlength\teff_length\ttot_counts\tuniq_counts\test_counts\teff_counts\tambig_distr_alpha\tambig_distr_beta\tfpkm\tfpkm_conf_low\tfpkm_conf_high\tsolvable\n";
				my $fasta_obj =
				  Bio::SeqIO->new( -file => $fasta_file, -format => 'fasta' );
				my $bundle_counter = 1;
				my $timer          = new Time::Progress;
				my $fasta_count    = `grep -c '^>'  < $fasta_file`;
				chomp($fasta_count);
				$timer->attr( min => 0, max => $fasta_count );

				while ( my $seq_obj = $fasta_obj->next_seq() ) {
					print $timer->report( "eta: %E min, %40b %p\r",
						$bundle_counter )
					  if ( $bundle_counter % 100 == 0 );
					my $id          = $seq_obj->id();
					my $seq_md5hash = &sqlite_get_md5($id);
					my $expression_stats =
					  &sqlite_get_expression_statistics( $seq_md5hash,
						$readset );
					next if $expression_stats->{'express_fpkm'};
					my $len = $seq_obj->length();
					my $seq = $seq_obj->seq();
					my $tmp_fasta_file =
					  $express_dir . '/' . $seq_md5hash . '.fasta.tmp';
					my $tmp_sam_file =
					  $express_dir . '/' . $seq_md5hash . '.tmp';
					my $tmp_express_file = "$tmp_sam_file.dir/results.xprs";
					open( FASTA, ">$tmp_fasta_file" ) || die($!);
					print FASTA ">$id\n$seq\n";
					close FASTA;

					#create sam
					open( SAM, ">$tmp_sam_file" ) || die($!);
					print SAM
					  "\@HD\tVN:1.0\tSO:unsorted\n\@SQ\tSN:$id\tLN:$len\n";
					close SAM;

					#sort for express
					&process_cmd(
"$samtools_exec view -F4 -o $tmp_sam_file $original_bam $id  |sort -k1  >>$tmp_sam_file  "
					);
					if ( -s $tmp_sam_file < 200 ) {
						warn "No alignments for $id. Skipping\n" if $debug;
						print LOG "No alignments for $id. Skipping\n";
						unlink($tmp_fasta_file) unless $debug;
						unlink($tmp_sam_file)   unless $debug;
						next;
					}
					&process_cmd(
"$express_exec --no-update-check -o $tmp_sam_file.dir -B 2  -p 4 -m 250 $tmp_fasta_file $tmp_sam_file >/dev/null "
					) unless -s $tmp_express_file;
					sleep(3);
					if ( -s $tmp_express_file < 200 ) {
						warn "No post-express coverage for $id. Skipping\n"
						  if $debug;
						print LOG
						  "No post-express coverage for $id. Skipping\n";
						unlink($tmp_fasta_file)   unless $debug;
						unlink($tmp_sam_file)     unless $debug;
						unlink($tmp_express_file) unless $debug;
						next;
					}
					open( IN, $tmp_express_file || die($!) );
					my $h = <IN>;
					while ( my $ln = <IN> ) {
						$ln =~ s/^\d+/$bundle_counter/;
						print EXPR $ln;
					}
					close IN;
					unlink($tmp_express_file)        unless $debug;
					unlink($tmp_sam_file)            unless $debug;
					unlink($tmp_fasta_file)          unless $debug;
					remove_tree("$tmp_sam_file.dir") unless $debug;
					$bundle_counter++;
				}
				close EXPR;
			}
		}
		else {
			if ($contextual_alignment) {
				if ( -s $namesorted_bam ) {
					&process_cmd(
"$express_exec --no-update-check -o $express_dir --output-align-samp -B 2  -p 4 -m 250 $fasta_file $namesorted_bam >/dev/null "
					) unless -s "$express_dir/results.xprs";
					die "Express failed to produce output\n"
					  unless -s "$express_dir/results.xprs";
					sleep(30);
					&process_cmd(
"$samtools_exec view -u $express_dir/hits.1.samp.bam | samtools sort -m 18000000000 - $express_bam_base"
					) unless -s "$express_bam_base.bam";
				}
				elsif ( -s $namesorted_sam ) {
					&process_cmd(
"$express_exec --no-update-check -o $express_dir --output-align-samp -B 2  -p 4 -m 250 $fasta_file $namesorted_sam >/dev/null "
					) unless -s "$express_dir/results.xprs";
					die "Express failed to produce output\n"
					  unless -s "$express_dir/results.xprs";
					sleep(30);
					&process_cmd(
"$samtools_exec view -S -u $express_dir/hits.1.samp.sam | samtools sort -m 18000000000 - $express_bam_base"
					) unless -s "$express_bam_base.bam";
				}
				rename( "$express_dir/results.xprs", $express_results );

			 # rename( "$express_dir/varcov.xprs", $express_results.".varcov" );
			}
			else {
				my $timer       = new Time::Progress;
				my $fasta_count = `grep -c '^>' < $fasta_file`;
				chomp($fasta_count);
				$timer->attr( min => 0, max => $fasta_count );
				open( EXPR, ">$express_results" ) || die($!);
				print EXPR
"bundle_id\ttarget_id\tlength\teff_length\ttot_counts\tuniq_counts\test_counts\teff_counts\tambig_distr_alpha\tambig_distr_beta\tfpkm\tfpkm_conf_low\tfpkm_conf_high\tsolvable\n";
				my $fasta_obj =
				  Bio::SeqIO->new( -file => $fasta_file, -format => 'fasta' );
				my $bundle_counter = 1;

				while ( my $seq_obj = $fasta_obj->next_seq() ) {
					print $timer->report( "eta: %E min, %40b %p\r",
						$bundle_counter )
					  if ( $bundle_counter % 100 == 0 );
					my $id          = $seq_obj->id();
					my $seq_md5hash = &sqlite_get_md5($id);
					my $expression_stats =
					  &sqlite_get_expression_statistics( $seq_md5hash,
						$readset );
					next if $expression_stats->{'express_fpkm'};
					my $len = $seq_obj->length();
					my $seq = $seq_obj->seq();
					my $tmp_fasta_file =
					  $express_dir . '/' . $seq_md5hash . '.fasta.tmp';
					my $tmp_sam_file =
					  $express_dir . '/' . $seq_md5hash . '.tmp';
					my $tmp_express_file = "$tmp_sam_file.dir/results.xprs";
					open( FASTA, ">$tmp_fasta_file" ) || die($!);
					print FASTA ">$id\n$seq\n";
					close FASTA;

					#create sam
					open( SAM, ">$tmp_sam_file" ) || die($!);
					print SAM
					  "\@HD\tVN:1.0\tSO:unsorted\n\@SQ\tSN:$id\tLN:$len\n";
					close SAM;

					#sort for express
					&process_cmd(
"$samtools_exec view -F4 $original_bam $id |sort -k1  >>$tmp_sam_file "
					);
					if ( !-s $tmp_sam_file || -s $tmp_sam_file < 200 ) {
						warn "No alignments for $id. Skipping\n" if $debug;
						print LOG "No alignments for $id. Skipping\n";
						unlink($tmp_fasta_file) unless $debug;
						unlink($tmp_sam_file)   unless $debug;
						next;
					}
					&process_cmd(
"$express_exec --no-update-check -o $tmp_sam_file.dir --output-align-samp -B 2  -p 4 -m 250 $tmp_fasta_file $tmp_sam_file >/dev/null "
					) unless -s $tmp_express_file;
					sleep(3);
					if ( !-s $tmp_express_file || -s $tmp_express_file < 200 ) {
						warn "No post-express coverage for $id. Skipping\n"
						  if $debug;
						print LOG
						  "No post-express coverage for $id. Skipping\n";
						unlink($tmp_fasta_file)          unless $debug;
						unlink($tmp_sam_file)            unless $debug;
						unlink($tmp_express_file)        unless $debug;
						remove_tree("$tmp_sam_file.dir") unless $debug;
						next;
					}
					open( IN, $tmp_express_file ) || die($!);
					my $h = <IN>;
					while ( my $ln = <IN> ) {
						$ln =~ s/^\d+/$bundle_counter/;
						print EXPR $ln;
					}
					close IN;
					unlink($tmp_express_file)        unless $debug;
					unlink($tmp_sam_file)            unless $debug;
					unlink($tmp_fasta_file)          unless $debug;
					remove_tree("$tmp_sam_file.dir") unless $debug;
					$bundle_counter++;
				}
				close EXPR;
			}
		}
		unlink($namesorted_sam) unless $debug;
		&remove_tree($express_dir) unless $debug;
		if ($contextual_alignment) {
			if ($debug) {
				unlink( $original_bam . '.orig' ) if -s $original_bam;
				rename( $original_bam, $original_bam . '.orig' );
			}
			else {
				unlink($original_bam) if !$debug;
			}
			chdir($result_dir);
			link( basename($express_bam), basename($original_bam) );
			chdir($cwd);
		}
	}

	&process_express_bias( $express_results, $readset );

	if ($contextual_alignment) {
		return ($express_bam);
	}
	else {
		return ($original_bam);
	}
}

sub perform_stats() {
	if ( !$no_graphs ) {
		print "\n"
		  . &mytime()
		  . "Stats n graphs: Calculating per gene stats and creating images in $result_dir/graphs \n";
		print LOG &mytime()
		  . "Calculating per gene stats and creating images in $result_dir/graphs \n";
	}
	else {
		print "\n" . &mytime() . "Calculating per gene stats in $result_dir \n";
		print LOG &mytime() . "Calculating per gene stats in $result_dir \n";
	}
	my $timer_counter = int(0);
	my $timer         = new Time::Progress;
	$timer->attr( min => 0, max => scalar(@reference_sequence_list) );
	open( STATS,       ">$result_dir/$uid.raw_stats.tsv" )   || die($!);
	open( STATS_RATIO, ">$result_dir/$uid.ratio.stats.tsv" ) || die($!);
	print STATS_RATIO
	  "Checksum\tGene alias\tReadset1\tReadset2\tRaw RPKM\tKANGADE";
	print STATS_RATIO "\tExpress-corrected FPKM\tExpress-corrected counts"
	  if $perform_bias_correction;
	print STATS_RATIO "\n";
	print STATS
"Checksum\tReadset\tGene alias\tRPKM\tMean\tstd. dev\tMedian\tUpper limit\n";

	foreach my $id (@reference_sequence_list) {
		&process_main_stats($id);
		$timer_counter++;
		if ( $timer_counter % 100 == 0 ) {
			print $timer->report( "eta: %E min, %40b %p\r", $timer_counter );
			$dbh->do("PRAGMA shrink_memory") if ( $timer_counter % 100 == 0 );
		}
	}
}

sub stats_with_graph_js($$$) {
	## TODO to avoid BioPerl.
}

sub stats_with_graph($$$) {
	my ( $seq_md5hash, $seq_id, $seq_size ) = @_;
	my $graph_file = "$result_dir/graphs/$seq_md5hash.graph.png";
	return if -s $graph_file;

	# I think there is a memory leak somewhere here....!

	# build a panel with one track per readset
	my $graph = Bio::Graphics::Panel->new(
		-length    => $seq_size,
		-width     => 800,
		-pad_left  => 40,
		-pad_right => 40,
	);
	$graph->add_track(
		arrow => Bio::SeqFeature::Generic->new(
			-start        => 1,
			-end          => $seq_size,
			-display_name => $seq_id
		),
		-description => 1,
		-tick        => 2,
		-fgcolor     => 'black',
		-double      => 1,
		-label       => 1,
		-bump        => 0,
	);
	foreach my $readset (@readsets) {
		my $readset_metadata_ref = &sqlite_get_readset_metadata($readset);
		my $readset_name         = $readset_metadata_ref->{'alias'};
		my $expression_statistics_ref =
		  &sqlite_get_expression_statistics( $seq_md5hash, $readset );
		die "Cannot find any coverage statistics for $seq_id vs $readset !\n"
		  unless $expression_statistics_ref;
		my $mean =
		  $expression_statistics_ref->{'mean_hits'}
		  ? int( $expression_statistics_ref->{'mean_hits'} )
		  : int(0);
		my $median =
		  $expression_statistics_ref->{'median_hits'}
		  ? int( $expression_statistics_ref->{'median_hits'} )
		  : int(0);
		my $rpkm =
		  $expression_statistics_ref->{'rpkm'}
		  ? int( $expression_statistics_ref->{'rpkm'} )
		  : int(0);
		my $express_fpkm =
		    $expression_statistics_ref->{'express_fpkm'}
		  ? $expression_statistics_ref->{'express_fpkm'}
		  : int(0);
		my $effective_counts =
		    $expression_statistics_ref->{'express_eff_counts'}
		  ? $expression_statistics_ref->{'express_eff_counts'}
		  : int(0);
		my $stats_description = " mean: $mean median: $median RPKM: $rpkm";
		$stats_description .=
		  " FPKM: $express_fpkm Effective counts: $effective_counts "
		  if $expression_statistics_ref->{'express_fpkm'};
		my $readset_feature = Bio::SeqFeature::Generic->new(
			-display_name => $readset_name,
			-start        => 1,
			-end          => $seq_size,
			-source_tag   => $readset_name . $stats_description,
			-type         => 'expression',
		);

		if ( $expression_statistics_ref->{'total_hits'} ) {
			my $hit_max = $expression_statistics_ref->{'max_hits'};
			my $hit_sd  = $expression_statistics_ref->{'hit_sd'};
			my $depth_hash_ref =
			  &sqlite_get_depth_data( $seq_md5hash, $readset );
			for ( my $pos = 0 ; $pos < $seq_size ; $pos += 10 ) {
				last if $pos >= $seq_size;

#if ( $pos % 10 == 0 ) {
#          my $depth = $depth_hash_ref->{$pos} ? $depth_hash_ref->{$pos} : int(0);
				my ( $xl, $depth );

				#sliding window for average
				for ( my $x = $pos ; $x < $pos + 10 ; $x++ ) {
					$depth +=
					  $depth_hash_ref->{$x} ? $depth_hash_ref->{$x} : int(0);
					$xl++;
				}
				$depth = ( $depth && $depth > 0 ) ? $depth / $xl : int(0);
				my $t = Bio::SeqFeature::Lite->new(
					-start => $pos + 1,
					-stop  => $pos + 1,
					-type  => 'depth',
					-score => $depth
				);
				$readset_feature->add_SeqFeature($t);

				#}
			}
			my $pivot =
			    $expression_statistics_ref->{'median_hits'}
			  ? $expression_statistics_ref->{'median_hits'} - ( 2 * $hit_sd )
			  : 1;
			my $track = $graph->add_track(
				xyplot       => $readset_feature,    #work!
				-description => 1,
				-label       => 1,
				-max_score     => $hit_max ? $hit_max + ( $hit_max * 0.1 ) : 1,
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
		elsif ($expression_statistics_ref->{'total_hits'}
			&& !$expression_statistics_ref->{'median_hits'}
			&& $expression_statistics_ref->{'total_hits'} > 500 )
		{
			warn Dumper $expression_statistics_ref if $debug;
			warn
"Something really bad happened: maybe the statistics were not properly processed...\n";
		}
		else {
			warn "No coverage of $seq_id ($seq_md5hash) vs $readset !\n"
			  if $debug;
			my $track = $graph->add_track(
				xyplot       => $readset_feature,    #work!
				-description => 1,
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
	print GRAPH $graph->png;
	close GRAPH;
	undef($graph);
}

sub process_main_stats($) {

	# i don't think there is a memory leak here.

	#having both bp data and rpkm, we can post process
	my $seq_id           = shift;
	my $seq_md5hash      = &sqlite_get_md5($seq_id);
	my $seq_size         = &sqlite_get_seq_length($seq_md5hash);
	my $readsets_covered = int(0);
	$demand_all_readsets = scalar @readsets if $demand_all_readsets;
	for ( my $i = 0 ; $i < scalar(@readsets) ; $i++ ) {
		my $readset_metadata_ref =
		  &sqlite_get_readset_metadata( $readsets[$i] );
		my $readset_name        = $readset_metadata_ref->{'alias'};
		my $readset_total_reads = $readset_metadata_ref->{'total_reads'};
		my $hit_stat            = Statistics::Descriptive::Full->new();
		my $stats_ref =
		  &sqlite_get_expression_statistics( $seq_md5hash, $readsets[$i] );
		my $rpkm = $stats_ref->{'rpkm'} ? $stats_ref->{'rpkm'} : int(0);
		my $total_hit_reads =
		  $stats_ref->{'total_hits'} ? $stats_ref->{'total_hits'} : int(0);
		my ( $no_bp_covered, @hits );
		my $depth_hash_ref =
		  &sqlite_get_depth_data( $seq_md5hash, $readsets[$i] );

		if ( !$depth_hash_ref ) {
			warn "No depth data for $seq_id vs $readset_name. Skipping\n"
			  if $debug;
			print LOG "No depth data for $seq_id vs $readset_name. Skipping\n";
			next;
		}

		# parsing depth
		for ( my $pos = 0 ; $pos < $seq_size ; $pos++ ) {
			my $depth =
			    $depth_hash_ref->{$pos}
			  ? $depth_hash_ref->{$pos}
			  : int(0);
			push( @hits, $depth );
			$no_bp_covered++ if $depth == 0;
		}
		$hit_stat->add_data( \@hits );
		my $hit_median = $hit_stat->median() ? $hit_stat->median() : int(0);
		my $hit_sd =
		  $hit_stat->standard_deviation()
		  ? sprintf( "%.2f", $hit_stat->standard_deviation() )
		  : int(0);
		my $hit_mean =
		  $hit_stat->mean() ? sprintf( "%.2f", $hit_stat->mean() ) : int(0);
		my $hit_max = $hit_stat->max() ? $hit_stat->max() : int(0);
		my $mean_reads =
		  $total_hit_reads
		  ? sprintf( "%.2f", ( $total_hit_reads / $seq_size ) )
		  : int(0);
		if ( $hit_median == 0
			|| ( $median_cutoff && $hit_median < $median_cutoff ) )
		{
			print LOG
"Median depth of coverage of $seq_id vs $readset_name didn't pass cutoff or was zero ($hit_median). Skipping\n";
			warn
"Median depth of coverage of $seq_id vs $readset_name didn't pass cutoff or was zero ($hit_median). Skipping\n"
			  if $debug;
			print STATS
"$seq_md5hash\t$readset_name\t$seq_id\t$rpkm\t$hit_mean\t$hit_sd\t$hit_median\t$hit_max\n";
			next;
		}
		$readsets_covered++;
		&sqlite_add_main_expression_statistics(
			$seq_md5hash, $readsets[$i], $hit_mean, $no_bp_covered,
			$mean_reads,  $hit_median,   $hit_max,  $hit_sd
		);

		print STATS
"$seq_md5hash\t$readset_name\t$seq_id\t$rpkm\t$hit_mean\t$hit_sd\t$hit_median\t$hit_max\n";
	}
	if ( $readsets_covered == 0
		|| ( $demand_all_readsets && $readsets_covered < $demand_all_readsets )
	  )
	{
		warn
"Skipping $seq_id as it didn't align to sufficient readsets ($readsets_covered)\n"
		  if $debug;
		print LOG
"Skipping $seq_id as it didn't align to sufficient readsets ($readsets_covered)\n";
		$skipped_references{$seq_md5hash} = 1;
		return;
	}

	&stats_with_graph( $seq_md5hash, $seq_id, $seq_size ) if !$no_graphs;

}

sub process_expression_level() {

	#DEBUG HERE IS A PROBLEM... SOMEWHERE

# given an array of genes with expression levels, sort it, and give me a distribution so i can see how it looks.
	my $expression_level_tsv = "$result_dir/$uid.expression_levels.stats.tsv";
	my $fpkm_expression_level_matrix =
	  "$result_dir/$uid.expression_levels.matrix.fpkm";
	open( OUT,           ">$expression_level_tsv" )           || die($!);
	open( COUNTS_MATRIX, ">$counts_expression_level_matrix" ) || die($!);
	open( FPKM_MATRIX,   ">$fpkm_expression_level_matrix" )   || die($!);
	print OUT
	  "Checksum\tGene alias\tReadset\tRPKM\tFPKM\teff.counts\tKANGADE\n";
	my $count_matrix_print = "#transcript\t";

	for ( my $i = 0 ; $i < scalar(@readsets) ; $i++ ) {
		my $readset_metadata = &sqlite_get_readset_metadata( $readsets[$i] );
		$count_matrix_print .= $readset_metadata->{'alias'} . "\t";
	}
	chop($count_matrix_print);
	my $fpkm_matrix_print = $count_matrix_print;
	print COUNTS_MATRIX $count_matrix_print . "\n";
	print FPKM_MATRIX $fpkm_matrix_print . "\n";
	foreach my $seq_id (@reference_sequence_list) {
		my $seq_md5hash = &sqlite_get_md5($seq_id);
		next if $skipped_references{$seq_md5hash};
		$count_matrix_print = $seq_md5hash . "\t";
		$fpkm_matrix_print  = $seq_md5hash . "\t";
		for ( my $i = 0 ; $i < scalar(@readsets) ; $i++ ) {
			my $readset_metadata_C =
			  &sqlite_get_readset_metadata( $readsets[$i] );
			my $readset_name_C = $readset_metadata_C->{'alias'};
			my $stats_ref_C =
			  &sqlite_get_expression_statistics( $seq_md5hash, $readsets[$i] );
			$stats_ref_C->{'express_eff_counts'} = int(0)
			  if !$stats_ref_C->{'express_eff_counts'};
			$stats_ref_C->{'express_fpkm'} = int(0)
			  if !$stats_ref_C->{'express_fpkm'};
			$stats_ref_C->{'rpkm'} = int(0) if !$stats_ref_C->{'rpkm'};
			$stats_ref_C->{'kangade_counts'} = int(0)
			  if !$stats_ref_C->{'kangade_counts'};
			$count_matrix_print .= $stats_ref_C->{'express_eff_counts'} . "\t";
			$fpkm_matrix_print .=
			  sprintf( "%.2f", $stats_ref_C->{'express_fpkm'} ) . "\t";
			print OUT $seq_md5hash . "\t" 
			  . $seq_id . "\t"
			  . $readset_name_C . "\t"
			  . $stats_ref_C->{'rpkm'} . "\t"
			  . $stats_ref_C->{'express_fpkm'} . "\t"
			  . $stats_ref_C->{'express_eff_counts'} . "\t"
			  . $stats_ref_C->{'kangade_counts'} . "\n";

			# now parse each pair of readsets
			for ( my $k = $i + 1 ; $k < scalar(@readsets) ; $k++ ) {
				next if $readsets[$i] eq $readsets[$k];
				my $readset_metadata_E =
				  &sqlite_get_readset_metadata( $readsets[$k] );
				my $readset_name_E = $readset_metadata_E->{'alias'};
				my $stats_ref_E =
				  &sqlite_get_expression_statistics( $seq_md5hash,
					$readsets[$k] );
				$stats_ref_E->{'express_eff_counts'} = int(0)
				  if !$stats_ref_E->{'express_eff_counts'};
				$stats_ref_E->{'express_fpkm'} = int(0)
				  if !$stats_ref_E->{'express_fpkm'};
				$stats_ref_E->{'rpkm'} = int(0) if !$stats_ref_E->{'rpkm'};
				$stats_ref_E->{'kangade_counts'} = int(0)
				  if !$stats_ref_E->{'kangade_counts'};
				my $raw_rpkm_ratio =
				  !$stats_ref_E->{'rpkm'} || $stats_ref_E->{'rpkm'} == 0
				  ? 'Inf'
				  : sprintf( "%.2f",
					$stats_ref_C->{'rpkm'} / $stats_ref_E->{'rpkm'} );
				$add_raw_fold_change->execute(
					$raw_rpkm_ratio, $seq_md5hash,
					$readsets[$i],   $readsets[$k]
				) unless $contextual_alignment;
				$fold_changes{$seq_md5hash}{'raw'}{$readset_name_C}
				  {$readset_name_E} = $raw_rpkm_ratio;
				my $kangade_ratio =
				  $fold_changes{$seq_md5hash}{'kangade'}{$readset_name_C}
				  {$readset_name_E}
				  ? $fold_changes{$seq_md5hash}{'kangade'}{$readset_name_C}
				  {$readset_name_E}
				  : 'N/A';

				if ($perform_bias_correction) {
					my $fpkm_ratio =
					    !$stats_ref_E->{'express_fpkm'}
					  || $stats_ref_E->{'express_fpkm'} == 0
					  ? 'Inf'
					  : sprintf( "%.2f",
						$stats_ref_C->{'express_fpkm'} /
						  $stats_ref_E->{'express_fpkm'} );
					my $counts_ratio =
					    !$stats_ref_E->{'express_eff_counts'}
					  || $stats_ref_E->{'express_eff_counts'} == 0
					  ? 'Inf'
					  : sprintf( "%.2f",
						$stats_ref_C->{'express_eff_counts'} /
						  $stats_ref_E->{'express_eff_counts'} );
					$add_express_fold_change->execute(
						$fpkm_ratio,   $counts_ratio, $seq_md5hash,
						$readsets[$i], $readsets[$k]
					) unless $contextual_alignment;
					print STATS_RATIO $seq_md5hash . "\t" 
					  . $seq_id . "\t"
					  . $readset_name_C . "\t"
					  . $readset_name_E . "\t"
					  . $raw_rpkm_ratio . "\t"
					  . $kangade_ratio . "\t"
					  . $fpkm_ratio . "\t"
					  . $counts_ratio . "\n";
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
		chop($count_matrix_print);
		chop($fpkm_matrix_print);
		print COUNTS_MATRIX $count_matrix_print . "\n";
		print FPKM_MATRIX $fpkm_matrix_print . "\n";
	}
	close OUT;
	close COUNTS_MATRIX;
	close FPKM_MATRIX;
}

sub prepare_edgeR_graphs() {
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
			&sqlite_set_get_as_housekeeping( $checksum, 1, $readset1,
				$readset2 );
		}
		close IN;
	}
	if (%edgeR_housekeeping_genes) {
		open( OUT1,   ">$unchanged_gene_list" )    || die($!);
		open( OUTALL, ">$housekeeping_gene_list" ) || die($!);
		foreach my $checksum ( sort keys %edgeR_housekeeping_genes ) {
			&sqlite_set_get_as_housekeeping( $checksum, 1 )
			  if $edgeR_housekeeping_genes{$checksum} == $comparisons;
			print OUT1 $checksum . "\t"
			  . $edgeR_housekeeping_genes{$checksum} . "\n";
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
	else {
		&process_edgeR_graphs( \%edgeR_diff_expressed_genes );
		if ( $ps2pdf_exec && $pdfcrop_exec ) {
			my @ps_files = glob( $edgeR_dir . '*ps' );
			foreach my $ps (@ps_files) {
				my $pdf = $ps;
				$pdf =~ s/\.(\w+)$/.pdf/;
				&process_cmd(
					$ps2pdf_exec . "$pdf.t $ps >/dev/null 2>/dev/null" );
				&process_cmd(
					$pdfcrop_exec . " $pdf.t $pdf >/dev/null  2>/dev/null" );
				unlink("$pdf.t");
			}
			if ($inkscape_exec) {
				my @pdf_files = glob( $edgeR_dir . '*.pdf' );
				foreach my $pdf (@pdf_files) {
					my $png = $pdf;
					my $svg = $pdf;
					$svg =~ s/\.(\w+)$/.svg/;
					copy( $pdf, $svg );
					$png =~ s/\.(\w+)$/.png/;

		#inkscape converts input to svg! also does not work with multipage docs.
					&process_cmd(
"$inkscape_exec --vacuum-defs --export-background=white -z -D -d 90 --export-png=$png -f $svg >/dev/null 2>/dev/null"
					) if $inkscape_exec;
				}
			}
		}
	}
	return ( \%edgeR_diff_expressed_genes, \%edgeR_housekeeping_genes );
}

sub process_housekeeping_stats() {

	# skip $skipped_references{$seq_md5hash};
	# %fold_changes

=pod
 ok what do i want here? i want to use one of the matrix files (mine or that after TMM normalization)
 and find which genes are invariable. i think that it would be best to be conservative and do the following:
 
 
 
=cut

}

sub process_completion() {

# N.B. Remember: R/FPKM is for looking at differences between libraries and mean/median_hits for looking at differences within a library.
# RPKM is appropriate for Single-End libraries and FPKM for Paired-End libraries
	my $elapsed = $total_timer->report("%L");
	print "\nCompleted in $elapsed!\n";
	print LOG "\nCompleted in $elapsed!\n";
	print LOG "#Ended: " . &mytime . "\n";
	close LOG;
	&sqlite_destroy();
	## cleanup
	unless ($debug) {
		foreach my $readset (@readsets) {
			my $baseout = $result_dir . $uid . '_' . basename($readset);

			#unlink("$baseout.bam");
			#unlink("$baseout.bam.bai");
		}
	}
}

sub perform_TMM_normalization_edgeR() {
	print "\n" . &mytime
	  . "Performing trimmed mean of fold-change normalization\n";
	my $matrix_file = shift;
	my $matrix_base = basename($matrix_file);
	my $TMM_file    = $edgeR_dir . 'TMM_info.txt';

 # AP: accounting for groups (column 2 [1] in target.dat and TMM_info.txt files)
	my ( $fpkm_TMM_matrix_file, %trans_lengths );
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
			my $dat_file = $edgeR_dir . "$readset_cf.dat";

			$readset_R_filenames{$readset} = $dat_file;
			open( my $ofh, ">$dat_file" )
			  || die( "Cannot write to $dat_file\n" . $! );
			$ofhs{$readset} = $ofh;
		}
		# create a tab delimited file which is gene tab count
		while ( my $ln = <IN> ) {
			chomp($ln);
			my @x = split( "\t", $ln );
			my $gene_name = shift @x;
			foreach my $readset (@readsets_from_matrix) {
				my $ofh = $ofhs{$readset};
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
	$fpkm_TMM_matrix_file =
	  &write_normalized_fpkm_file( $matrix_file, $TMM_file );
	return $fpkm_TMM_matrix_file;
}

sub perform_edgeR_pairwise() {

# this will now estimate differential expression using the counts for each group
	my @groups = sort keys %groups_readsets;
	print &mytime
	  . "Running pair-wise edgeR comparisons for "
	  . scalar( keys %groups_readsets )
	  . " groups...\n";
	foreach my $g (@groups) {
		print "Group $g readsets: "
		  . join( ", ", keys %{ $groups_readsets{$g} } ) . "\n";
	}

	for ( my $i = 0 ; $i < @groups - 1 ; $i++ ) {
		print $groups[$i] . " vs ";
		for ( my $k = $i + 1 ; $k < @groups ; $k++ ) {
			print $groups[$k] . "\n";
			&run_edgeR(
				$groups[$i], $groups[$k]
			);
		}
	}
}

sub process_edgeR_graphs() {

	#from b.haas
	my $edgeR_diff_expressed_genes_ref = shift;
	my $norm_matrix_file = basename($fpkm_expression_level_matrix_TMM);
	chdir($edgeR_dir);
	my $diff_expr_matrix = "$norm_matrix_file.diff_expressed";
	open( OUT, ">$diff_expr_matrix" ) || die($!);
	open( IN,  $norm_matrix_file )    || die($!);
	my $header = <IN>;
	chomp($header);
	print OUT $header . "\tFDR\n";

	while ( my $ln = <IN> ) {
		chomp($ln);
		my @data = split( "\t", $ln );
		if ( $edgeR_diff_expressed_genes_ref->{ $data[0] } ) {
			print OUT $ln . "\t"
			  . $edgeR_diff_expressed_genes_ref->{ $data[0] } . "\n";
		}
	}
	close IN;
	my $R_script = "$diff_expr_matrix.R";
	open( OUTR, ">$R_script" ) or die "Error, cannot write to $R_script. $!";
	print OUTR "
source('$RealBin/R/edgeR_funcs.R')
edgeR_heatmaps_all('$diff_expr_matrix',$tree_clusters)
    ";
	print OUTR "save.image(file='$diff_expr_matrix.Rdata');\n" if $debug;
	close OUTR;
	&process_cmd("R --vanilla --slave -f $R_script >/dev/null 2>/dev/null");

	# we now have the heatmaps and cluster data. Produce JSON objects for web
	# NEWICK data
	my $hcgenes   = $diff_expr_matrix . '.hcgenes.newick';
	my $hcsamples = $diff_expr_matrix . '.hcsamples.newick';
	my $html_file = &prepare_heatmap_for_canvas(
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
	my $htmlfile = $diff_expr_matrix;
	$htmlfile .= '.heatmap.html';
	print &mytime()
	  . "canvas: Preparing HTML/JS code for interactive graph of $htmlfile\n";
	my $js_dir = $edgeR_dir . 'js/';
	mkdir($js_dir) unless -d $js_dir;
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
		die "NEWICK of samples has negative values!\n";
	}
	close IN;
	open( IN, $hcgenes_newick ) || die($!);
	$for_json_array{'t'}{'vars'} = <IN>;
	chomp( $for_json_array{'t'}{'vars'} );
	$for_json_array{'t'}{'vars'} =~ s/\;$//;

	#negative values not allowed!
	if ( $for_json_array{'t'}{'vars'} =~ /-\d\.\d\d+/ ) {
		die "NEWICK of genes has negative values!\n";
	}
	close IN;

	# add var data / ancillary data / properties
	open( IN, $data_file ) || die($!);
	my $data = decode_json(<IN>);
	close IN;
	foreach my $readset ( @{ $for_json_array{'y'}{'smps'} } ) {
		if ( $library_metadata{$readset} && %{ $library_metadata{$readset} } ) {
			foreach my $key ( keys %{ $library_metadata{$readset} } ) {
				push(
					@{ $for_json_array{'x'}{$key} },
					$library_metadata{$readset}{$key}
				);
			}
		}
	}
	for ( my $i = 0 ; $i < @{ $for_json_array{'y'}{'vars'} } ; $i++ ) {
		my $seq_md5sum = $for_json_array{'y'}{'vars'}->[$i];

		#  foreach my $seq_md5sum ( @{ $for_json_array{'y'}{'vars'} } ) {
		my @aliases = &sqlite_get_seq_aliases($seq_md5sum);
		@aliases = qw/none/ unless @aliases;

		#debug trial
		my $alias =
		  $user_alias{$seq_md5sum} ? $user_alias{$seq_md5sum} : $seq_md5sum;
		$for_json_array{'y'}{'vars'}->[$i] = $alias;

# kind of have to do this because R has no access to the aliases and aliases are important as they are the only way to ensure that
# non-identical sequences are used...
		$for_json_array{'t'}{'vars'} =~ s/\b$seq_md5sum\b/$alias/;
		push( @{ $for_json_array{'z'}{'checksum'} }, $seq_md5sum );
		push( @{ $for_json_array{'z'}{'gene'} },     join( ', ', @aliases ) );
		push( @{ $for_json_array{'z'}{'fdr'} },      $fdr{$seq_md5sum} );
		push( @{ $for_json_array{'y'}{'data'} },     $data->{$seq_md5sum} );

# we can add a secondary dataset; dunno how to depict it though if we have a lot of data (or even what to depict); so here it is for the purpose of demonstration
#    foreach my $readset (@{$for_json_array{'y'}{'smps'}}){
#      my $expression_hashref = &sqlite_get_expression_statistics($seq_md5sum,$readset,'alias');
#      # 'readset/file/name' => {'express_fpkm' => '202.231576','hit_sd' => '11.26','max_hits' => 46,'mean_hits' => '16.78','express_eff_counts' => '296.952865','rpkm' => 124,'no_coverage' => 186,'mean_reads' => '0.25','total_hits' => 406,'median_hits' => 18,'kangade_counts' => 381},
#      my $count = $expression_hashref->{'express_fpkm'} ? $expression_hashref->{'express_fpkm'} :$expression_hashref->{'rpkm'} ;
#      $count = int(0) if !$count;
#      push(@{$for_json_array{'y'}{'data2'}},$count);
#    }
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
  <script type="text/javascript"  src="js/canvasXpress.min.js" ></script>
  <script>
    var show_heatmap = function() {
     var cx = new CanvasXpress(
                "canvas",
  
  ' . encode_json( \%for_json_array ) . '
  , {
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
	copy(
		$RealBin . "/js/canvasXpress.min.js",
		$js_dir . '/canvasXpress.min.js'
	);
	return $htmlfile;
}

sub prepare_scatter_for_canvas() {
	my $diff_expr_matrix = shift;

	#2D: Y:logFC; X:logCPM; two series: for FDR <= $fdr_pval_cutoff
	#3D: Y:logFC; X:logCPM; Z:FDR
	my @edgeR_results = glob( $edgeR_dir . "*_vs_*.edgeR.all.txt" );

	# preliminary
	my ( %comparisons, %aliases );
	my $htmlfile_2d = $diff_expr_matrix . '.scatterplots_2D.html';
	my $htmlfile_3d = $diff_expr_matrix . '.scatterplots_3D.html';
	my $js_dir      = $edgeR_dir . 'js/';
	mkdir($js_dir) unless -d $js_dir;
	copy(
		$RealBin . "/js/canvasXpress.min.js",
		$js_dir . '/canvasXpress.min.js'
	);
	my $print_2d = '
  <!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">
  <html>
  <head>
  <title>2D Scatterplots</title>
  <meta http-equiv="X-UA-Compatible"  content="chrome=1">
  <meta http-equiv="CACHE-CONTROL"  CONTENT="NO-CACHE">
  <meta http-equiv="Content-Type"  content="text/html; charset=utf-8"/>
  <script type="text/javascript"  src="js/canvasXpress.min.js" ></script>
  <script> var show_Scatter2D = function() {';
	my $print_3d = '
  <!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">
  <html>
  <head>
  <title>3D Scatterplots</title>
  <meta http-equiv="X-UA-Compatible"  content="chrome=1">
  <meta http-equiv="CACHE-CONTROL"  CONTENT="NO-CACHE">
  <meta http-equiv="Content-Type"  content="text/html; charset=utf-8"/>
  <script type="text/javascript"  src="js/canvasXpress.min.js" ></script>
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
		my $res            = basename($edgeR_result);
		my $single_2d_file = $edgeR_dir . $res . '.scatterplots_2D.html';
		my $single_3d_file = $edgeR_dir . $res . '.scatterplots_3D.html';
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
			my @aliases = &sqlite_get_seq_aliases($seq_md5sum);
			@aliases = qw/none/ unless @aliases;
			my $alias =
			  $user_alias{$seq_md5sum} ? $user_alias{$seq_md5sum} : $seq_md5sum;
			push( @{ $for_json_array{'y'}{'vars'} },     $alias );
			push( @{ $for_json_array{'z'}{'checksum'} }, $seq_md5sum );
			push( @{ $for_json_array{'z'}{'aliases'} },
				join( ", ", @aliases ) );
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

	# get the table of counts, either via express or from survey.
}

=pod 
deseq script from Robles/Alexie
  my $rscript = '
    library(DESeq);
    args <- commandArgs(TRUE);
    InFile  <- args[1];
    OuFile  <- args[2];
    Nreps   <- args[3];
    RsltsPath  <- args[4];
    countsTable <- read.delim(InFile, header=TRUE, stringsAsFactors=TRUE );
    rownames( countsTable ) <- countsTable$X;
    countsTable <- countsTable[ , -1 ];
    #Identify ALL conditions:
    AllConds <- c(rep("ref",Nreps),rep("exp",Nreps));
    cat("Performing DESeq DE Calculation for n =",length(AllConds)/2,"replicates...\n");
    cds <- newCountDataSet( countsTable, AllConds );
    cds <- estimateSizeFactors( cds );
    cds <- estimateDispersions( cds ); #method is blind without reps...
    res <- nbinomTest( cds, "ref","exp" );
    # Write output to a csv file (All isoforms):
    write.table(res[ order(res$pval), ] ,file=sprintf("%s%s",OuFile,"_DESeq.csv"),sep=",");
    ';

=cut

####
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
		print OUT1
		  join( "\t", $R_filename, $library_metadata{$name}{'group'}, $name )
		  . "\n";
	}
	close OUT1;
	my $tmm_norm_script = $edgeR_dir . "tmp_runTMM.R";
	open( OUT2, ">$tmm_norm_script" )
	  or die "Error, cannot write to $tmm_norm_script";
	print OUT2 "
  source('$RealBin/R/edgeR_funcs.R')
  TMM_normalize('$data_file_list','$output')
    ";
	close OUT2;
	&process_cmd( "R --vanilla --slave -f $tmm_norm_script >/dev/null",
		$edgeR_dir );

	unless ( -s $output ) {
		print LOG "TMM normalization failed.\n";
		die "TMM normalization failed.\n";
	}
	unlink($tmm_norm_script) unless $debug;
}

####
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
	my $R_script    = $base . "_run.R";
	my $target_file = $base . "_target.dat";

	open( OUT1, ">$target_file" )
	  or die "Error, cannot write to $target_file: $!";
	print OUT1 "files\tgroup\tdescription\n";
	
	foreach my $group ( $group_A, $group_B ) {
		foreach my $readset ( keys %{$groups_readsets{$group}} ) {
			my $file = $groups_readsets{$group}{$readset} || die "Cannot find file: Group $group; readset $readset\n";
			print OUT1 join( "\t", $file, $group, $readset ) . "\n";
		}
	}
	close OUT1;

	open( OUTR, ">$R_script" ) or die "Error, cannot write $R_script: $!";
	print OUTR "
  source('$RealBin/R/edgeR_funcs.R')
  edgeR_DE_analysis(targetsFile='$target_file',datafile='$base' ,dispersion=$edgeR_dispersion,FDR=$fdr_pval_cutoff,kclusters=$tree_clusters);
  ";
	close OUTR;
	&process_cmd( "R --vanilla --slave -f $R_script >/dev/null 2>/dev/null",
		$edgeR_dir );
	my $html_file = &prepare_heatmap_for_canvas(
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
}

sub write_normalized_fpkm_file {

#init. from b.haas
# this will take the matrix file, from e.g. express or raw, and convert to FPKM by
# applying the TMM normalization previously done with edgeR.
	my ( $matrix_file, $tmm_info_file ) = @_;
	my $normalized_fpkm_file =
	  $edgeR_dir . basename($matrix_file) . ".normalized.FPKM";
	my ( %eff_lib_sizes, %fpkm_hash );

	# get effective library sizes
	open( IN, $tmm_info_file )
	  or die "Error, cannot open file $tmm_info_file: $!";
	my $header = <IN>;
	while ( my $ln = <IN> ) {
		chomp($ln);
		my @data = split( "\t", $ln );
		next unless $data[2] && $data[5];
		$eff_lib_sizes{ $data[2] } = $data[5];
	}
	close IN;

	open( OUT, ">$normalized_fpkm_file" ) || die($!);
	open( IN,  $matrix_file )             || die($!);
	$header = <IN>;
	print OUT $header;
	chomp $header;
	my @readsets_from_matrix = split( "\t", $header );

	while ( my $ln = <IN> ) {
		chomp($ln);
		my @data    = split( "\t", $ln );
		my $gene    = $data[0];
		my $seq_len = &sqlite_get_seq_length($gene)
		  or die "Error, no seq length for $gene";
		my $print = $gene;

		for ( my $i = 1 ; $i < scalar(@data) ; $i++ ) {
			my $readset      = $readsets_from_matrix[$i];
			my $eff_lib_size = $eff_lib_sizes{$readset}
			  or die "Error, no eff lib size for $readset";
			my $read_count = $data[$i];
			my $fpkm =
			  $read_count / ( $seq_len / 1e3 ) / ( $eff_lib_size / 1e6 );
			$fpkm = sprintf( "%.2f", $fpkm );
			$print .= "\t$fpkm";
			$fpkm_hash{$gene}{$readset} = $fpkm;
		}
		print OUT $print . "\n";
	}
	close OUT;

	unless ( $normalized_fpkm_file && -s $normalized_fpkm_file ) {
		warn "EdgeR failed to produce the output file\n";
		warn $normalized_fpkm_file . "\n" if $normalized_fpkm_file;
		return;
	}
	unlink( $matrix_file . '.normalized.FPKM' );
	link( $normalized_fpkm_file, $matrix_file . '.normalized.FPKM' );

	# append FPKM values in the main text outfile and also do ratios
	print "Storing results...\n";
	open( STATSIN, "$result_dir/$uid.expression_levels.stats.tsv" ) || die($!);
	my $header_S = <STATSIN>;
	my @headers = split( "\t", $header_S );
	if ( scalar(@headers) < 8 ) {
		open( STATSOUT, ">$result_dir/$uid.expression_levels.stats.tsv.t" )
		  || die($!);
		chomp($header_S);
		print STATSOUT $header_S . "\tTMM FPKM\n";
		while ( my $ln = <STATSIN> ) {
			chomp($ln);
			my @data = split( "\t", $ln );
			next unless $data[2];
			print STATSOUT $ln . "\t"
			  . $fpkm_hash{ $data[0] }{ $data[2] } . "\n";
		}
		close STATSIN;
		close STATSOUT;
		unlink("$result_dir/$uid.expression_levels.stats.tsv");
		rename(
			"$result_dir/$uid.expression_levels.stats.tsv.t",
			"$result_dir/$uid.expression_levels.stats.tsv"
		);
	}
	open( STATSIN, "$result_dir/$uid.ratio.stats.tsv" ) || die($!);
	my $headerR = <STATSIN>;
	@headers = split( "\t", $headerR );
	if ( scalar(@headers) < 9 ) {
		open( STATSOUT, ">$result_dir/$uid.ratio.stats.tsv.t" ) || die($!);
		chomp($headerR);
		print STATSOUT $headerR . "\tTMM FPKM\n";
		while ( my $ln = <STATSIN> ) {
			chomp($ln);
			my @data = split( "\t", $ln );
			next unless $data[2];
			my $ratio =
			     $fpkm_hash{ $data[0] }{ $data[3] }
			  && $fpkm_hash{ $data[0] }{ $data[3] } > 0
			  ? sprintf( "%.2f",
				$fpkm_hash{ $data[0] }{ $data[2] } /
				  $fpkm_hash{ $data[0] }{ $data[3] } )
			  : int(0);
			print STATSOUT $ln . "\t" . $ratio . "\n";
		}
		close STATSIN;
		close STATSOUT;
		unlink("$result_dir/$uid.ratio.stats.tsv");
		rename(
			"$result_dir/$uid.ratio.stats.tsv.t",
			"$result_dir/$uid.ratio.stats.tsv"
		);
	}

	return $normalized_fpkm_file;
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
		die "$fq: Not in illumina format!\n" unless $ids =~ /^@/;
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
