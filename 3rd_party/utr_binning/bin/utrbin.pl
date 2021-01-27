#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use Cwd;
use Pod::Usage;
use Getopt::Long;
use File::Basename;
use FindBin qw/$RealBin/;
$ENV{'PATH'} = "$RealBin:" . $ENV{'PATH'};

=head1 NAME

 utrbin.pl - Process UTR RNA-Seq data (e.g. Lexogen QuantSeq 3' mRNA-Seq)

=head1 USAGE

 Mandatory
	-gff_file	:s	GFF3 file with gene annotations in relation to genome
	-bam_files	:s{1,}	1+ BAM files with alignments of reads to genome. Must be co-ordinate sorted (SO:coordinate)
	-gene_fasta	:s	FASTA file with gene FASTA (cf. create_features_from_gff3.pl or other program)
	
 Optional
	-min_aln_length :i	Minimum read to genome alignment to consider (incl indels; defaults to 40)
	-max_aln_length :i	Maximum read to genome alignment to consider (incl indels; defaults to 80)
	-skip_prop      :f      0-1 float: Skip reads whose sequence contains a single nucleotide more than -skip_prop times; defaults to 0.7. Set to 0 to disable
	-limit_end      :i      Maximum number of bases from gene end for the (left-most) alignment to be considered (default 0, ie. not used)
				For Lexogen that would be ~ 700 depending on 3'UTR annotation quality
	-utr		:i	Which UTR to expect data on (5 or 3; defaults to 3)
	-direction	:s	Forward or Reverse (case insensitive or just f or r): 
				only keep reads that are in Forward or Reverse in relation to gene. Defaults to none, keep both
				For Lexogen quantseq 3mrna sequencing (FWD) that would be Forward
	-cpus|threads   :i	Number of parallel threads to use for samtools etc (defaults to 4)

=head1 AUTHORS

Sam El-Kamand and Alexie Papanicolaou

        Western Sydney University
        alexie@butterflybase.org

=head1 DISCLAIMER & LICENSE

Copyright 2021 the Western Sydney University
This software is released under the Mozilla Public License v.2.

It is provided "as is" without warranty of any kind.
You can find the terms and conditions at http://www.mozilla.org/MPL/2.0.


=head1 BUGS & LIMITATIONS

None known

=cut

my ($verbose);
my $cwd = getcwd;

my ($samtools_exec,$bedtools_exec,$gff2bed_exec,$Rscript_exec,$sort_exec ) = &check_program('samtools','bedtools','gff2bed','Rscript','sort');
my ($gff_file,@bam_files,$gene_fasta, %gff_index,$limit_end,$direction);
my $bedtool_direction="";
my ($min_aln_length,$max_aln_length,$utr_use,$skip_prop,$cpus) = (40,80,3,0.7,4);

pod2usage $! unless GetOptions(
	'gff_file=s' => \$gff_file,
	'bam_files=s{1,}' => \@bam_files,
	'gene_fasta=s' => \$gene_fasta,
	'min_aln_length:i' => \$min_aln_length,
	'max_aln_length:i' => \$max_aln_length, 
	'debug|verbose'	   => \$verbose,
	'limit_end:i'      => \$limit_end,
	'utr:i'		   => \$utr_use,
	'direction:s'      => \$direction,
	'skip_prop:f'	   => \$skip_prop,
	'cpus|threads:i'   => \$cpus
);

pod2usage "No annotation GFF file provided or not found" if !$gff_file || !-s $gff_file;
pod2usage "No alignment BAM files provided or not found" if !$bam_files[0] || !-s $bam_files[0];
pod2usage "No gene FASTA file provided or not found" if !$gene_fasta || !-s $gene_fasta;
&check_options();

&index_gff();
#Convert gff3 file to bed file
&process_cmd("$gff2bed_exec < $gff_file | "
              . 'awk \'{if($8=="gene") print }\' '
              . " | $sort_exec -k1,1 -k 2,2n  > $gff_file.genes.bed");
my $bedgenecount = `wc -l < $gff_file.genes.bed`;chomp($bedgenecount);
my $new_reference = &convert_uc($gene_fasta);

foreach my $bam_file (@bam_files){ 
	my $tsv_file = &process_bam($bam_file);
	next if -s "$tsv_file.bam" && !$verbose;
	print "Creating new alignment from $tsv_file\n";
	open (IN,$tsv_file);
	my %tsv_hash;
	my @headers = split("\t",<IN>);
	while (my $ln=<IN>){
		my @data = split("\t",$ln);
		chomp($data[-1]);
		next unless $data[4];
		my @reads = split(',',$data[0]);
		foreach my $read (@reads){
			$tsv_hash{$read} = $data[4];
		}
	}
	close IN;
	my @sam_headers = `$samtools_exec view -H $bam_file|grep -v '^\@SQ'`;
	open (BAM,"$samtools_exec view $bam_file |");
	open (OUT,">$tsv_file.sam");
	print OUT @sam_headers;
	while(my $record = <BAM>){
		next if $record=~/^@/ || $record =~/^\s*$/;
	        my @data = split("\t",$record);

		# check complexity
		next if $skip_prop && $skip_prop >0 && &estimate_complexity($data[9]) == int(0);

		if ($tsv_hash{$data[0]}){
			my $gene_id = $tsv_hash{$data[0]};
			die "Unknown gene $gene_id" unless $gff_index{$gene_id}{'strand'};
			my $gene_length =  $gff_index{$gene_id}{'end'} - $gff_index{$gene_id}{'start'};
			my $aln_length = &cigar_to_alnlen($data[5]);

			$data[2] = $gene_id;
			# new start; 1- based leftmost mapping POSition
			# calculate how many bases from the left of the gene, when orientated 5'->3' (even if antisense)
			if ($gff_index{$gene_id}{'strand'} eq '-'){
				$data[1] = $data[1] & 16 ? 0 : 16;
				my $gene_end = $gff_index{$gene_id}{'end'}; # rightmost bp of gene relation to genome, i.e. gene start if strand eq -
				$data[3] = $data[3] < $gene_end ? $gene_end - $data[3] : int(0);
			}else{
				my $gene_start = $gff_index{$gene_id}{'start'};
				$data[3] = $gene_start < $data[3] ? $data[3] - $gene_start : int(0);
			}

			# alignment was outside gene start
			next if $data[3] == 0;
			# if read alignment is > gene length, then we have to skip it too
			next if ( $data[3] + $aln_length > $gene_length);

			# if alignment > $limit_end then skip
			if ($utr_use && $limit_end && $limit_end >0){
  			  if ($utr_use == 3){
				next unless ($data[3] > ($gene_length - $limit_end));
			  }elsif($utr_use == 5){
				next unless ($data[3] < $limit_end);
			  }
			}

			print OUT join("\t",@data);
		}
		# else we are not printing it
	}
	close BAM;
	close OUT;
	&process_cmd("$samtools_exec view -u -T $new_reference $tsv_file.sam 2>/dev/null| $samtools_exec calmd -u -A /dev/stdin $new_reference 2>/dev/null | $samtools_exec sort -l 9 --threads $cpus -m 1G -o $tsv_file.bam && $samtools_exec index $tsv_file.bam");
	unlink("$tsv_file.sam") unless $verbose;
}



#################################
sub check_program() {
 my @paths;
 foreach my $prog (@_) {
  my $path = `which $prog`;
  die "Error, path to required $prog cannot be found\n"
    unless $path =~ /^\//;
  chomp($path);
  push( @paths, $path );
 }
 return @paths;
}

sub process_cmd {
 my ( $cmd, $dir ) = @_;
 print "CMD: $cmd\n" if $verbose;
 chdir($dir) if $dir;
 my $ret = system($cmd);
 if ( $ret && $ret != 256 ) {
  die  "Error, cmd died with ret $ret\n";
 }
 chdir($cwd) if $dir;
 return;
}


#################################
sub process_bam(){
	my $bam_file = shift || return;
	print "Processing $bam_file\n";
	my $base_name = fileparse($bam_file,('.bam'));
	my $tsv_file = $base_name.".ambiguity_resolved.tsv";
	return $tsv_file if -s $tsv_file;
	# convert bam to bed file with cigar string in 7th column. Note score is the mapping quality but we could change that if we want
	# Make sure bed is sorted by chromosome and start pos
	#filter out long joints (Cigar string, N > X) #(fix to use shell variables)
	#For the moment just filtering by total length of alignment (<= MAX_ALIGNMENT_LENGTH and >= MIN_ALIGNMENT_LENGTH), note read length
	&process_cmd("$bedtools_exec bamtobed -i $bam_file -cigar| $sort_exec -k1,1 -k2,2n > $base_name.bed");

	# ALEXIE TODO samtools view -F 4 -m $minquery -q 40 --threads $cpus


	die "Conversion to BED failed\n" unless -s "$base_name.bed";

	&process_cmd("awk -v max='$max_aln_length' -v min='$min_aln_length'"
                     .' \'{if($3-$2 <= max && $3-$2 >= min) print}\' '
                     . "< $base_name.bed > $base_name.filtered.bed");

	my $prefiltercount = `wc -l < $base_name.bed`;chomp($prefiltercount);
	my $postfiltercount = `wc -l < $base_name.filtered.bed`;chomp($postfiltercount);
	my $readsfiltered = $prefiltercount-$postfiltercount;
	my $percentsurvived = sprintf( "%.2f",$postfiltercount/$prefiltercount);
	#merge overlapping reads (on strand according to -direction). Then re-arrange columns to ensure it stays in bed format. Note that score column contains the number of reads merged to create the interval

	#SAM TODO: bedtools merge -c 1,1 -o collapse,count -s -d -1 -i $BAM > merged.bed
	&process_cmd("$bedtools_exec merge $bedtool_direction -d -1 -i $base_name.filtered.bed -c 6,4,4, -o distinct,collapse,count | "
		     . ' awk \'{OFS="\t"; print $1,$2,$3,$5,$6,$4}\' '
		     . " > $base_name.filtered.merged.bed");
	my $mergedintervalcount=`wc -l < $base_name.filtered.merged.bed`;chomp($mergedintervalcount);
	
	#intersect read data with genes
	&process_cmd("$bedtools_exec intersect $bedtool_direction -wa -wb -a $base_name.filtered.merged.bed -b $gff_file.genes.bed > $base_name.intersecting.genes.bed");
	&process_cmd("cat $base_name.intersecting.genes.bed | "
                    .'awk \'{print $1,$2,$3,$10,$4}\' '
		    ." | $sort_exec -k5,5 | uniq -d -f4 > $base_name.ambiguous.genes.tsv");
	my $geneoverlappedbymergedalignment = `wc -l < $base_name.intersecting.genes.bed`;chomp($geneoverlappedbymergedalignment);
	
#	&process_cmd("$Rscript_exec --vanilla $RealBin/removeDuplicateIntervals.R $base_name.intersecting.genes.bed $base_name.ambiguity_resolved");
	&replace_R_script($base_name);
	unless ($verbose){
		unlink("$base_name.bed");
		unlink("$base_name.filtered.bed");
		unlink("$base_name.intersecting.genes.bed");
		unlink("$base_name.filtered.merged.bed");
		unlink("$base_name.ambiguous.genes.tsv");
		unlink("$base_name.ambiguity_resolved.bed");
	}

	warn "Procedure failed for $base_name. Skipping\n" unless -s $tsv_file;
	return $tsv_file;
}

################################
sub index_gff(){
print "Indexing GFF $gff_file\n";
open (GFF,$gff_file);
while (my $ln=<GFF>){
	my @data = split("\t",$ln);
	next unless $data[2] && $data[2] eq 'gene';
	chomp($data[-1]);
	my $id = $data[8];
	$id=~s/^.*ID=([^;]+).+$/$1/;
	$gff_index{$id} = {
		'reference' => $data[0],
		'start'     => $data[3], # 1- based leftmost mapping POSition
		'end'	    => $data[4],
		'strand'    => $data[6]
	};
}
close GFF;
}


########################
sub convert_uc(){
  # this is needed because samtools view breaks with lowercase refences
  my $infile = shift || die "No file";
  my $outfile = $infile.'.uc';
  return $outfile if -s $outfile;
  open (IN,$infile);
  open (OUT,">$outfile");

  while (my $line=<IN>){
	if ($line=~/^>/){print OUT $line;}
	else {print OUT uc($line);}
 }
 close (IN);
 close (OUT);
 return $outfile;
}

###########################
sub cigar_to_alnlen(){
#M   alignment match (can be a sequence match or mismatch)
#I   insertion to the reference
#D   deletion from the reference
#N   skipped region from the reference
#S   soft clipping (clipped sequences present in SEQ)
#H   hard clipping (clipped sequences NOT present in SEQ)
#P   padding (silent deletion from padded reference)
#=   sequence match
#X   sequence mismatch

# 72M4S
	my $cigar = shift;
	my $aln_length = int(0);
	$aln_length += $1 while ($cigar=~/(\d+)M/g);
	$aln_length += $1 while ($cigar=~/(\d+)I/g);
	$aln_length -= $1 while ($cigar=~/(\d+)D/g);
	$aln_length -= $1 while ($cigar=~/(\d+)S/g);
	
	return $aln_length;

}


sub replace_R_script(){
	my $basefile = shift || return;
	my $file = $basefile .'.intersecting.genes.bed';
	my $outfile = $basefile .'.ambiguity_resolved';

	open (IN,$file);
	
	my %hash;
	while (my $ln=<IN>){
#naturalis.1007  16556   16632   NS500799:498:H5LWNBGXB:3:11610:22454:15542      1       +       naturalis.1007  12060   17204   JAMg_model_1084 .       +       JAMg    gene    .       ID=JAMg_model_1084;Name=JAMg_model_1084
#naturalis.1     2582276 2582347 NS500799:498:H5LWNBGXB:4:13401:23705:12897      1       -       naturalis.1     2582331 2583350 JAMg_model_773  .       -       JAMg    gene    .       ID=JAMg_model_773;Name=JAMg_model_773
		chomp($ln);
		next unless $ln;
		my @data = split("\t",$ln);
		next unless $data[11];
		next unless $data[0] eq $data[6]; # just in case
		#TODO Add a check on strand - user defined.
 		if ($direction){
			next if $direction eq 'f' && $data[5] ne $data[11];
			next if $direction eq 'r' && $data[5] eq $data[11];
		}

		my @reads = split(",",$data[3]);
		foreach my $read (@reads){

		  #skipping if found in more than one good gene
		  next if $hash{$read} && $hash{$read}{'gene'} ne $data[9];
		  #proceed. other checks on which UTR and region are dealt with in the BAM processing ($limit_end)

		  if (  !$hash{$read} || (
			$hash{$read} && 
			  ($hash{$read}{'genomestrand'} eq '+' && $hash{$read}{'leftmostgenealign'} > $data[7])
			||($hash{$read}{'genomestrand'} eq '-' && $hash{$read}{'leftmostgenealign'} < $data[7])
			)){

		  $hash{$read} = {
			'chrom' => $data[0],
			'start' => $data[1],
			'end' => $data[2],
			'gene' => $data[9],
			'reads_mapped' => $data[4],
			'genomestrand' => $data[5],
			'leftmostgenealign' => $data[7]
		  };
		 }
		}
	}
	close IN;
        open (OUT1,">$outfile.tsv"); #read_names chrom start end allocated_feature_name reads_mapped strand
	print OUT1 "read_name\tchrom\tstart\tend\tgene\tgenomestrand\n";
	foreach my $read (keys %hash){
		print OUT1 "$read\t"
			   .$hash{$read}{'chrom'}."\t"
			   .$hash{$read}{'start'}."\t"
			   .$hash{$read}{'end'}."\t"
			   .$hash{$read}{'gene'}."\t"
			   .$hash{$read}{'genomestrand'}."\n";
	}
	close OUT1;
}
#############################
sub check_options(){
	pod2usage "-skip_prop needs to between 0 and 1" unless $skip_prop >= 0 && $skip_prop <= 1;
	foreach my $f (@bam_files){ die "BAM file $f not found\n" unless -s $f; }
	if ($direction){
		if ($direction=~/^f/i){
			$direction = 'f';
			 $bedtool_direction = ' -s ';
		}elsif($direction=~/^r/i){
			$direction = 'r';
			$bedtool_direction = ' -S ';
		}else{
			pod2usage "Direction can only be forward or reverse\n"
		}
	}
	pod2usage "-utr can only be 3 or 5\n" unless !$utr_use || $utr_use == 5 || $utr_use == 3;
	if ($cpus && $cpus > 1){
		$sort_exec .= " --parallel=$cpus ";
	}
}


#############################
sub estimate_complexity(){
	# return 1 (pass) or 0 (not pass)
	my $sequence = shift;
	return 1 if $skip_prop == 0;
	$sequence = uc($sequence);
	$sequence=~s/[^ATCG]+//g;
	my $len = length($sequence);
	my $As = ( $sequence =~ tr/A// );
	my $Ts = ( $sequence =~ tr/T// );
	my $Cs = ( $sequence =~ tr/C// );
	my $Gs = ( $sequence =~ tr/G// );

	return int(0) if  ($As/$len) > $skip_prop;
	return int(0) if  ($Ts/$len) > $skip_prop;
	return int(0) if  ($Cs/$len) > $skip_prop;
	return int(0) if  ($Gs/$len) > $skip_prop;
	return 1;
}
