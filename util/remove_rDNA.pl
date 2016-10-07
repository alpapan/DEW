#!/usr/bin/env perl

use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use FindBin qw/$RealBin/;

my (@infiles,$is_paired);

my $database = "$RealBin/../databases/rDNA_nt_inv.fsa_nr.dust";

pod2usage $! unless GetOptions(
	'in|fastq:s{1,}' => \@infiles,
	'db|database|rDNA:s' => \$database,
	'paired' => \$is_paired,	# two and only two
);
push(@infiles,@ARGV);
my $kraken_out = "kraken.res";

my $kraken_cmd = "kraken --db $database --threads 4 --fastq-input --preload --only-classified-output --check-names ";
if ($is_paired){
	$kraken_cmd .= " --paired " if $is_paired;
	my $files = join(" ",@infiles);
	my $hash_ref = &kraken_file($files);
	&purge_file($infiles[0],$hash_ref);
	&purge_file($infiles[1],$hash_ref);
}else{
	foreach my $infile (@infiles){
		my $hash_ref = &kraken_file($infile);
		&purge_file($infile,$hash_ref);
	}
}


####################
sub purge_file(){
	my ($infile,$hash_ref) = @_;
	open (IN,$infile);
	open (OUT,">$infile.cut");
	while (my $id = <IN>){
		my $seq = <IN>;my $id2 = <IN>;my $qlt = <IN>;
		if ($id=~/^@(\S+)/){
			my $check = $1;
			chomp($check);
			if ($is_paired){
				chop($check);
				chop($check);
			}
			next if $hash_ref->{$check};
			print OUT $id . $seq . $id2 . $qlt;
		}
	}	
	close IN;
	close OUT;
}



sub kraken_file(){
	my %kraken_hash;
	my $infile = shift;
	system($kraken_cmd." $infile|cut -f 2 > $kraken_out");
	open (IN,$kraken_out);
	while (my $ln=<IN>){
		chomp($ln);
		$kraken_hash{$ln} = 1 if $ln;
	}
	close IN;
	unlink($kraken_out);
	return \%kraken_hash;
}
