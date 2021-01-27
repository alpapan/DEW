#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;

# merge
#T9Leaf1_S53_L001_R1_001.fastq.trimmomatic
#T9Leaf1_S53_L001_R2_001.fastq.trimmomatic
#T9Leaf1_S53_L002_R1_001.fastq.trimmomatic
#T9Leaf1_S53_L002_R2_001.fastq.trimmomatic


my $suffix = shift;
$suffix = 'fastq' unless $suffix;
$suffix =~s/^\.+//;
my @files = glob("./*$suffix");

my %hash;

for (my $i=0;$i<scalar(@files);$i++){
	my $file = $files[$i];	
	if ($file=~/^(\S+)_S\d+_L(\d+)_R(\d)_(.+)$/){
		my ($prename,$lane,$direction,$ending) = ($1,$2,$3,$4);
		$hash{$prename}{$direction}{$lane} = $file;
	}
}

foreach my $prename (sort keys %hash){
	foreach my $direction (sort keys %{$hash{$prename}}){
		my $move_str = join(" ",sort (values(%{$hash{$prename}{$direction}})));
		if (scalar(keys %{$hash{$prename}{$direction}}) > 1){
			my $cmd = "cat $move_str "
				. " > $prename"
				.'_R'
				.$direction
				.'.'
				.$suffix ;
			print $cmd."\n";
			sleep(1);
			system ($cmd);
		}else{
			my $cmd = "mv $move_str"
				." "
				.$prename
				.'_R'
				.$direction
				.'.'
				.$suffix ;
			print $cmd."\n";
			sleep(1);
			system ($cmd);
		}
	}

}
