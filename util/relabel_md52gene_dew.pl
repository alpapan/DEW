#!/usr/bin/env perl

use strict;
use warnings;

my $md5file = glob("*md5");
die unless $md5file;
my @files_to_edit = @ARGV;
die unless @files_to_edit;

my %hash;
open (IN,$md5file);
while (my $ln=<IN>){
	next if $ln=~/^\s*$/;
	chomp($ln);
	my @data= split("\t",$ln);
	$hash{$data[0]} = $data[1];
}
close IN;

foreach my $file (sort @files_to_edit){
	open (IN,$file);
	open (OUT,">$file.with_ids");
	while (my $ln=<IN>){
		if ($ln!~/^#/ && $ln=~/^(\S+)/){
			my $id = $hash{$1} if $hash{$1};
			if ($id){
				my @data = split("\t",$ln);
				$data[0] = $id;
				$ln = join("\t",@data);
			}
		}
		print OUT $ln;
	}
	close IN;
	close OUT;

}
