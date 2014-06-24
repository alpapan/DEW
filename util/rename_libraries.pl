#!/usr/bin/env perl

use strict;
use warnings;

my $csv = shift ||die;
open (IN,$csv);

while (my $ln=<IN>){
	my @data=split("\t",$ln);
	next if !$data[1] || $data[1] eq 'name';
	print "UPDATE readsets SET alias='".$data[1]."' WHERE readset_file='".$data[0]."';\n";
}

close IN;
