#!/usr/bin/perl -w
use strict;
my $md5sum = shift;
my $txtfile = shift || die;
my %hash;
open( IN, $txtfile );
while ( my $ln = <IN> ) {
  chomp($ln);
  $hash{$ln} = 1 if $ln;
}
close IN;
open( IN,  $md5sum );
open( OUT, ">$txtfile.genes" );
while ( my $ln = <IN> ) {
  my @data = split( "\t", $ln );
  next unless $data[1];
  print OUT $data[1] if $hash{ $data[0] };
}
close IN;
close OUT;
