#!/usr/bin/perl -w

use strict;
use File::Basename;

my @files = @ARGV;
my $md5sum = shift(@files);
die unless @files && -s $md5sum;

my %hash;
open( IN, $md5sum );
while ( my $ln = <IN> ) {
  chomp($ln);
  my @data = split( "\t", $ln );
  next unless $data[1];
  $hash{$data[0]} = $data[1] if $ln;
}
close IN;

foreach my $file (@files){
  next unless -s $file;
  next if -d $file;
  
  my $filename = basename($file);
  my $dirname = dirname($file);
  $filename =~/^(\w+)(.+)$/;
  my $md5 = $1 ||next;
  my $xtn = $2;
  next unless $hash{$md5};
  mkdir('gene_names') unless -d 'gene_names';
  my $new_filename = $dirname .'/gene_names/'. $hash{$md5}.$xtn;
  symlink('../'.$file,$new_filename);
}