#!/usr/bin/perl -w

=pod


=head1 AUTHORS

 Alexie Papanicolaou

        CSIRO Ecosystem Sciences
        alexie@butterflybase.org

=head1 DISCLAIMER & LICENSE

Copyright 2012-2014 the Commonwealth Scientific and Industrial Research Organization. 
This software is released under the Mozilla Public License v.2.

It is provided "as is" without warranty of any kind.
You can find the terms and conditions at http://www.mozilla.org/MPL/2.0.

=cut


use strict;
use FindBin qw/$RealBin/;
use lib $RealBin. '/PerlLib';
use Digest::SHA qw/sha1/;
use Fasta_reader;

my ( %hash, %redundant );
my $count = int(0);
my $file = shift || die;
open( OUT, ">$file.noredundant.fsa" );

my $in_obj = new Fasta_reader($file);
while ( my $object = $in_obj->next() ) {
 my $seq = $object->get_sequence();
 my $id  = $object->get_accession();
 my $md5 = sha1($seq);
 if ( !$hash{$md5} ) {
  print OUT ">$id\n" . $seq . "\n";
  $hash{$md5}++;
 }else {
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
print "Excluded $count identically redundant sequences from $file\n";