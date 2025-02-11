#!/usr/bin/perl -w
# 3bases.pl - 
# Given an input fasta file of DNA sequences, output the sequences as codons sets of 3 nucleotides
# ATGCTAGCATGCATCG -> ATG CTA GCA TGC ATC G
# output filename: inputfile.3bases

$filename = $ARGV[0];
open F, $filename;
@lines = <F>;
close F;

open W, ">".$filename.".3bases"; 
foreach $line (@lines){
  if ($line =~ /^>/){
    print W $line;
  } else {
    my $seq = "";
    for($i=0; $i<length($line)-3;$i+=3){
      $seq .= substr $line, $i, 3;
      $seq .= "{}";
    }
    print W $seq, "\n";
  }
}
close W;
