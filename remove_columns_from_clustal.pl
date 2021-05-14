#!/usr/bin/perl -w

open FILE, $ARGV[0];
@aln = <FILE>;
close FILE;

open FILE, $ARGV[1];
@pos = <FILE>;
close FILE;

$ref= {};
$longest = 1;
$header = "";

@order = ();
foreach $line (@aln){
  if( $line =~ /^CLUSTAL/ ){
    $header = $line;
  }
  next if $line =~ /^CLUSTAL/ || $line =~ /^$/ || $line =~ /^\s/;

  @s = split(/\s+/, $line);
  #chomp $s[1];
  $ref->{$s[0]} .= $s[1];
  if( length($s[0]) > $longest){ 
    $longest = length($s[0]);
  }
  push @order, $s[0];
}

$output = $ARGV[0];
$output =~ s/\.aln//;
$output .= "_trimmed.aln";
open OUT, ">$output";
print OUT $header,"\n\n";
for($j = 0; $j < scalar(@order); $j++){
  $org = $order[$j];
  for($i = scalar(@pos)-1; $i >= 0; $i--){
    next if $pos[$i] =~ '^$';
    substr($ref->{$org}, $pos[$i]-1, 1, "");
  }
  #print join("\t", $org, $ref->{$org}), "\n";
  printf OUT "%-${longest}s", $org;
  print OUT "     ", $ref->{$org}, "\n";
}
close OUT;
