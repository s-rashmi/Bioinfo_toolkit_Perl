#!/usr/bin/perl -w

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
