# Given a string of DNA sequence
# 1) check whether given sequence is a DNA
# 2) count the nucleotide frequency
# 3) Calculate GC & AT content

#!/usr/bin/perl

print "Enter DNA seq:\n";
chomp($dna = <>);
$dna = uc($dna);
@dna = split(//, $dna);

# verify whether DNA
if( $dna =~ /[^ATGC]/){
	print "\nGiven sequence is NOT DNA\n";
	exit;
} else {
	print "\nGiven sequence is verified as DNA\n";
}

# nucleotide freq
$a = $dna =~ tr/'A'/'A'/;
$t = $dna =~ tr/'T'/'T'/;
$g = $dna =~ tr/'G'/'G'/;
$c = $dna =~ tr/'C'/'C'/;

print "\nNucleotide frequency:\n";
print "A:$a\n";
print "T:$t\n";
print "G:$g\n";
print "C:$c\n";

# GC content/ AT content
$len = length($dna);
$gc = (($g + $c)/$len)*100;
print "\nGC content:\n$gc\n";

$at = (($a + $t)/$len)*100;
print "\nAT content:\n$at\n";
