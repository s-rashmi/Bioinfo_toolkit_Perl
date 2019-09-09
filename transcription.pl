# transcription
# transcribe given DNA sequence into RNA

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

# translation
$rna = $dna;
$rna =~ tr/'T'/'U'/;
print "\nRNA for given seq:\n$rna\n";
