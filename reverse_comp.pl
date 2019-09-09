# Print the reverse complement of the given DNA

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

# reverse complement
$cmp = $dna;
$cmp =~ tr/ATGC/TACG/;
$revC = reverse($cmp);
print "\nReverse complement:\n";
print "$revC\n";
