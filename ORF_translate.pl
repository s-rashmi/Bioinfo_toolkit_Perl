# find reading frames and ORF of given DNA
# count the number of codons and check if the codons are complete
# if complete translate

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

# find ORF of DNA, count no. of codons, check if codon complete and translate
@f = (1,2,3);
$count = 0;
%gencode = (
'AUG' => 'M', # Methionine
'UGG' => 'W', # Tryptophan
'UGC' => 'C', 'UGU' => 'C', # Cysteine
'GAC' => 'D', 'GAU' => 'D', # Aspartic Acid
'UUC' => 'F', 'UUU' => 'F', # Phenylalanine
'CAC' => 'H', 'CAU' => 'H', # Histidine
'AAC' => 'N', 'AAU' => 'N', # Asparagine
'UAC' => 'Y', 'UAU' => 'Y', # Tyrosine
'GAA' => 'E', 'GAG' => 'E', # Glutamic Acid
'AAA' => 'K', 'AAG' => 'K', # Lysine
'CAA' => 'Q', 'CAG' => 'Q', # Glutamine
'AUA' => 'I', 'AUC' => 'I', 'AUU' => 'I', # Isoleucine
'UAA' => '*', 'UAG' => '*', 'UGA' => '*', # Stop
'GCA' => 'A', 'GCC' => 'A', 'GCG' => 'A', 'GCU' => 'A', # Alanine
'GGA' => 'G', 'GGC' => 'G', 'GGG' => 'G', 'GGU' => 'G', # Glycine
'CCA' => 'P', 'CCC' => 'P', 'CCG' => 'P', 'CCU' => 'P', # Proline
'ACA' => 'T', 'ACC' => 'T', 'ACG' => 'T', 'ACU' => 'T', #Threonine
'GUA' => 'V', 'GUC' => 'V', 'GUG' => 'V', 'GUU' => 'V', # Valine
'CUA' => 'L', 'CUC' => 'L', 'CUG' => 'L', 'CUU' => 'L', 'UUA' => 'L', 'UUG' => 'L', # Leucine
'CGA' => 'R', 'CGC' => 'R', 'CGG' => 'R', 'CGU' => 'R', 'AGA' => 'R', 'AGG' => 'R', # Arginine
'UCA' => 'S', 'UCC' => 'S', 'UCG' => 'S', 'UCU' => 'S', 'AGC' => 'S', 'AGU' => 'S' # Serine
);

foreach $frame(@f){
	print "\nReading frame: $frame\n";
	$count = 0;

	for($i = $frame-1; $i < length($dna); $i=$i+3){
		$sub = substr($dna, $i, 3);
		print "$sub ";
		$count++ if length($sub) == 3;
	}
	print "\n";
	print "No. of codons: $count\n";

	print "\nReading frame: -$frame\n";
	$count = 0;

	for($i = $frame-1; $i <length($dna); $i=$i+3){
		$sub = substr($revC, $i, 3);
		print "$sub ";
		$count++ if length($sub) == 3;
	}
	
	print "\n";
	print "No. of codons: $count\n";
}

$start = index($dna, "ATG");
$orf = "ATG ";
$pseq = "M";

for($i = $start+3; $i < length($dna); $i=$i+3){
	$sub = substr($dna, $i, 3);
	if($sub eq "TAA" || $sub eq "TAG" || $sub eq "TGA"){
		$orf .= $sub;
		$sub =~ tr/T/U/;
		$pseq .= $gencode{$sub};
		last;
	}
	
	$orf .= $sub." ";
	$sub =~ tr/T/U/;
	$pseq .= $gencode{$sub};
}

print "\nORF:\n$orf\n";
print "\nProtein:\n$pseq\n";
