#!/usr/bin/perl

package rashmi; # using Perl v5.10.1

use 5; # requires atleast Perl 5 interpreter
use strict; # restricted - global vars referred by full package name
use warnings; # warnings for potential errors

require Exporter;
# allows to export names into caller's namespace

our(@ISA, @EXPORT, @EXPORT_OK);
# predeclares as global variables, to be used under 'strict'. Another way is to use vars, but is somewhat deprecated.

@ISA = qw(Exporter);
# if a method is not found in this package, Perl searches for it in the packages in the @ISA array

# Symbols, functions, variables, etc. to be exported by default
@EXPORT = qw(mysql_connect hyphy_basepath); 

# Symbols, functions, variables, etc. to be exported on request
@EXPORT_OK = qw(single_breakpoint run_gard run_gard_processor parse_gard run_NielsenYang run_YangNielsen_bs run_phyml rm_bootstrap full_tree_text tree_text_wo_bootstrap tree2num);

# Other modules
use File::Temp;
use Cwd;
use Slctree;
use slchen;

#--------------------------------------------------------------
# Sub-routines
# for arrays and hashes, pass them to the subroutine as references

# return my database credentials
sub mysql_connect {
  my $database = $_[0];
  if (!defined $database) { $database = "rashmi_test"; } 
  my $host = 'slchen.gis.a-star.edu.sg';
  my $user = 'rashmi';
  my $pass = 'rashmi';
  my $dbh = DBI->connect('DBI:mysql:database='.$database.';host='.$host, $user, $pass, { LongReadLen => 100000000 });

  return $dbh;
}
# END MYSQL_CONNECT


# get basepath of HYPHY
sub hyphy_basepath {
  my $hostname = `hostname`;
  chomp $hostname;
  my $basepath = "";
  if( $hostname =~ /compute-/){
    $basepath = "/home/gisv88/bin/hyphy/";
  } else {
    #my @location = `locate HYPHY`;
    #defined $location[0] 
    #or die "Can't locate HYPHY. Make sure HYPHY is installed!\n";
    #chomp $location[0];
    #$basepath = $location[0]."/";
    $basepath = "/usr/local/src/hyphy_2.0/hyphy/";
  }
  return $basepath;
}
# END HYPHY_LIBPATH



# given a reference to a hash, get variables to be passed to Hyphy-SBP
# run sbp and return the results
# Usage: single_breakpoint( $reference_to_seqhash, ref to options hash ); scalar value
sub single_breakpoint {
  my $basepath = rashmi::hyphy_basepath;
  my $batchfile = $basepath."lib/hyphy/TemplateBatchFiles/SingleBreakpointRecomb.bf";

  my $tempdir = File::Temp::tempdir(CLEANUP => 1);
  my $fasta = $tempdir."/fasta.temp";
  my $fasta_path = Cwd::abs_path($fasta);
  my $results = $tempdir."/results.temp";
  my $results_path = Cwd::abs_path($results);

  my $href = $_[0];
  my $oref = $_[1];
  my $rref = {};

  my $name_hash = slchen::to_phylip_names($href);
  my $sbp_run = "";

  # define defaults
  if ( !defined $href ){
    die "No alignment provided! Aborting!\n"; 
  } else {
    open (FF, ">".$fasta);
    foreach my $key (sort keys %$href){
      print FF ">$key\n";
      print FF $href->{$key}, "\n";
    }
    close FF;
  }

  if ( !defined $oref->{datatype} ) {
    $oref->{datatype} = "1"; } # Nucleotide
  if ( !defined $oref->{KH}->{testing} ) {
    #$oref->{KH}->{testing} = "3"; }
    # Verify with KH resampling, using joint tree as the null 
    $oref->{KH}->{testing} = "1"; }
    # use only AIC to measure 

  if ( $oref->{KH}->{testing} >= "2" ) { # verify with KH resampling
    if( !defined $oref->{KH}->{sampling} ) {
      $oref->{KH}->{sampling} = "50"; } # this should be reciprocal of p-value desired, thus for p-value=0.02, we need 50 replicates
  }
  if ( !defined $oref->{model}->{name} ) {
    $oref->{model}->{name} = "HKY85"; }
  if ( !defined $oref->{model}->{options}->[0] ) { 
    # must check element - just checking for array returns defined!
    $oref->{model}->{options} = [2]; } # Global parameters

  my $kernal = $basepath."HYPHYMP";
  if( $oref->{KH}->{testing} >= 2 ){
    $sbp_run = `(echo $oref->{datatype}; echo $fasta_path; echo $oref->{KH}->{testing}; echo $oref->{KH}->{sampling}; echo $oref->{model}->{name}; echo $oref->{model}->{options}->[0]; echo $results_path)|'$kernal' $batchfile`;
  } else {
    $sbp_run = `(echo $oref->{datatype}; echo $fasta_path; echo $oref->{KH}->{testing}; echo $oref->{model}->{name}; echo $oref->{model}->{options}->[0]; echo $results_path)|'$kernal' $batchfile`;
  }

  slchen::from_phylip_names($href, $name_hash);
  $sbp_run = &slchen::from_phylip_names($sbp_run, $name_hash);

  # check if SBP reports recombination
  if( $sbp_run =~ /\s*There are too few potential break points to support 1 recombination events\./ ){
    $rref->{result} = "NO";
    $rref->{output} = "There are too few potential break points to support 1 recombination events."; 
  } elsif( $sbp_run =~ /Save results to::/) { 
    my @text = split /Save results to::/, $sbp_run;
    if($text[1] ne "" && $text[1] =~ m/AIC-c\s+Best supported breakpoint is located at position/){
      $rref->{result} = "YES";
      $rref->{output} = $text[1];
      my $splits_file = $results_path."_cAIC.splits";
      my $splits_path = Cwd::abs_path($splits_file);
      my $read_splits = do {
        open (my $sf, $splits_path) or die $!;
        local $/ = undef;
        <$sf>;
      };
   
      $read_splits = &slchen::from_phylip_names($read_splits, $name_hash);
      $rref->{splits} = $read_splits; 
    } else {
      $rref->{result} = "NO";
      $rref->{output} = $text[1];
    }
  
    # parse for text to be displayed on web 
    #my @lines = split("\n\n", $text[1]);
    #my $parsed_text = "";
    #foreach my $l (0..$#lines){
    #  if( $lines[$l] !~ (/Breakpoint at position\s*|messages/)) { 
    #      $parsed_text .= $lines[$l]."\n"; 
    #  } 
    #}
    #$rref->{parsed} = $parsed_text;

  } else {
    $rref->{result} = "NO";
    $rref->{output} = "SBP did not run successfully!"; 
  }

  system "rm -rf $tempdir";
  return $rref;
}
# END SINGLE_BREAKPOINT



# given a reference to a hash, get variables to be passed to Hyphy-Gard
# run GARD and return the results
# Usage: run_gard( $reference_to_seqhash, ref to options hash ); scalar value
sub run_gard {
  my $basepath = rashmi::hyphy_basepath;
  my $batchfile = $basepath."lib/hyphy/TemplateBatchFiles/GARD.bf";

  my $tempdir = File::Temp::tempdir(CLEANUP => 1);
  my $fasta = $tempdir."/fasta.temp";
  my $fasta_path = Cwd::abs_path($fasta);
  my $results = $tempdir."/results.temp";
  my $results_path = Cwd::abs_path($results);

  my $href = $_[0];
  my $oref = $_[1];
  my $rref = {};

  my $name_hash = slchen::to_phylip_names($href);

  # define defaults
  if ( !defined $href ){
    die "No alignment provided! Aborting!\n";
  } else {
    open (FF, ">".$fasta);
    foreach my $key (sort keys %$href){
      print FF ">$key\n";
      print FF $href->{$key}, "\n";
    }
    close FF;
  }

  if ( !defined $oref->{model}->{char} ) {
    $oref->{model}->{char} = "010010"; } # for HKY85
  if ( !defined $oref->{rateVariation}->[0] ) {
    # must check element - just checking for array returns defined!
    $oref->{rateVariation} = [1]; } # Homogeneous rates across sites 
  if ( $oref->{rateVariation}->[0] ne "1" ) {
    if ( !defined $oref->{rateVariation}->[1] ) {
      $oref->{rateVariation}->[1] = "10"; # bins for Beta-Gamma and GDD
    }
  } 

  my $gard_run = "";
  my $mpi = "";
  my $hostname = `hostname`;
  if( $hostname =~ /compute-0/ ){
    $mpi = "/usr/mpi/gcc/openmpi-1.2.8/bin/mpirun";
  } else {
    $mpi = "mpirun";
  }

  my $kernal = $basepath."HYPHYMPI";
  if ( $oref->{rateVariation}->[0] eq "1" ) {
    $gard_run = `(echo $fasta_path; echo $fasta_path; echo $oref->{model}->{char}; echo $oref->{rateVariation}->[0]; echo $results_path)|'$mpi' -np 16 '$kernal' $batchfile`;
#    $gard_run = "(echo $fasta_path; echo $fasta_path; echo $oref->{model}->{char}; echo $oref->{rateVariation}->[0]; echo $results_path)|'$mpi' -np 16 '$kernal' $batchfile";
  } else {
    $gard_run = `(echo $fasta_path; echo $fasta_path; echo $oref->{model}->{char}; echo $oref->{rateVariation}->[0]; echo $oref->{rateVariation}->[1]; echo $results_path)|'$mpi' '-np 16' '$kernal' $batchfile`;
#    $gard_run = "(echo $fasta_path; echo $fasta_path; echo $oref->{model}->{char}; echo $oref->{rateVariation}->[0]; echo $oref->{rateVariation}->[1]; echo $results_path)|'$mpi' '-np 16' '$kernal' $batchfile";
  }
# system $gard_run;

  if( $gard_run =~ /ERROR: Too few sites for c-AIC inference/ ){
    $rref->{output} = "Too few sites for c-AIC inference";
  } elsif( $gard_run =~ /Total run time/ ){
    &slchen::from_phylip_names($href, $name_hash);
    $gard_run = &slchen::from_phylip_names($gard_run, $name_hash);
    $rref->{output} = $gard_run;
  
    my $splits_file = $results_path."_splits";
    my $splits_path = Cwd::abs_path($splits_file);
    my $read_splits = do {
      open (my $sf, $splits_path) or die $!;
      local $/ = undef;
      <$sf>;
    };
  
    $read_splits = &slchen::from_phylip_names($read_splits, $name_hash);
    $rref->{splits} = $read_splits; 
  } elsif( $gard_run !~ /Segmentation fault/ ){
    $rref->{output} = "Too few sites for c-AIC inference";
  } else {
    $rref->{output} = "GARD did not run successfully!";
  }

  system "rm -rf $tempdir";
  return $rref;
}
# END RUN_GARD



# given a reference to a hash, get variables 
# run GARDProcessor and return the results
# Usage: run_gard_processor( $reference to seq hash which also has ouptut from gard ); scalar value
sub run_gard_processor {
  my $basepath = rashmi::hyphy_basepath;
  my $batchfile = $basepath."lib/hyphy/TemplateBatchFiles/GARDProcessor.bf";

  my $tempdir = File::Temp::tempdir(CLEANUP => 1);
  my $fasta = $tempdir."/fasta.temp";
  my $fasta_path = Cwd::abs_path($fasta);
  my $splits = $tempdir."/splits.temp";
  my $splits_path = Cwd::abs_path($splits);

  my $href = $_[0];
  my $rref = $_[1];

  my $name_hash = slchen::to_phylip_names($href);

  # define defaults
  if ( !defined $href ){
    die "No alignment provided! Aborting!\n";
  } else {
    open (FF, ">".$fasta);
    foreach my $key (sort keys %$href){
      print FF ">$key\n";
      print FF $href->{$key}, "\n";
    }
    close FF;
  }

  if( !defined $rref->{splits} ){
    die "No splits results provided! Aborting!\n";
  } else {
    foreach my $key ( keys %$name_hash ){
      $rref->{splits} =~ s/$name_hash->{$key}/$key/g;
    }
    open (SF, ">".$splits);
    print SF $rref->{splits};
    close SF;
  }

  my $kernal = $basepath."HYPHYMP";
  my $gardP_run = `(echo $fasta_path; echo $splits_path)|'$kernal' $batchfile`;
  #my $gardP_run = "(echo $fasta_path; echo $splits_path)|'$kernal' $batchfile";
  #system $gardP_run;

  &slchen::from_phylip_names($href, $name_hash);
  $gardP_run = &slchen::from_phylip_names($gardP_run, $name_hash);
  if( $gardP_run =~ /Mean splits identify/ ){
    $rref->{gp_output} = $gardP_run;
  } else {
    $rref->{gp_output} = "Gard processor did not run successfully!";
  }

  system "rm -rf $tempdir";
  return $rref;
}
# END RUN_GARD_PROCESSOR



# given a reference to a hash, get variables
# Parse GARD & GARD processor results
# Usage: parse_gard( $reference to hash ); scalar value
sub parse_gard {

  my $rref = $_[0];

  my @data = ();
  my @bp_pvalue = ();
  my @partition_length = ();
  my @all_pvalue = ();
  my @all_bp = ();
  my $significant_bp = "";

  # define defaults
  if ( !defined $rref->{gp_output} ){
    die "Results to be parsed not provided! Aborting!\n";
  } else {
    @data = split ("\n", $rref->{gp_output});
  }

  if( !defined $rref->{cutoff} ){
    $rref->{cutoff} = 0.05;  
  }

  foreach my $d (0..$#data){
    # Full alignment length
    if($data[$d] =~ /^\s*Sites\s*:(\d+)/){
      $rref->{align_length} = $1;
    }

    # parse breakpoints and p-values
    if( $data[$d] =~ /^Breakpoint\s+\|\s+LHS Raw p\s+\|\s+LHS adjusted p\s+\|\s+RHS Raw p\s+\|\s+RHS adjusted p/ ) {
      $d++;
      while( $data[$d] !~ /At p/ && $data[$d] !~ /^$/){
        $data[$d] =~ s/\s+//g;
        my @digits = split(/\|/, $data[$d]);
        push @bp_pvalue, $digits[0];
        push @bp_pvalue, $digits[2];
        push @bp_pvalue, $digits[4];
        $d++;
      }
    }

    # length of each partition (match Partition 1 : 588 sites)
    #if($data[$d] =~ /^\s*Partition\s+\d+\s+.\s+(\d+)\s+sites/g){
    #  push @partition_length, $1;
    #}
    # P-values
    #if($data[$d] =~ /^\s*KH Testing partition/ ) {
    #  if($d+2 <= $#data && $data[$d+1] =~ /^\s*Tree\s+\d+.*p-value = (\d+\.?\d*)/){
    #    # to match Tree 1 base LRT = 24.8763. p-value = 0.0014
    #    push @all_pvalue, $1;
    #    if($data[$d+2] =~ /^\s*Tree\s+\d+.*p-value = (\d+\.?\d*)/){
    #      push @all_pvalue, $1;
    #    }
    #  }
    #}
  }

  # check if max of two p-value is < cutoff ( default: 0.05 )
  for ( my $i=0; $i < scalar(@bp_pvalue); $i+=3){
    if( $bp_pvalue[$i+1] < $bp_pvalue[$i+2] ){
      if( $bp_pvalue[$i+2] <= $rref->{cutoff} ){
        $significant_bp .= $bp_pvalue[$i]."\t".$bp_pvalue[$i+1]."\t".$bp_pvalue[$i+2]."\n";
      }
    } elsif( $bp_pvalue[$i+1] <= $rref->{cutoff} ){
        $significant_bp .= $bp_pvalue[$i]."\t".$bp_pvalue[$i+1]."\t".$bp_pvalue[$i+2]."\n";
    }
  }

  # get breakpoints
  #my $break = 0;
  #for ( my $i=0; $i < scalar(@partition_length); $i++){
  #  $break = $break + $partition_length[$i];
  #  if ( $break != $rref->{align_length} ) {
  #    push @all_bp, $break;
  #  }
  #}
  # check if max of two p-values is < $cutoff (default: 0.01)
  #for ( my $j=0; $j < scalar(@all_bp); $j++ ) {
  #  if ( $all_pvalue[$j*2] < $all_pvalue[($j*2)+1] ) {
  #    if ( $all_pvalue[($j*2)+1] <= $rref->{cutoff} ) {
  #      $significant_bp .= $all_bp[$j]."\t".$all_pvalue[$j*2]."\t".$all_pvalue[($j*2)+1]."\n";
  #    }
  #  } elsif ( $all_pvalue[$j*2] <= $rref->{cutoff} ) {
  #      $significant_bp .= $all_bp[$j]."\t".$all_pvalue[$j*2]."\t".$all_pvalue[($j*2)+1]."\n";
  #  } 
  #}

  $rref->{significant_bp} = $significant_bp;

  return $rref;
}
# END PARSE_GARD



# given a reference to a hash, get sequence and breakpoints
# and fragment the alignment at the breakpoints
# Usage: fragment {reference to seq hash, reference to breakpoints hash}
sub fragment {
  my $href = $_[0];
  my $bpref = $_[1];
  my @break = (0);
  my @fragments = ();
  my @frags_temp = ();
  my $alignment = "";

  # get breakpoints from GARD results
  if( $bpref->{significant_bp} ne "" ){
    my @bp_db = split (/\s/, $bpref->{significant_bp});
    if( scalar(@bp_db) < 3 ){
      foreach (@bp_db){
        push @break, $_ if $_ !~ /\./;
      }
    } else{
      for (my $b=0; $b<scalar(@bp_db)-2; $b+=3){
        push @break, $bp_db[$b]; 
      }
    }
    # get length of alignment
    my @names = sort keys %$href;
    my $alignment_length = length($href->{$names[0]});
    push @break, $alignment_length;

    @break = slchen::sortu(@break);

    # fragment the alignment
    foreach my $b (1..$#break){
#print "$break[$b]\n";
      next if $break[$b] > $alignment_length;
      next if $break[$b] == 0;
      $alignment = "";
      my $current = 0;

      # here we try to break such that, if a breakpoint is exactly b/w two codons
      # then we make a clean cut. However, if a bp is in b/w a codon, we eliminate
      # the whole codon because we cant know which fragment to put this codon with 
      if( ($break[$b]-$break[$b-1])%3 !=0 ){
#print $break[$b]-$break[$b-1], "\n";
        if( (($break[$b]+1)-$break[$b-1])%3 == 0 ){
#print $break[$b]+1-$break[$b-1], "\n";
          #$break[$b] = $break[$b]+1;
          $current = $break[$b]-2;
          $break[$b] = $break[$b]+1; 
        } elsif( (($break[$b]-1)-$break[$b-1])%3 == 0 ){
          #$break[$b] = $break[$b]-1;
          $current = $break[$b]-1;
          $break[$b] = $break[$b]+2; 
        }
      } else {
        $current = $break[$b];
      }
      foreach my $n (0..$#names){
#print "$names[$n]\n";
#print "$break[$b]\n";
        my $id = $names[$n];
        $alignment .= ">$id\n";

        if( $current == $break[$b] ){
          $alignment .= substr($href->{$names[$n]}, $break[$b-1], $break[$b]-$break[$b-1]);
          print STDERR "error: length ", ($break[$b]-$break[$b-1]), " from ",$break[$b-1]+1,
                     " to $break[$b] is not divisible by 3\n", join (" ", @break), "\n"
                     if ($break[$b]-$break[$b-1]) % 3 != 0;
        } else {
          $alignment .= substr($href->{$names[$n]}, $break[$b-1], $current-$break[$b-1]);
          print STDERR "error: length ", ($break[$b]-$break[$b-1]), " from ",$break[$b-1]+1,
                     " to $break[$b] is not divisible by 3\n", join (" ", @break), "\n"
                     if ($current-$break[$b-1]) % 3 != 0;
        }
        # breaks from 0->x-1, x->y-1, y->(length-1), x and y are breakpoints
        $alignment .= "\n";
      }
      push @fragments, $alignment;
    }

    @frags_temp = @fragments;
    foreach my $i (0..$#frags_temp){
      my @seqs = split("\n", $frags_temp[$i]);
      if( $seqs[0] =~ />/ && ($seqs[1] =~ />/ || $seqs[1] eq "")){
        splice(@fragments, $i, 1);
      }
    }
  }
  return @fragments;
}
# END FRAGMENT



# given a reference to a hash, get variables
# run NielsenYang.bf - for site and branch analysis (Model 1a and 2 of PAML)
# Usage: run_NielsenYang( $seq_hash, $tree_hash, $options_hash ); scalar value
sub run_NielsenYang {
  my $basepath = rashmi::hyphy_basepath;
  my $batchfile = $basepath."lib/hyphy/TemplateBatchFiles/NielsenYang.bf";

  my $tempdir = File::Temp::tempdir(CLEANUP => 1);
  my $fasta = $tempdir."/fasta.temp";
  my $fasta_path = Cwd::abs_path($fasta);
  my $tree = $tempdir."/tree.temp";
  my $tree_path = Cwd::abs_path($tree);
  my $results = $tempdir."/results.temp";
  my $results_path = Cwd::abs_path($results);

  my $href = $_[0];
  my $tref = $_[1];
  my $oref = $_[2];
  my $rref = {};

  my $name_hash = slchen::to_phylip_names($href);

  # define defaults
  if ( !defined $href ){
    die "No alignment provided! Aborting!\n";
  } else {
    open (FF, ">".$fasta);
    foreach my $key (sort keys %$href){
      print FF ">$key\n";
      print FF $href->{$key}, "\n";
    }
    close FF;
  }

  if( !defined $tref->{tree} ){
    die "No tree provided! Aborting!\n";
  } else {
    foreach my $key (keys %$name_hash){
      $name_hash->{$key} =~ s/\s+//g; # this is because if the name has space, phyml will ignore it and put the words together but 
                                      # the name_hash still keeps the space - thus, will not substitute it with the phylip name`
      $tref->{tree} =~ s/$name_hash->{$key}/$key/g;
    }
    open (TF, ">".$tree);
    print TF $tref->{tree};
    close TF;
  }

  if( !defined $oref->{genetic_code} ) {
    $oref->{genetic_code} = "1"; } 
  if ( !defined $oref->{distributions}->[0] ) {
    $oref->{distributions} = [2, 2, 3]; } 
    # 2->to chose distr., 2->for neutral distr., 3->for selection distr.
  if ( !defined $oref->{initial_values} ) {
    $oref->{initial_values} = "1"; } 
  if( !defined $oref->{model} ) {
    $oref->{model} = "4"; } 
  if( ($oref->{model} eq "5") || ($oref->{model} eq "6") ) {
    if ( !defined $oref->{custom_model} ) {
      $oref->{custom_model} = "010010"; } # for HKY85
  }
  if( !defined $oref->{cutoff} ) {
    $oref->{cutoff} = "0.5"; } 
  if( !defined $oref->{done} ) {
    $oref->{done} = "d"; }

  # run Hyphy - NielsenYang
  my $run_ny = "";
  my $echo_string = "echo $oref->{genetic_code}; echo $fasta_path; echo $tree_path";
  foreach my $d (0..$#{$oref->{distributions}}){
    $echo_string .= "; echo $oref->{distributions}->[$d]";
  }
  $echo_string .= "; echo $oref->{done}; echo $oref->{initial_values}; echo $oref->{model}";
  if( ($oref->{model} ne "5") || ($oref->{model} ne "6") ) {
    $echo_string .= "; echo $oref->{cutoff}; echo $results_path";
  } else {
    $echo_string .= "; echo $oref->{custom_model}; echo $oref->{cutoff}; echo $results_path";
  }

  my $kernal = $basepath."HYPHYMP";
  $run_ny = `($echo_string)|'$kernal' $batchfile`;

  &slchen::from_phylip_names($href, $name_hash);
  $tref->{tree} = &slchen::from_phylip_names($tref->{tree}, $name_hash);
  $run_ny = &slchen::from_phylip_names($run_ny, $name_hash);

  $rref = parse_NielsenYang($run_ny);

  system "rm -rf $tempdir";

  return $rref;
}
# END RUN_NIELSENYANG



# given a results string parse output from NielsenYang.bf
sub parse_NielsenYang {
  my $data = $_[0];
  my $mref = {};
  my $results_path = "";
  my @lines;
  my @read_file;

  #---------------------------------------------------
  # parse terminal output
  @lines = split ("\n", $data);

  foreach my $l (0..$#lines) {
    # parse for results filepath
    if ( $lines[$l] =~ /^\s*Detailed results including sites with dN\/dS>1 will be written to$/ ) {
      chomp $lines[$l+1];
      $results_path = $lines[$l+1];
    }

    # parse for lnL, dN/dS and Nparams for each model
    if( $lines[$l] =~ /\s*Log likelihood\s*/ ) {
      foreach my $n ($l+1..$#lines) { # from next line till end
        if ( $lines[$n] =~ /\|/ ) {
          my @variables = split(/\|/, $lines[$n]); # the first element is empty -> before the first pipe
          # get rid of all extra spaces
          foreach my $v (0..$#variables) {
            $variables[$v] =~ s/\s//g;
          }
          # get model no and desc
          my @model_no = split(/\./, $variables[1]); 
          # get log likelihood
          $mref->{$model_no[1]}->{LnL} = $variables[2];
          # get dN/dS
          $mref->{$model_no[1]}->{dNdS} = $variables[3];
          # get no of parameters
          $mref->{$model_no[1]}->{Nparams} = $variables[4];
        }
      }
    }
  }
  #--------------------------------------------------------
  # parse main results file
  if ( !defined $results_path ) {
    print "Bad results! Cannot find results file path! Aborting!\n";
  } else {
    open RF, $results_path or die $!;
    @read_file = <RF>;
    close RF;
  }

  my @model_results = split (/\n\*\*\*/, join("", @read_file)); # split results for each model
  my @rates;
  my @weights;
  foreach my $m (1..$#model_results) {
    foreach my $model ( keys %$mref ) {
      if ( $model_results[$m] =~ /$model/ ) {
        $mref->{$model}->{output} = $model_results[$m];

        @lines = split ("\n", $model_results[$m]);
        foreach my $l (0..$#lines) {
          # get Ntime
          if( $lines[$l] =~ /\>Done in (\d+)/ ){
            $mref->{$model}->{Ntime} = $1;
          } 
          if( $lines[$l] =~ /^Rate\[\d+\]=\s*(\d+\.\d+)\s*\(weight=(\d+\.\d+)\)$/ ){
            push @rates, $1;
            push @weights, $2;  
          } 
        }
        foreach my $r (0..$#rates) {
          my $omega = "w_$r";
          my $p = "P_$r";
          $mref->{$model}->{$omega} = shift (@rates);
          $mref->{$model}->{$p} = shift (@weights);
        }
      }
    }
  }
  #------------------------------------------------------------
  # parse nex file for each model
  my $count = 0;
  foreach my $model ( sort keys %$mref ){
    if ( !defined $results_path ) {
      print "Bad results! Cannot find results file path! Aborting!\n";
    } else {
      open MF, $results_path."_MODEL_".$count.".nex" or die $!;
      @read_file = <MF>;
      close MF;
    }

    my @hyphy_results = split (/BEGIN HYPHY;/, join("", @read_file)); # split results for each model 
    $mref->{$model}->{output} = $mref->{$model}->{output}.$hyphy_results[1];
    @lines = split ("\n", $hyphy_results[1]);
    foreach my $l (0..$#lines) {
      if( $lines[$l] =~ /\s*global\s*kappa=(-?\d+\.?\d*.*);/ ){
        $mref->{$model}->{kappa} = $1;
        last;
      }
    }
    $count++;
  }

  return $mref;
}
# END PARSE_NIELSENYANG



# given a reference to a hash, get variables
# run YangNielsenBranchSite2005.bf - for branch site analysis 
# Usage: run_YangNielsen_bs( $reference to hash ); scalar value
sub run_YangNielsen_bs {
  my $basepath = rashmi::hyphy_basepath;
  my $batchfile = $basepath."lib/hyphy/TemplateBatchFiles/YangNielsenBranchSite2005.bf";

  my $tempdir = File::Temp::tempdir(CLEANUP => 1);
  my $fasta = $tempdir."/fasta.temp";
  my $fasta_path = Cwd::abs_path($fasta);
  my $tree = $tempdir."/tree.temp";
  my $tree_path = Cwd::abs_path($tree);
  my $results = $tempdir."/results.temp";
  my $results_path = Cwd::abs_path($results);

  my $href = $_[0];
  my $tref = $_[1];
  my $oref = $_[2];
  my $rref = {};

  my $name_hash = slchen::to_phylip_names($href);

  # define defaults
  if ( !defined $href ){
    die "No alignment provided! Aborting!\n";
  } else {
    open (FF, ">".$fasta);
    foreach my $key (sort keys %$href){
      print FF ">$key\n";
      print FF $href->{$key}, "\n";
    }
    close FF;
  }

  if( !defined $tref->{tree} ){
    die "No tree provided! Aborting!\n";
  } else {
    foreach my $key ( keys %$name_hash ){
      $tref->{tree} =~ s/$name_hash->{$key}/$key/g;
    }
    open (TF, ">".$tree);
    print TF $tref->{tree};
    close TF;
  }

  if( !defined $oref->{genetic_code} ) {
    $oref->{genetic_code} = "1"; }
  if( !defined $oref->{model} ) {
    $oref->{model} = "1"; }
  if( $oref->{model} ne "2" ) {
    if ( !defined $oref->{branch}->[0] ) {
      die "No foreground branch chosen! Aborting!";
    }
  } 
  if( !defined $oref->{done} ) {
    $oref->{done} = "d"; }

  # run Hyphy - YangNielsenBranchSite
  my $run_bs = "";
  my $kernal = $basepath."HYPHYMP";
  if( $oref->{model} eq "2" ) {
    $run_bs = `(echo $oref->{genetic_code}; echo $fasta_path; echo $tree_path; echo $oref->{model}; echo $results_path)|'$kernal' $batchfile`;
  } else {
    my $echo_string = "echo $oref->{genetic_code}; echo $fasta_path; echo $tree_path; echo $oref->{model}";
    foreach my $b (0..$#{$oref->{branch}}){
      $echo_string .= "; echo $oref->{branch}->[$b]"; 
    }
    $echo_string .= "; echo $oref->{done}; echo $results_path";

    $run_bs = `($echo_string)|'$kernal' $batchfile`;
  }

  &slchen::from_phylip_names($href, $name_hash);
  $tref->{tree} = &slchen::from_phylip_names($tref->{tree}, $name_hash);
  $run_bs = &slchen::from_phylip_names($run_bs, $name_hash);

  $rref = parse_YangNielsen_bs($run_bs, $results_path);

  return $rref;

  system "rm -rf $tempdir";
}
# END RUN_YANGNIELSEN_BS



# given a results string parse output from YangNielsenBranchSite2005.bf
sub parse_YangNielsen_bs {
  my $data = $_[0];
  my $results_path = $_[1];
  my $pref = {};
  my @lines;
  my @read_file;

  #---------------------------------------------------
  # parse terminal output
  my (@wr, @rates, @weights);
  my $concat = "";

  @lines = split ("\n", $data);

  foreach my $l (0..$#lines) {
    # get base compositions
    if( $lines[$l] =~ /Base composition:/ ){
      $concat .= $lines[$l].$lines[$l+1].$lines[$l+2];
      $concat .= $lines[$l+3].$lines[$l+4];
    }

    # parse for kappa 
    if( $lines[$l] =~ /Obtaining nucleotide branch lengths and kappa to be used as starting values/ ){
      $concat .= "Obtaining nucleotide branch lengths and kappa to be used as starting values";
      $concat .= $lines[$l+2].$lines[$l+3];
    } 
    if( $lines[$l] =~ /^\s*kappa\s*=\s*(\d+.\d+)/ ) {
      $pref->{kappa} = $1;
    }

    if( $lines[$l] =~ /\d+ foreground branch\(es\) set to:/ ){
      foreach my $n ($l..$#lines-3){
        $concat .= $lines[$n];
      }
    }

    # parse for lnl
    if( $lines[$l] =~ /^\s*Log Likelihood\s*=\s*(\-?\d+.\d+)/ ) {
      $pref->{LnL} = $1;
    }
    # parse for omega and P
    if( $lines[$l] =~ /\s*Inferred rate distribution:/ ){
      foreach my $n ($l+1..$#lines-3) { # from next line till end
        @wr = split(/\s+/, $lines[$n]);
        if( $wr[2] =~ /\d.$/){ # for class 0 and 1
          push @rates, $wr[5];
          push @weights, $wr[8];
        }
        if( $wr[2] =~ /2a.$/){ # for class 2a
          push @rates, $wr[10];
          push @weights, $wr[13];
        }
      }
    }
  }
  foreach my $r (0..$#rates) {
    my $omega = "w_$r";
    my $p = "P_$r";
    $pref->{$omega} = shift (@rates);
    $pref->{$p} = shift (@weights);
  }
  
  $pref->{output} = $concat;
  #--------------------------------------------------------
  # parse main results file
  if ( !defined $results_path ) {
    print "Bad results! Cannot find results file path! Aborting!\n";
  } else {
    open RF, $results_path or die $!;
    @read_file = <RF>;
    close RF;
  }

  my $sites = join("", @read_file);
  $pref->{output} = $pref->{output}.$sites;

  return $pref;
}
# END PARSE_YANGNIELSEN_BS



# given a reference to a sequence hash and options hash, get variables to run phyml
sub run_phyml {
  my $href = $_[0];
  my $rref = {};
  my $name_hash = {};
  my $tree_array = "";

  my $tempdir = File::Temp::tempdir(CLEANUP => 1);
  my $fasta = $tempdir."/fasta.temp";
  my $phylip_file = $tempdir."/phylip.temp";

  # define defaults
  if ( !defined $href ){
    die "No alignment provided! Aborting!\n";
  } else {
    # convert names of sequences to phylip compatible format
    $name_hash = slchen::to_phylip_names($href);

    open (FF, ">".$fasta);
    foreach my $key (sort keys %$href){
      print FF ">$key\n";
      print FF $href->{$key}, "\n";
    }
    close FF;
  }

  # convert from fasta to phylip format
  my $phylip_format = `fasta2phy.pl $fasta`;

  open (PF, ">".$phylip_file);
  print PF $phylip_format;
  close PF; 

  my $hostname =`hostname`;
  chomp $hostname;
  my $phyml_output = "";
  if( $hostname =~ /compute-0/){
    $phyml_output = `/home/gisv88/Phyml/phyml-20100720/bin/phyml -i $phylip_file -p`;
  } else {
    $phyml_output = `phyml -i $phylip_file -p`;
    #$phyml_output = "phyml -i $phylip_file -p";
    #system $phyml_output;
  }
  #----------------------------------------------------------------
  # parse results file for tree
  my $tree_file = $phylip_file."_phyml_tree.txt";
  my $read_tree = do {
    open (my $tf, $tree_file) or die $!;
    local $/ = undef;
    <$tf>;
  };
  chomp $read_tree;

  # convert the original seq hash and the tree back to original names
  &slchen::from_phylip_names($href, $name_hash);
  my $tree_fasta = slchen::from_phylip_names($read_tree, $name_hash);

  $tree_array = [$tree_fasta];
  my @tree_parsed = rm_bootstrap($tree_array);
 
  $rref->{tree} = $tree_fasta;
  $rref->{tree_parsed} = $tree_parsed[0];

  #----------------------------------------------------------------
  # parse stats file for variables like parsimony, tree size etc
  my $stats_file = $phylip_file."_phyml_stats.txt";
  my @read_stats;
  open (SF, $stats_file) or die $!;
  @read_stats = <SF>;
  close SF;

  foreach my $rs (0..$#read_stats){
    if( $read_stats[$rs] =~ /Parsimony:\s+\t+(\d+)/ ){
      $rref->{parsimony} = $1; 
    }
    if( $read_stats[$rs] =~ /Tree size:\s+\t+(\d+.\d+)/){
      $rref->{tree_size} = $1; 
    }
  }
  $rref->{type} = "ML";
  return $rref;

  system "rm -rf $tempdir";
}
# END RUN_PHYML


# remove bootstrap values from phyml tree
# hyphy segfaults if bootstrap values are present in tree
# makes use of procedures from Slctree.pm
sub rm_bootstrap {
  my ($tree_ref) = @_;
  my $treedata = {};
  my $probability = {};
  my $optiondata = {};

  &Slctree::parse_tree($tree_ref, $treedata, $probability, $optiondata);
  my @tree = full_tree_text($treedata);
    
  return @tree;
}
# END RM_BOOTSTRAP


# adapted from Slctree.pm to get tree string without bootstrap values
sub full_tree_text {
  # this calls tree_text_wo_bootstrap
  # expects a top level tree data hash where first level keys are tree number
  # adds on the outside parentheses and probabilities and semicolon
  my ($dref, $pref) = @_;
  my $key;
  my $tree;
  my @tree;
  foreach $key (keys %$dref) {
    $tree = "(" . tree_text_wo_bootstrap($dref->{$key}) . ")";
    if (defined $pref->{$key} && $pref->{$key} > 0) {
      $tree .= "[" . $pref->{$key} . "]";
    }
    $tree .= ";";
    push @tree, $tree;
  }
  return @tree;
}


# adapted from Slctree.pm to get tree string without bootstrap values
sub tree_text_wo_bootstrap {
  # print out tree text from tree data hash
  my $nodeprefix = "Nodenumber";
  my ($dref) = @_;
  my $return;
  my $key;
  my $node;
  my @node;
  # each level should have at least special keys - branchlength, sortvalue, bootstrap
  foreach $key (keys %$dref) {
    next if $Slctree::special_key{$key};
    if ($key =~ /$nodeprefix/) {
      $node = "(" . tree_text_wo_bootstrap($dref->{$key}) . ")";
      #if (defined $dref->{$key}->{bootstrap} && $dref->{$key}->{bootstrap} > 0) {
        #$node .= $dref->{$key}->{bootstrap};
      #}
      if (defined $dref->{$key}->{branchlength} && $dref->{$key}->{branchlength} >= 0) {
        $node .= ":" . $dref->{$key}->{branchlength};
      }
    } else {
      # we are at a terminal leaf
      $node = $key;
      # only look for branchlength, since we can't have a bootstrap for a terminal node
      if (defined $dref->{$key}->{branchlength} && $dref->{$key}->{branchlength} >= 0) {
        $node .= ":" . $dref->{$key}->{branchlength};
      }
    }
    push @node, $node;
  }

  return join (",", @node);
}
# END TREE_TEXT_WO_BOOTSTRAP



# given a tree, find the number(s) assigned to its leaves and nodes by hyphy
# return the number list
# Usage: tree2num { reference to hash with tree }
sub tree2num {
  my $basepath = rashmi::hyphy_basepath;
  my $batchfile = "hyphy-ps/parse_tree.bf";
  my $bf_path = Cwd::abs_path($batchfile);

  my $tref = $_[0];

  my $tempdir = File::Temp::tempdir(CLEANUP => 1);
  my $tree = $tempdir."/tree.temp";
  my $tree_path = Cwd::abs_path($tree);

  # define defaults
  if( !defined $tref->{tree} ){
    die "No tree provided! Aborting!\n";
  } else {
    open (TF, ">".$tree);
    print TF $tref->{tree};
    close TF;
  }

  my $output = `(echo $tree_path)|'HYPHYMP' $bf_path`; 
  
  my @parse = split (/::\n/, $output);       # remove "please select tree file" line in output from hyphy
  my @num_list = split (/Check/, $parse[1]); # remove "Check message.log" line in output from hyphy 
  chomp $num_list[0];                        # remove extra newline created by hyphy 

  return ($num_list[0]);
} 
# END TREE2NUM


# given a fasta, output interleaved nex format compatible with MR.bayes 
# Usage: fasta2nex_bayes { fasta hash }
sub fasta2nex_bayes {

  my $seq = shift @_;
  my $nex = "";
  $nex .= "#NEXUS\n";
  $nex .= "BEGIN DATA;\n";

  my @keys = sort keys %$seq;
  my $taxnum = scalar(@keys);
  my $chars = length($seq->{$keys[0]});
  $nex .= "\tDIMENSIONS NTAX=$taxnum NCHAR=$chars;\n";

  $nex .= "\tFORMAT DATATYPE=DNA INTERLEAVE=yes GAP=-;\n";
  $nex .= "MATRIX\n";
  for (my $i =0; $i <= $chars; $i=$i+80){
    foreach my $key (@keys){
      $nex .= "$key ";
      $nex .= substr($seq->{$key}, $i, 80);
      $nex .= "\n";
    }
    $nex .= "\n";
  }
  $nex .= ";\n";
  $nex .= "END;\n";

  return ($nex);
}
# END FASTA2NEX_BAYES


# given an alignment, remove duplicates and pass back unique set
sub rm_duplicates {
  my $seq = shift @_;
  my $href = {};
  $href = slchen::fasta2hash($seq);

  my @uniq = ();
  foreach my $key (sort keys %$href ){
    if( $href->{$key} =~ /(A|T|G|C)/ ){
      if(scalar(@uniq) == 0){
        push @uniq, ">$key";
        push @uniq, $href->{$key};
      } elsif ( grep {$_ eq  $href->{$key}} @uniq ){
      } else {
        push @uniq, ">$key";
        push @uniq, $href->{$key};
      }
    }
  } 
  $seq = join ("\n", @uniq);
  return $seq;
}
# END RM_DUPLICATES

# given an alignment, remove duplicates and pass back unique set
# here duplicates include the ones with gaps as well
# ATG--CA is a duplicate of ATGGTCA
sub rm_gap_duplicates {
  my $seq = shift @_;
  my $href = {};
  my @uniq = ();
  # remove identical duplicates
  $seq = rm_duplicates($seq);
  $href = slchen::fasta2hash($seq);

  # sequences with no gaps at all, pass through
  # record and delete them from the hash
  foreach my $key (sort keys %$href){
    $href->{$key} =~ s/N/-/g;
    if( $href->{$key} !~ /-/ ){ 
      push @uniq, ">$key";
      push @uniq, $href->{$key};
      delete $href->{$key};
    }
  }

  # check the remaining gapped sequences for a duplicate
  SEQ: foreach my $key (sort keys %$href){
    my $char = "-";
    my $offset = 0;
    my @index = ();

    $href->{$key} =~ s/N/-/g;
    # get positions of gaps
    my $gin = index($href->{$key}, $char, $offset);
    while( $gin != -1 ){
      if( scalar(@index) == 0 ){
        push @index, $gin;
      } elsif( $gin != $index[scalar(@index)-1]+1 ){
        push @index, $gin+1;
      }
      $offset = $gin+1;
      $gin = index($href->{$key}, $char, $offset);
    }
    push @index, length($href->{$key});

    $offset = 0;
    INDEX: foreach my $in (@index){
      my $substr = substr($href->{$key}, $offset, $in-$offset);
      $substr =~ s/N//g;
      $substr =~ s/-//g;

      if($substr ne ""){
        my @ref = ();
        foreach my $uniq (@uniq){
          next if $uniq =~ />/;

          my $ref_sub = substr($uniq, $offset, $in-$offset);
          $ref_sub =~ s/N/-/g;
          push @ref, $ref_sub;
        }

        if( grep{ $_ =~ $substr } @ref ){
        } else {
          push @uniq, ">$key";
          push @uniq, $href->{$key};
          next SEQ;
        }
      }
      $offset = $in;
    }
  }
  
  $seq = join("\n", @uniq);
  return $seq;
}
# END RM_GAP_DUPLICATES



#-------------------------------------------------------------
# Perl uses 'eval' to process the code. Thus for 'eval' to
# evaluate to TRUE and not fail, we provide a 1
1;
