#!/usr/bin/perl
use warnings;
use strict;


foreach my $infile (@ARGV){
    if ($infile =~ /gz$/){
	open (IN,"zcat $infile |") or die $!;
    }
    else{
	open (IN,$infile) or die $!;
    }
    warn "Reading in file  >>> $infile <<<\n";

    my $outfile = $infile;
    $outfile =~ s/$/_multiple_counts.txt/;
    open (OUT,'>',$outfile) or die $!;
    warn "Writing output to  >>> $outfile <<<\n\n";

  my %counts_TTAGGG; # this corresponds to CCCTAAA on the forward strand
  my %counts_CCCTAA; # this is the TTTAGGG sequence on the forward strand
  my %counts_TTTTAA; # this is the reverse complement of TTTTAA

  my $total_base_count = 0;
  my $total_hexamer_count = 0;
  my $count=0;
  while (1){
    my $line1 = <IN>;
    my $seq = <IN>;
    my $line3 = <IN>;
    my $line4 = <IN>;
    last unless ($line1 and $seq and $line3 and $line4);
    chomp $seq;
    # my @bases = split //,$seq;

    my $TTAGGG_count = 0;
    my $CCCTAA_count = 0;
    my $TTTTAA_count = 0;

    ## counting the actual telomere repeat
    while ($seq =~ /TTAGGG/g){
      ++$TTAGGG_count;
    }
    ## counting the reverse complement of the telomere repeat (which does contain Cs!)
    while ($seq =~ /CCCTAA/g){
      ++$CCCTAA_count;
    }
    while ($seq =~ /TTTTAA/g){
      ++$TTTTAA_count;
    }
    $counts_TTAGGG{$TTAGGG_count}++;
    $counts_CCCTAA{$CCCTAA_count}++;
    $counts_TTTTAA{$TTTTAA_count}++;

    ++$count;
    if ($count%5000000 ==0){
      warn "Processed $count lines so far\n";
    }
  }
    #  last if ($count == 49999);
    
    print "Total number of sequences analysed: $count\n";
    print OUT "Total number of sequences analysed: $count\n";
    
    print  "Summary for TTAGGG\n";
    print OUT "Summary for TTAGGG\n";
    print OUT "\t$infile\n";
    foreach my $occurrence(sort {$a<=>$b} keys %counts_TTAGGG){
	print "$occurrence\t$counts_TTAGGG{$occurrence}\n";
	print OUT "$occurrence\t$counts_TTAGGG{$occurrence}\n";
    }
    
    print  "Summary for CCCTAA\n";
    print OUT "Summary for CCCTAA\n";
    print OUT "\t$infile\n";
    foreach my $occurrence(sort {$a<=>$b} keys %counts_CCCTAA){
	print "$occurrence\t$counts_CCCTAA{$occurrence}\n";
	print OUT "$occurrence\t$counts_CCCTAA{$occurrence}\n";
    }
    
    print  "Summary for TTTTAA\n";
    print OUT "Summary for TTTTAA\n";
    print OUT "\t$infile\n";
    foreach my $occurrence(sort {$a<=>$b} keys %counts_TTTTAA){
	print "$occurrence\t$counts_TTTTAA{$occurrence}\n";
	print OUT "$occurrence\t$counts_TTTTAA{$occurrence}\n";
    }
    
}
