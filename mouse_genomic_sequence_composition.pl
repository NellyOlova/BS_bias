#!/usr/bin/perl
use warnings;
use strict;

my $word_size = 1;
my $seq = '';
my %freqs;
my $total_number_of_words_counted = 0;

chdir ("/data/public/Genomes/Mouse/NCBIM37/") or die "Can't move to folder /data/public/Genomes/Mouse/NCBIM37: $!";
while (my $filename = <*.fa>){
  my $chrom = chromosome_number($filename);
  $seq = read_sequence_file($filename);
  print "The length of the concatenated chromosome $chrom is ",length($seq)," bp\n";
  process_sequence($seq);
}

my $outfile = 'mouse_genome_sequence_composition_single_bases.txt';
chdir ("/data/Projects/Miguel/") or die "Can't move to folder /data/Projects/Miguel/: $!";
open (OUT,'>',$outfile) or die $!;

print "\nNumber of different words found in the mouse genome for the word size ## $word_size ##: ",scalar(keys %freqs),"\n";
print OUT "Number of different words found in the mouse genome for the word size ## $word_size ##: ",scalar(keys %freqs),"\n";
print "Total number of words counted (discarded N containing words): $total_number_of_words_counted\n";
print OUT "Total number of words counted (discarded N containing words): $total_number_of_words_counted\n";
print "\nword\tcount\tpercent total\n";
print OUT "\nword\tcount\tpercent total\n";

calculate_averages();

sub calculate_averages {
  foreach my $word (sort { $freqs{$b} <=> $freqs{$a} } keys %freqs) {
    my $percentage = sprintf ("%.2f",100*$freqs{$word}/$total_number_of_words_counted);
    print "$word\t$freqs{$word}\t$percentage\n";
    print OUT "$word\t$freqs{$word}\t$percentage\n";
  }
}

sub process_sequence {
  my $sequence = shift;
  my $index = 0;
  my $counted = 0;
  my $stop = length($sequence);
  $stop -= $word_size;
  while ($index < $stop) {
    if ($index%5000000==0){
      warn "Current index number is $index\n";
    }
    my $word = substr($sequence,$index,$word_size);
    if (index($word,'N') < 0) {
      $freqs{$word}++;
      ++$counted;
    }
    ++$index;
  }
  print "Finished counting different words (checked $index different possibilities)\n";
  print "Number of sequences counted and entered into the hash: $counted\n";
  $total_number_of_words_counted += $counted;
}

sub read_sequence_file{
  my $filename = shift;
  my $sequence;
  warn "Now reading sequence data from $filename\n";
  open (MOUSE_CHROMOSOME,$filename)or die "Can't open $filename: $!";
  $_ = <MOUSE_CHROMOSOME>; #fasta header
  while (<MOUSE_CHROMOSOME>){
    chomp;
    $sequence .= uc$_;
  }
  close MOUSE_CHROMOSOME or die $!;
  return $sequence;
}
sub chromosome_number{
  my $filename = shift @_;
  if ($filename =~ /\.([^\.]+)\.fa$/){
    return $1;
  }
  else{
    die "Unable to extract the chromosome number: $filename!";
  }
}
