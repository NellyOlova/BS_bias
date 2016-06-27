#!/usr/bin/perl
use warnings;
use strict;

my $infile = shift @ARGV;
unless ($infile =~ /bam$/){
    die "Please provide a BAM file to continue!\n";
}
open (IN,"samtools view -h $infile |") or die $!;

my $outfile = $infile;
$outfile =~ s/\.bam$/_nonCpG_filtered.bam/;
open (OUT,"| samtools view -bS - > $outfile") or die $!;

my $removed = $infile;
$removed =~ s/$/_removed_seqs.txt/;
open (REMOVED,'>',$removed) or die $!;

my $count = 0;
my $kicked = 0;

while (<IN>){
    ++$count;
    if ($_ =~ /^\@/){
	print OUT;
	next;
    }
#    last if ($count == 50);
    my ($meth_call) = $1 if ($_ =~ /XM:Z:(.+?)\s/);
    # warn "$meth_call\n"; sleep(1);
    my $sequence_fails = 0;
    
    my $nonCpG_count = 0;
    
    my @chars = split //,$meth_call;
    
    foreach my $char(@chars){
	if ( ($char eq 'H') or ($char eq 'X') ){
	    ++$nonCpG_count;
	}
	elsif ( ($char eq 'z') or ($char eq 'Z') ){
	    # $consecutive_nonCpG_count = 0; # resetting the counter if there is any kind of CpG on the way
	}
	
	if ($nonCpG_count >= 3){
	    $sequence_fails = 1;
	    last;
	}
    }

    if ($sequence_fails){
	++$kicked;
	print REMOVED;
    }
    else{
	print OUT;
    }

}
warn "Analysed $count sequences in total\n";
my $percent = sprintf ("%.1f",$kicked/$count*100);
warn "Sequences removed because of too much non-CpG methylation: $kicked ($percent%)\n\n";

