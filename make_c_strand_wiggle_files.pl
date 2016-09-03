#!/usr/bin/perl
use warnings;
use strict;

my @files = </bi/scratch/Genomes/Mouse/NCBIM37/*fa>;


open(C,'>','c_content_NCBIM37.wig') or die $!;
open(G,'>','g_content_NCBIM37.wig') or die $!;
open(N,'>','n_content_NCBIM37.wig') or die $!;

print C 'track type=wiggle_0 name="C content" description="Top strand C content"';
print G 'track type=wiggle_0 name="G content" description="Top strand G content"';
print N 'track type=wiggle_0 name="N content" description="Top strand N content"';

print C "\n";
print G "\n";
print N "\n";

foreach my $file (@files) {
		process_file($file);
}

close C or die $!;
close G or die $!;
close N or die $!;

sub process_file {
		my ($file) = @_;

		my $chr;
		if ($file =~ /chromosome\.(\S+)\.fa/) {
				$chr = $1;
		}
		else {
				die "Couldn't get chr name from $file";
		}

		# UCSC can't tranlate MT names
		$chr = 'M' if ($chr eq 'MT');

		warn "Processing chr $chr\n";

		print C "fixedStep chrom=chr${chr} start=1 step=200 span=200\n";
		print G "fixedStep chrom=chr${chr} start=1 step=200 span=200\n";
		print N "fixedStep chrom=chr${chr} start=1 step=200 span=200\n";

		my $seq = read_seq($file);

		my @bins_c;
		my @bins_g;
		my @bins_n;

		for (0..((int(length($seq)/200))-1)) {

				my $subseq = substr($seq,200*$_,200);

				$bins_c[$_] = $subseq =~ tr/Cc//;
				$bins_g[$_] = $subseq =~ tr/Gg//;
				$bins_n[$_] = $subseq =~ tr/Nn//;;
		}
		
		print C join("\n",@bins_c),"\n";
		print G join("\n",@bins_g),"\n";
		print N join("\n",@bins_n),"\n";


}

sub read_seq {
		my ($file) = @_;
		open (IN,$file) or die $!;
		$_=<IN>;

		my $seq;
		while (<IN>) {
				chomp;
				$seq .= $_;
		}

		close IN;

		return $seq;
}
