#!/usr/bin/perl
use warnings;
use strict;

my @files = @ARGV;

foreach my $file (@files) {
		process_file($file);
}

sub process_file {
		my ($file) = @_;

		warn "Processing $file\n";

		my $chr = "";
		my $seq;

		my $outfile = $file;
		$outfile =~ s/\.txt//;

		$outfile .= "_composition.txt";

		open (IN,$file) or die $!;
		open (OUT,'>',$outfile) or die $!;

		my $header = <IN>;
		chomp $header;
		$header .= "\tG\tA\tT\tC\tN\n";
		print OUT $header;

		while (<IN>) {
				chomp;
				next unless ($_);
				my (undef,$this_chr,$start,$end) = split(/\t/);

				if ($chr ne $this_chr) {
						warn "Processing $file Chr $this_chr\n";
						$seq = read_seq($this_chr);
						$chr = $this_chr;
				}
				
				my $subseq = substr($seq,$start-1,($end-$start)+1);

				my $c = $subseq =~ tr/Cc//;
				my $g = $subseq =~ tr/Gg//;
				my $a = $subseq =~ tr/Aa//;
				my $t = $subseq =~ tr/Tt//;
				my $n = $subseq =~ tr/Nn//;

				print OUT $_."\t$g\t$a\t$t\t$c\t$n\n";

		}

		close IN;
		close OUT or die $!;

}

sub read_seq {
		my ($chr) = @_;

		my $file = "/bi/scratch/Genomes/Mouse/NCBIM37/Mus_musculus.NCBIM37.52.dna.chromosome.${chr}.fa";

		open (SEQ,$file) or die $!;
		my $header =<SEQ>;

		my $seq;
		while (my $line = <SEQ>) {
				chomp $line;
				$seq .= $line;
		}

		close SEQ;

		return $seq;
}
