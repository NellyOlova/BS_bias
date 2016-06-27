#!/usr/bin/perl
use warnings;
use strict;

my @files = @ARGV;

my $pseudo_mappings = read_aliases();


foreach my $file (@files) {
		process_file($file);
}

sub read_aliases {
		open (IN,'/bi/scratch/Rebecca/Genomes/Harpegnathus_saltator/H3.3/aliases.txt') or die $!;

		my %pseudo;
		while (<IN>) {
				chomp;
				my ($scaffold,$pseudo,$offset) = split(/\t/);
				push @{$pseudo{$pseudo}},[lc($scaffold),$offset];
		}

		return \%pseudo;
}

sub get_scaffold {
		my ($pseudo,$position) = @_;

		unless (exists $pseudo_mappings->{$pseudo}) {
				die "No mapping for pseudo chr '$pseudo'";
		}

		foreach my $index (0..$#{$pseudo_mappings->{$pseudo}}) {
				if ($pseudo_mappings->{$pseudo}->[$index]->[1] > $position) {


#						warn "Scaffold for $pseudo was ".
#								$pseudo_mappings->{$pseudo}->[$index-1]->[0] .
#								" with offset ".
#								$pseudo_mappings->{$pseudo}->[$index-1]->[1].
#								" next was ".
#								$pseudo_mappings->{$pseudo}->[$index]->[0] .
#								" with offset ".
#								$pseudo_mappings->{$pseudo}->[$index]->[1].
#								" diff was ".
#								($pseudo_mappings->{$pseudo}->[$index]->[1]-$pseudo_mappings->{$pseudo}->[$index-1]->[1]) .
#								"\n";

						return ($pseudo_mappings->{$pseudo}->[$index-1]->[0],$pseudo_mappings->{$pseudo}->[$index-1]->[1]);;
				}
		}

		# If we get here then we want to return the last segment
		return (${$pseudo_mappings->{$pseudo}}[-1]->[0],${$pseudo_mappings->{$pseudo}}[-1]->[1]);
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
#				next unless ($_);
				my (undef,$pseudo,$start,$end) = split(/\t/);

#				warn "Position was $pseudo $start $end\n";
				my ($this_chr,$offset) = get_scaffold($pseudo,$start);

				$start -= $offset;
				$end -= $offset;

#				warn "Corrected position was $this_chr $start $end\n";

				if ($chr ne $this_chr) {
#						warn "Processing $file Chr $this_chr\n";
						$seq = read_seq($this_chr);
						$chr = $this_chr;
				}

				if ($end >= length($seq)) {
						warn "Fell off end of scaffold by ".($end-length($seq))."\n";
						next;
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

		my $file = "/bi/home/andrewss/scratch/Genomes/Harpegnathos_saltator/HS3.3/Split_Scaffolds/${chr}.fasta";

		open (SEQ,$file) or die "Can't find '$file':$!";
		my $header =<SEQ>;

		my $seq;
		while (my $line = <SEQ>) {
				chomp $line;
				$seq .= $line;
		}

		close SEQ;

		return $seq;
}
