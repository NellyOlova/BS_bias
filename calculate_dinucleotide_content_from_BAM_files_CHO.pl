#!/usr/bin/perl
use warnings;
use strict;
use Cwd;
$|++;

my $parent_dir = getcwd();
my %chromosomes;
my $genome_folder = "/bi/scratch/Genomes/Chinese_Hamster_PreEnsembl/";
# my $genome_folder = "/bi/scratch/Genomes/Chinese_Hamster_PreEnsembl/";
read_genome_into_memory ($genome_folder);

chdir $parent_dir or die;

my %positions;
my %nuc;


foreach my $in (@ARGV){
    %nuc = (); # resetting the dinucleotide counter
    open (IN,"samtools view -h $in |") or die $!;
    warn "Reading $in\n"; sleep(1);
    
    my $out = $in;
    die "Input file has to end in .bam!\n" unless ($out =~ s/bam$/dinucleotide_content.txt/);
    open (OUT,'>',$out) or die "Failed to write to output file $out\n";
    my $count = 0;
    
    my ($single,$paired) = test_file($in);
    # warn "returned $single and $paired\n";
    
    if ($single){
	warn "Determined the file to be single-end!\n";
    }
    elsif ($paired){
	warn "Determined the file to be paired-end!\n";
    }
    else{
	die "Failed to figure out SE or PE...\n";
    }
    
    while (<IN>){
	chomp;
	if (/^\@/) {
	    #  warn "$_\n";
	    next;
	}
	# last if ($count == 30);
	++$count;
	if ($count%500000 == 0){
	    warn "Processed $count lines\n";
	}
	my ($flag,$chr,$start,$cigar,$sequence) = (split(/\t/))[1,2,3,5,9];
	# warn "$flag\t$chr\t$start\t$cigar\t$sequence\n"; 
	
	### checking CIGAR string for insertions or deletions, and chucking the sequence if they are no linear matches
	if ($cigar =~ /[ID]/){
	    # warn "ignoring sequence:\n $_\n"; sleep(1);
	}
	
	### extracting the genomic sequence instead
 	my $extracted_sequence = substr($chromosomes{$chr},$start - 1,length$sequence);
	# warn "Old sequence: $sequence\nExt sequence: $extracted_sequence\n";
	
	if ($single){
	    calc_single_end($extracted_sequence,$flag);
	}
	else{
	    calc_paired_end($extracted_sequence,$flag);
	}
	
    }

    foreach my $dinuc(sort keys %nuc){
	next if ($dinuc =~/N/);
	print OUT "$dinuc\t$nuc{$dinuc}\n";
    } 

}

sub calc_paired_end{
    my ($sequence,$flag) = @_;
    # warn "FLAG: $flag\n$sequence\n"; sleep(1);    
    
    if ($flag == 99){ # OT or CTOT
	# warn "flag 99. Don't need to do anything\n"; # fine, don't need to do anything
    }
    elsif ($flag == 83){ # OB or CTOB
	# reverse complementing
	$sequence = reverse $sequence;
	$sequence =~ tr/GATC/CTAG/;
    }
    elsif ($flag == 147 or $flag == 163){
	# warn "Read 2, skipping\n";
    }
    else{
	# die "failed to detect flag\n";
    }
    
    foreach my $start (0..(length$sequence)-2){ 
	my $dinuc = substr($sequence,$start,2);
	$nuc{$dinuc}++;
    }
    # warn "$sequence\n~~~~~\n";
}

sub calc_single_end{
    my ($sequence,$flag) = @_;
    
    if ($flag == 0){
	# warn "flag 0\n"; # fine, don't need to do anything
    }
    elsif ($flag == 16){
	# reverse complementing
	$sequence = reverse $sequence;
	$sequence =~ tr/GATC/CTAG/;
    }
    else{
	die "failed to detect flag\n";
    }
    
    foreach my $start (0..(length$sequence)-2){ 
	my $dinuc = substr($sequence,$start,2);
	$nuc{$dinuc}++;
    }
    # warn "$sequence\n~~~~~\n";
    #  sleep(4);    
}


sub test_file{
    my $in = shift;   
    open (TEST,"samtools view -H $in |") or die $!;
    
    while (<TEST>){
	chomp;
	if (/^\@PG/) {
	    if (/ -1 / and / -2 /){ # paired-end
		close TEST or warn $!;
		return (0,1);
	    }
	    else{ # single-end
		close TEST or warn $!;
		return (1,0);	
	    }
	}
    }

}

sub print_all_fastA_sequences{
    my $printed = 0;
    foreach my $key(keys %positions){
	#   print ">Chr$positions{$key}->{chrom}:$positions{$key}->{start}-$positions{$key}->{end}\n";
	my $sequence = substr ($chromosomes{$positions{$key}->{chrom}},$positions{$key}->{start},($positions{$key}->{end}-$positions{$key}->{start}+1));
	#  print "$sequence\n";
	++$printed;
	sleep (1);
    }
    warn "printed $printed positions\n";
}



#sub read_genome_into_memory{
#  ### working directoy
#  my $cwd = getcwd();
#  ### reading in and storing the mouse genome in the %chromosomes hash
#  # chdir ("/bi/scratch/Genomes/Chinese_Hamster_PreEnsembl/") or die "Can't move to folder /bi/scratch/Genomes/Chinese_Hamster_PreEnsembl: $!";
#  chdir ("/bi/scratch/Genomes/Chinese_Hamster_PreEnsembl/") or die; 
#  warn "Now reading in and storing sequence information of the hamster genome (build Chinese_Hamster_PreEnsembl)\n\n";
#  while (my $chromosome_filename = <*.fa>){
 #     my $chromosome_number = chromosome_number($chromosome_filename);
 #     my $sequence = read_chromosomal_sequence($chromosome_filename);
 #     $chromosomes{$chromosome_number}= $sequence;
 #     warn "chromosome $chromosome_number\n";
 # }
 # warn "\n";
 # chdir $cwd or die "Failed to move to directory $cwd\n";
 # warn "Moved back to working directory $cwd\n\n";
#}

#sub read_chromosomal_sequence{
#  my $filename = shift @_;
#  my $sequence;
#  # warn "Reading sequence data from $filename\n";
#  open (CHROMOSOME,$filename)or die "Can't open $filename: $!";
#  $_ = <CHROMOSOME>;
#  while (<CHROMOSOME>){
#    chomp;
#    $sequence .= uc$_;
#  }
#  close CHROMOSOME or die "Failed to close filehandle\n";
#  return $sequence;
#}

#sub chromosome_number{
#    my $filename = shift @_;
#    if ($filename =~ /\.([^\.]+)\.fa$/){
#warn "extracted chromosome name $1\n"; sleep(1);
#	return $1;
 #   }
  #  else{
#	die "Unable to extract the chromosome number: $filename!";
 #   }
#



sub read_genome_into_memory{
    ## working directoy
    my $cwd = shift;
    ## reading in and storing the specified genome in the %chromosomes hash
    chdir ($genome_folder) or die "Can't move to $genome_folder: $!";
    warn "Now reading in and storing sequence information of the genome specified in: $genome_folder\n\n";

    my @chromosome_filenames =  <*.fa>;

    ### if there aren't any genomic files with the extension .fa we will look for files with the extension .fasta
    unless (@chromosome_filenames){
      @chromosome_filenames =  <*.fasta>;
    }

    unless (@chromosome_filenames){
	die "The specified genome folder $genome_folder does not contain any sequence files in FastA format (with .fa or .fasta file extensions)\n";
    }
    
    foreach my $chromosome_filename (@chromosome_filenames){

	open (CHR_IN,$chromosome_filename) or die "Failed to read from sequence file $chromosome_filename $!\n";
	### first line needs to be a fastA header
	my $first_line = <CHR_IN>;
	chomp $first_line;
	$first_line =~ s/\r//;
	### Extracting chromosome name from the FastA header
	my $chromosome_name = extract_chromosome_name($first_line);
	my $sequence;

	while (<CHR_IN>){
	  chomp;
	  $_ =~ s/\r//; # removing carriage returns if present
	  if ($_ =~ /^>/){
	
	    ### storing the previous chromosome in the %chromosomes hash, only relevant for Multi-Fasta-Files (MFA)
	    if (exists $chromosomes{$chromosome_name}){
	      print "chr $chromosome_name (",length $sequence ," bp)\n";
	      die "Exiting because chromosome name already exists. Please make sure all chromosomes have a unique name!\n";
	    }
	    else {
	      if (length($sequence) == 0){
		  warn "Chromosome $chromosome_name in the multi-fasta file $chromosome_filename did not contain any sequence information!\n";
	      }
	      print "chr $chromosome_name (",length $sequence ," bp)\n";
	      $chromosomes{$chromosome_name} = $sequence;
	    }
	    ### resetting the sequence variable
	    $sequence = '';
	    ### setting new chromosome name
	    $chromosome_name = extract_chromosome_name($_);
	  }
	  else{
	    $sequence .= uc$_;
	  }
	}
	
 	### Processing last chromosome of a multi Fasta File or the only entry in case of single entry FastA files

	if (exists $chromosomes{$chromosome_name}){
	    print "chr $chromosome_name (",length $sequence ," bp)\t";
	    die "Exiting because chromosome name already exists. Please make sure all chromosomes have a unique name.\n";
	}
	else{
	    if (length($sequence) == 0){
		warn "Chromosome $chromosome_name in the file $chromosome_filename did not contain any sequence information!\n";
	    }

	    print "chr $chromosome_name (",length $sequence ," bp)\n";
	    $chromosomes{$chromosome_name} = $sequence;
	}
    }
    print "\n";
    chdir $cwd or die "Failed to move to directory $cwd\n";
}

sub extract_chromosome_name {
    ## Bowtie seems to extract the first string after the inition > in the FASTA file, so we are doing this as well
    my $fasta_header = shift;
    if ($fasta_header =~ s/^>//){
	my ($chromosome_name) = split (/\s+/,$fasta_header);
	return $chromosome_name;
    }
    else{
	die "The specified chromosome ($fasta_header) file doesn't seem to be in FASTA format as required!\n";
    }
}
