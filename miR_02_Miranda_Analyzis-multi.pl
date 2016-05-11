#!/usr/bin/perl -w
#-------------------------------------------------------------------------------
# Name: miR_02_Miranda_Analyzis
#-------------------------------------------------------------------------------
# Desc: Parse the FilteredMiranda result and add the GENENAME (instead of the
# Ensembl ID).
#
# 1.00: initial release.
#-------------------------------------------------------------------------------
# Author: Frédéric lecerf on May., 2016
#-------------------------------------------------------------------------------
use strict;
use Getopt::Long;
use Error qw(:try);
use Bio::FeatureIO;
use Bio::SeqFeature::Generic;
use Bio::Tools::GFF;
use File::Util;

#-------------------------------------------------------------------------------
# CONSTANT value definition
#-------------------------------------------------------------------------------
$|=1;
my $soft = "miR_02_Miranda_Analyzis.pl";
my $VERSION = "1.00";
my $year = "2016";

my $file_gtf = undef;

# Benchmark 
my $begin_time = times();


my %Options = (
		'gtf=s' 	=> \$file_gtf,
	      );
my %OptionsHelp = (
			'gtf            ' => '-g [file], GTF annotation file from Ensembl',
		    );

#-------------------------------------------------------------------------------
sub usage ( $ ) {
	my ($msg)=@_;
	print STDERR "Error: $msg\n";
	print "--------------------------------------------------------------------------------\n";
	print "$soft $VERSION ($year)\n";
	print "--------------------------------------------------------------------------------\n";
	print "\n";
	print "--------------------------------------------------------------------------------\n";
	print STDERR "Usage: $0 [options]\n";
	print STDERR "see README for specific configuration parameters\n";
	print STDERR "Options:\n";
	map {printf STDERR "\t$_: %s\n",$OptionsHelp{$_};} keys %OptionsHelp; 
	print STDERR "Please cite: Lecerf F., $year\n";
	exit(1);
}


#-------------------------------------------------------------------------------
# Check parameters
#-------------------------------------------------------------------------------
GetOptions(%Options);
defined $file_gtf	|| &usage('GTF file required (-g)');


print "--------------------------------------------------------------------------------\n";
print "$soft $VERSION ($year)\n";
print "--------------------------------------------------------------------------------\n";



#-------------------------------------------------------------------------------
# Loading data from GTF annotation file using GTF BioPerl module
#-------------------------------------------------------------------------------

my $parser = new Bio::Tools::GFF->new(-file=> $file_gtf, -gff_version => 2);

my %Annotation = ();

print "Loading data from Ensemble GTF annotation file: $file_gtf\n";
while( my $result = $parser->next_feature ) {
	#print "############################################################################################################################\n";
	try {
		my @gene_id       = $result->get_tag_values("gene_id");
		my @transcript_id = $result->get_tag_values("transcript_id");
		my @gene_name     = $result->get_tag_values("gene_name");
		
		#print "ID: $gene_id[0] Tid: $transcript_id[0] gene_name: $gene_name[0]\n";
		$Annotation{$transcript_id[0]} = $gene_name[0];
	}
	catch Bio::Root::Exception with {
		my $err = shift;
		if ($err =~ /gene_name/) {
			
		}
		else {
			#print "A Bioperl exception occurred:\n$err\n";
		}
	};
}

print "... done\n";
print "... ".(keys %Annotation)." annotation loaded\n";


#-------------------------------------------------------------------------------
# Read results file
#-------------------------------------------------------------------------------

my ($f) = File::Util->new();
my (@Files) = $f->list_dir('.', qw /--files-only --pattern=\.filteredMiranda$/);

open (OUT, ">Miranda_PROCESSED") or die 'cannot create output file!';

print OUT "#miR_Family\tGeneID\tTranscriptID\tHGNC\tScore\tEnergy\n";

foreach my $file_result (@Files) {
	open (IN,$file_result)	|| die "cannot open file $file_result";
	
	print "Reading Miranda results: $file_result\n";
	
	my %ResultsTargetScan = ();
	my $nb = 0;
	

	while (my $line = <IN>) {
		$line =~ s/\s+$//;
		next if ($line =~ /^#/);
		
		my @T = split(/\t/, $line);
		
		my $miR = $T[0];
		$miR =~ s/>gga-//;
		$miR =~ s/>hsa-//;
		$miR =~ s/>mmu-//;
		
		my @GeneID = split (/\|/,$T[1]);
		my $ensID = $GeneID[0];
		my $transcriptID = $GeneID[1];
		
		my $score = $T[2];
		my $energy = $T[3];
		
		my $genename = 'NA';
		$genename = $Annotation{$transcriptID} if (exists $Annotation{$transcriptID});
		
		print OUT "$miR\t$ensID\t$transcriptID\t$genename\t$score\t$energy\n";
	}
	close (IN);	
}

close (OUT);

print "Results written to file: Miranda_PROCESSED\n";


my $end_time=times();





print "\n--\n";
printf "Execution time: %.2f seconds CPU user-time\n",($end_time-$begin_time);