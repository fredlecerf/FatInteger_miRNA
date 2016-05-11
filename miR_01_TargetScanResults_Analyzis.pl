#!/usr/bin/perl -w
#-------------------------------------------------------------------------------
# Name: miR_01_TargetScanResults_Analyzis
#-------------------------------------------------------------------------------
# Desc: Analyze the Targetscan results based on the miR input file. The strategy
# of the analyses is the following :
#	- for miR specific to only one species, the script gets the results and
#	  sorts them based on context score and match site type,
#	- for multi-species miR, the selection is based on the presence of a miR in
#	  GGA species and another one (+ the criterion used for species specific
#	  miR).
#
# 1.00: initial release
# 1.01: the multispecies analysis is irrelevant with GGA since, there's only
# results for this species (even if there are many species in input : 2188).
# for mir-451a, there's no GGA in results.
# >> Code 'hack' to filter only the results.
# 1.02: add support for Ensemble annnotation file (GTF)
#-------------------------------------------------------------------------------
# Author: Frédéric lecerf on May., 2016
#-------------------------------------------------------------------------------
use strict;
use Getopt::Long;
use Error qw(:try);
use Bio::FeatureIO;
use Bio::SeqFeature::Generic;
use Bio::Tools::GFF;

#-------------------------------------------------------------------------------
# CONSTANT value definition
#-------------------------------------------------------------------------------
$|=1;
my $soft = "miR_01_TargetScanResults_Analyzis.pl";
my $VERSION = "1.02";
my $year = "2016";

my $file_result = undef;
my $file_miR = undef;
my $verbose = undef;
my $file_gtf = undef;

# Benchmark 
my $begin_time = times();


my %Options = (
        'results=s' => \$file_result,
		'input=s' 	=> \$file_miR,
		'gtf=s' 	=> \$file_gtf,
		'verbose' 	=> \$verbose,
	      );
my %OptionsHelp = (
		    'TargetScan results' => '-r [file]',
			'miR input file    ' => '-i [file]',
			'gtf               ' => '-g [file], GTF annotation file from Ensembl',
			'verbose           ' => '-v, verbose (optionnal)',
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
defined $file_result || &usage('TargetScan results file (-r) required');
defined $file_miR	 || &usage('miR input file (-i) required');
defined $file_gtf	|| &usage('GTF file required (-g)');


print "--------------------------------------------------------------------------------\n";
print "$soft $VERSION ($year)\n";
print "--------------------------------------------------------------------------------\n";


#-------------------------------------------------------------------------------
# Read miR input file
#-------------------------------------------------------------------------------
open (IN,$file_miR)	|| die "cannot open file $file_miR";

print "Reading targetscan input file: $file_miR\n";

my %Input_miR = ();
my $nb_miR = 0;

while (my $line = <IN>) {
	$line =~ s/\s+$//;
	my @T = split (/\t/, $line);
	my $miRfamily = $T[0];
	my @SpeciesID = split (';', $T[2]);
	
	$Input_miR{$miRfamily} = \@SpeciesID;
	$nb_miR++;
}

print "... $nb_miR miR loaded\n";

close (IN);


#-------------------------------------------------------------------------------
# Loading data from GTF annotation file using GTF BioPerl module
#-------------------------------------------------------------------------------

my $parser = new Bio::Tools::GFF->new(-file=> $file_gtf, -gff_version => 2);

my %Annotation = ();

print "Loading data from Ensemble GTF annotation file: $file_gtf\n";
while( my $result = $parser->next_feature ) {
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
open (IN,$file_result)	|| die "cannot open file $file_result";

print "Reading targetscan results: $file_result\n";

my %ResultsTargetScan = ();
my $nb = 0;

open(OUT, ">RESULTS_TargetScan") or die 'cannot create output file!';

print OUT "#miR_Family\tTranscriptID\tHGNC\tSpecies\tmiR_ID\tSite_Type\tContextScore\tContextScoreWeighted\tContextScorePercent\n";
while (my $line = <IN>) {
	$line =~ s/\s+$//;
	next if ($line =~ /^#/);
	
	my @T = split(/\t/, $line);
	
	my @GeneID = split (/\./,$T[0]);
	my $transcriptID = $GeneID[0];
	
	my $speciesID = $T[1];
	my $mirID = $T[2];
	my $siteType = $T[3];
	my $contextScore = $T[27];
	my $contextScoreWeighted = $T[30];
	my $contextScorePercent = $T[28];
	my $miR_Family = $T[35];
	
	next if ($speciesID ne '9031');
	
	my $genename = 'NA';
	$genename = $Annotation{$transcriptID} if (exists $Annotation{$transcriptID});
	
	print OUT "$miR_Family\t$transcriptID\t$genename\t$speciesID\t$mirID\t$siteType\t$contextScore\t$contextScoreWeighted\t$contextScorePercent\n";
}

print "Results written to file: RESULTS_TargetScan\n";

close (IN);
close (OUT);



my $end_time=times();





print "\n--\n";
printf "Execution time: %.2f seconds CPU user-time\n",($end_time-$begin_time);