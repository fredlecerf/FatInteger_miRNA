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
my $VERSION = "1.01";
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
	#next if (($siteType ne '7mer-m8')||($siteType ne '8mer-1a'));
	
	my $genename = 'NA';
	$genename = $Annotation{$transcriptID} if (exists $Annotation{$transcriptID});
	
	print OUT "$miR_Family\t$transcriptID\t$genename\t$speciesID\t$mirID\t$siteType\t$contextScore\t$contextScoreWeighted\t$contextScorePercent\n";
	#$ResultsTargetScan{$miR_Family}{$geneID}{$speciesID}{$siteType}{mirID} = $mirID;
	#$ResultsTargetScan{$miR_Family}{$geneID}{$speciesID}{$siteType}{CS} = $contextScore;
	#$ResultsTargetScan{$miR_Family}{$geneID}{$speciesID}{$siteType}{CSweight} = $contextScoreWeighted;
	#$ResultsTargetScan{$miR_Family}{$geneID}{$speciesID}{$siteType}{CDPercent} = $contextScorePercent;
	#$nb++;
	#print "$line\n" if ($geneID eq 'ENST00000370440.1');
}

#print "... $nb entries loaded\n";

close (IN);
close (OUT);


exit(1);

foreach my $fam (keys %ResultsTargetScan) {
	foreach my $gene (keys %{$ResultsTargetScan{$fam}}) {
		next if ($gene ne 'ENST00000370440.1');
		foreach my $species (keys %{$ResultsTargetScan{$fam}{$gene}}) {
			print "$fam - $gene - $species\n";
		}
	}
}

#-------------------------------------------------------------------------------
# Analyze results file
#-------------------------------------------------------------------------------

print "\nAnalyze TargetScan prediction file\n";

foreach my $miR_family (keys %Input_miR) {
	print "Dealing with $miR_family\n";
	my $nb_species = scalar @{$Input_miR{$miR_family}};
	print "... $nb_species species loaded for this miR\n";
	
	# multi-species analysis mode
	if ($nb_species > 1) {
		print "... switching to multispecies analyzis mode\n";
		foreach my $geneID (keys %{$ResultsTargetScan{$miR_family}}) {
			my $nb_species_inHash = keys %{$ResultsTargetScan{$miR_family}{$geneID}};
			
			# First, test if only one species and if equal to chicken
			# if many species, check if chicken is included in list
			my $onespecies = undef;
			my $multispecies = undef;
			my $chickenfound = undef;
			
			if ($nb_species_inHash == 1) {
				$onespecies = defined;
			}
			else {
				$multispecies = defined;
			}
			
			
			foreach my $species (keys %{$ResultsTargetScan{$miR_family}{$geneID}}) {
				#if ($miR_family eq 'miR-451a') {
				print " >>> nb species: $nb_species_inHash\n";
				print " >>> $species\n";
				#}
				$chickenfound = defined if ($species eq '9031');
				
			}
			$onespecies = 'yes' if (defined $onespecies);
			$onespecies = 'no' if (!defined $onespecies);
			$multispecies = 'yes' if (defined $multispecies);
			$multispecies = 'no' if (!defined $multispecies);
			$chickenfound = 'yes' if (defined $chickenfound);
			$chickenfound = 'no' if (!defined $chickenfound);
			
			print "$geneID = one species: $onespecies multispecies: $multispecies GGA found: $chickenfound\n";
			
		}
	}
	# only species
	else {
		print "... switching to unique species analyzis mode\n";
	}
	
}



my $end_time=times();





print "\n--\n";
printf "Execution time: %.2f seconds CPU user-time\n",($end_time-$begin_time);