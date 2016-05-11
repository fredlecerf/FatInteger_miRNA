#!/usr/bin/perl -w
#-------------------------------------------------------------------------------
# Name: miR_03_Joint_Analyzis
#-------------------------------------------------------------------------------
# Desc: Analyze results from TargetScan and miranda together and provide
# the common (or not) results)
#
# 1.00: initial release.
#-------------------------------------------------------------------------------
# Author: Frédéric lecerf on May., 2016
#-------------------------------------------------------------------------------
use strict;
use Getopt::Long;


#-------------------------------------------------------------------------------
# CONSTANT value definition
#-------------------------------------------------------------------------------
$|=1;
my $soft = "miR_03_Joint_Analyzis";
my $VERSION = "1.00";
my $year = "2016";

my $target_file = undef;
my $miranda_file = undef;

# Benchmark 
my $begin_time = times();


my %Options = (
        'target=s' => \$target_file,
		'miranda=s' 	=> \$miranda_file,
	      );
my %OptionsHelp = (
		    'Miranda results   ' => '-m [file], miranda results',
			'TargetScan results' => '-t [file], TargetScan results',
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
defined $miranda_file 	|| &usage('Miranda results file (-m) required');
defined $target_file	|| &usage('Targetscan results file (-s) required');


print "--------------------------------------------------------------------------------\n";
print "$soft $VERSION ($year)\n";
print "--------------------------------------------------------------------------------\n";


#-------------------------------------------------------------------------------
# Read TargetScan results file
#-------------------------------------------------------------------------------
open (IN, $target_file)	|| die "cannot open file $target_file";

print "Reading Targetscan results: $target_file\n";

my %Target = ();

while (my $line = <IN>) {
	$line =~ s/\s+$//;
	next if ($line =~ /^#/);
	
	my @T = split(/\t/, $line);
	
	my $miR_Family = $T[0];
	my $transcriptID = $T[1];
	my $genename = $T[2];
	my $siteType = $T[5];
	my $contextScore = $T[6];
	my $CSPercent = $T[8];
	
	# create a list to manipulate whole data easier
	# order is : ENST SiteType ContextScore CSPercent
	my @String = ("$transcriptID", "$siteType", "$contextScore", "$CSPercent");
	
	# Test first if entry exists (to filter redundant results)
	# If yes, keep only the best contextscore (the lowest)
	if (!exists $Target{$miR_Family}{$genename}) {
		$Target{$miR_Family}{$genename}{transcriptID} = $transcriptID;
		$Target{$miR_Family}{$genename}{CS} = $contextScore;
		$Target{$miR_Family}{$genename}{string} = \@String;
	}
	else {
		if ($contextScore < $Target{$miR_Family}{$genename}{CS}) {
			$Target{$miR_Family}{$genename}{transcriptID} = $transcriptID;
			$Target{$miR_Family}{$genename}{CS} = $contextScore;
			$Target{$miR_Family}{$genename}{string} = \@String;
		}
	}
}

close (IN);


#-------------------------------------------------------------------------------
# Read Miranda results file
#-------------------------------------------------------------------------------
open (IN, $miranda_file)	|| die "cannot open file $miranda_file";

print "Reading Miranda results: $miranda_file\n";

my %Miranda = ();

while (my $line = <IN>) {
	$line =~ s/\s+$//;
	next if ($line =~ /^#/);
	
	my @T = split(/\t/, $line);
	
	my $miR = $T[0];
	my $ensID = $T[1];
	my $transcriptID = $T[2];
	my $genename = $T[3];
	my $score = $T[4];
	my $energy = $T[5];
	
	# order is : ENSG ENST Score Energy
	my @String = ("$ensID", "$transcriptID", "$score", "$energy");
	
	if (!exists $Miranda{$miR}{$genename}) {
		$Miranda{$miR}{$genename}{energy} = $energy;
		$Miranda{$miR}{$genename}{string} = \@String;
	}
	else {
		if ($energy < $Miranda{$miR}{$genename}{energy}) {
		$Miranda{$miR}{$genename}{energy} = $energy;
		$Miranda{$miR}{$genename}{string} = \@String;
		}
	}
}



#-------------------------------------------------------------------------------
# Generate FINAL results
#-------------------------------------------------------------------------------


open (OUT, ">MixedRESULTS-TargetScan_Miranda") or die 'cannot create output file!';


print OUT "#miR_Family\tHGNC\tts-ENST\tts-SiteType\tts-ContextScore\tts-CSPercent\tm-ENSGALG\tm-ENSGALT\tm-Score\tm-Energy\n";
print "Generating joint results output\n";
foreach my $miR (keys %Miranda) {
	print "... processing with miR: $miR\n";
	foreach my $genename (keys %{$Miranda{$miR}}) {
		if (exists $Target{$miR}{$genename}) {
			print OUT "$miR\t$genename\t".join("\t", @{$Target{$miR}{$genename}{string}})."\t".join("\t", @{$Miranda{$miR}{$genename}{string}})."\n";
		}
		
		
	}
}

close (OUT);

print "Results written to file: MixedRESULTS-TargetScan_Miranda\n";


my $end_time=times();





print "\n--\n";
printf "Execution time: %.2f seconds CPU user-time\n",($end_time-$begin_time);