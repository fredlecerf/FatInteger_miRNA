#!/usr/bin/perl -w
#-------------------------------------------------------------------------------
# Name: miR_00_CleanSequences.pl
#-------------------------------------------------------------------------------
# Desc: 3'UTR sequences fetched from Biomart may be unavailable (no 3'UTR seq).
# In this case, we have to clean the whole FASTA file.
#
#-------------------------------------------------------------------------------
# VERSION HISTORY:
# 1.0 : initial release
#-------------------------------------------------------------------------------
# Author: Frédéric lecerf on february, 2013
#-------------------------------------------------------------------------------
use strict;
use Getopt::Long;

use Bio::Perl;
use Bio::SeqIO;


#-------------------------------------------------------------------------------
# SUB definition
#-------------------------------------------------------------------------------
sub smart_sort ( $ $ );

#-------------------------------------------------------------------------------
# CONSTANT value definition
#-------------------------------------------------------------------------------
$|=1;
my $soft = "miR_00_CleanSequences";
my $VERSION = "1.0";
my $year = "2016";


# Benchmark 
my $begin_time = times();
my $verbose = undef;
my $inputFile = undef;


my %Options = (
        'input=s'	=> \$inputFile,
		'verbose'	=> \$verbose,
	      );

my %OptionsHelp = (
		    'input    ' => '-i [file], fasta file',
			'verbose  ' => '-v, verbose mode (optional)',
		    );


#-------------------------------------------------------------------------------
sub usage ( $ ) {
	my ($msg)=@_;
	print STDERR "Error: $msg\n";
	print "--------------------------------------------------------------------------------\n";
	print "$soft $VERSION ($year)\n";
	print "--------------------------------------------------------------------------------\n";
	print "Desc: 3'UTR sequences fetched from Biomart may be unavailable (no 3'UTR seq).\n";
	print "In this case, we have to clean the whole FASTA file.\n";
	print "see source code description for more details\n";
	print "--------------------------------------------------------------------------------\n";
	print STDERR "Usage: $0 [options]\n";
	print STDERR "see README for specific configuration parameters\n";
	print STDERR "Options:\n";
	map {printf STDERR "\t$_: %s\n",$OptionsHelp{$_};} keys %OptionsHelp; 
	print STDERR "Please cite: Lecerf F., $year\n";
	exit(1);
}


#-------------------------------------------------------------------------------
# Check parameters & read data file
#-------------------------------------------------------------------------------
GetOptions(%Options);
defined $inputFile	|| &usage('Fasta file required (-i)');

print "--------------------------------------------------------------------------------\n";
print "$soft $VERSION ($year)\n";
print "--------------------------------------------------------------------------------\n";


#-------------------------------------------------------------------------------
# Load FASTA file
#-------------------------------------------------------------------------------

my $sequence_stream = Bio::SeqIO->new (-file => "$inputFile" , -format => 'fasta');

my $unavailableSeq = 0;
my $totalSeq = 0;
my $exportSeq = 0;

my $outFileName = "Clean_".$inputFile;

my $stream_out = Bio::SeqIO->new (-format => 'fasta', -file => ">$outFileName");

# parse FASTA file
while (my $seqObject = $sequence_stream->next_seq()) {
	my $header = $seqObject->primary_id();
	my $seq    = $seqObject->seq();
	
	$totalSeq++;
	
	# Test if sequence is marked as 'unavailable'
	if ($seq =~ /^Sequence/) {
		$unavailableSeq++;
		next;
	}
	
	# Else export sequence to new file
	my $seqObjectOut = new_sequence($seq, $header);
	$stream_out->write_seq($seqObjectOut);
	
	$exportSeq++;
}

print "$totalSeq sequence(s) read\n";
print "... $unavailableSeq sequence(s) unavailable\n";
print "... $exportSeq sequence(s) exported to $outFileName\n";


my $end_time=times();
print "\n--\n";
printf "Execution time: %.2f seconds CPU user-time\n",($end_time-$begin_time);




#-------------------------------------------------------------------------------
# "smart" sort : whatever numerical or alpha
#-------------------------------------------------------------------------------
sub smart_sort ( $ $ ) {
	my ($a,$b) = @_;
	# Numerical sort
	return ($a <=> $b) if (($a =~ /^\d+$/)&&($b =~ /^\d+$/)); 
	# Alpha sort
	return ($a cmp $b) if ($a.$b =~ /^\w+$/);
	return ($a cmp $b) if ($a.$b =~ /^\d+\w+$/);
	return ($a cmp $b) if ($a.$b =~ /^\w+\d+$/);
}