#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);

if ($#ARGV == -1) {
	print STDERR "usage: $0 <input file1> <input file2> <output file> <header T or F> <omit edge character from key T|F> <fastq version>\n";
	exit 1;
}

my $ifn1 = $ARGV[0];
my $ifn2 = $ARGV[1];
my $ofn = $ARGV[2];
my $header = $ARGV[3];
my $omit = $ARGV[4] eq "T";
my $fastq_version = $ARGV[5];

my %ids;
my $count = 0;

open(IN, $ifn1) || die "cannot open map $ifn1\n";
print STDERR "Reading input file $ifn1 into hash...\n";
while (my $line = <IN>) {
	$count++;

	chomp $line;
	my @f = split("\t", $line);
	my $id = $f[0];
	my $strand = $f[1];
	my $chr = $f[2];
	my $coord = $f[3];
	my $seq = $f[4];
	
	if ($omit)
	{
		if (index("SXIJ", $fastq_version) >= 0)
		{
			$id =~ s/[0-9]$//;
		}
		elsif ($fastq_version eq "L")
		{
			$id =~ s/ [0-9].*//;
		}
		else
		{
			die "Unknown fastq version: $fastq_version\n";
		}
	}
	$chr =~ s/chr//;

	$ids{$id} = {};

	$ids{$id}->{chr} = $chr;
	$ids{$id}->{coord} = $coord;
	$ids{$id}->{strand} = $strand;
	$ids{$id}->{seq} = $seq;

	print STDERR "line: $count\n" if ($count % 1000000 == 0);
}
close(IN);

print STDERR "line: $count\n";
print STDERR "Reading input file $ifn2 ...\n";
$count = 0;

open(IN, $ifn2) || die "cannot open map $ifn2\n";

# hash to discard pair ended reads that appear more than once
my %reads_found;
my $discarded_reads = 0;

open(OUT, ">", $ofn) || die "cannot open output file $ofn\n";
if ($header eq "T") 
{
	print OUT "chr1\tcoord1\tstrand1\tchr2\tcoord2\tstrand2\n";
}

$count = 0;
while (my $line = <IN>) {
	$count++;

	chomp $line;
	my @f = split("\t", $line);
	my $id = $f[0];
	my $strand = $f[1];
	my $chr = $f[2];
	my $coord = $f[3];
	my $seq = $f[4];

	if ($omit)
	{
		if (index("SXIJ", $fastq_version) >= 0)
		{
			$id =~ s/[0-9]$//;
		}
		elsif ($fastq_version eq "L")
		{
			$id =~ s/ [0-9].*//;
		}
		else
		{
			die "Unknown fastq version: $fastq_version\n";
		}
	}
	$chr =~ s/chr//;

	next if (!defined($ids{$id}));

	my $seq_key = uc($seq." - ".$ids{$id}->{seq});
	if (defined($reads_found{$id})) {
		$reads_found{$id}++;
		$discarded_reads++;
		print "discarded read: $seq_key\n";
		next;
	}
	$reads_found{$id} = 1;

	print OUT $chr, "\t", $coord, "\t", $strand, "\t";
	print OUT $ids{$id}->{chr}, "\t", $ids{$id}->{coord}, "\t", $ids{$id}->{strand}, "\n";

	print STDERR "line: $count\n" if ($count % 1000000 == 0);
}
print STDERR "line: $count\n";
print STDERR "discarded identical pair ended reads: $discarded_reads\n";

close(OUT);
close(IN);
