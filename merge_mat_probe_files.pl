#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);
use Data::Dumper;

if ($#ARGV == -1) {
	print STDERR "usage: $0 <fends table> <output mat prefix>  <type>\n";
	exit 1;
}

my $in_fends_fn = $ARGV[0];
my $mat_fn_prefix = $ARGV[1];
my $type = $ARGV[2];


#print STDERR "Input files: ", join(",", @ifns), "\n";
print STDERR "Output prefix: $mat_fn_prefix\n";

##########################################################################################
# read fends file
##########################################################################################

my %fends;

open(IN, $in_fends_fn) || die $in_fends_fn;
print STDERR "Reading input file $in_fends_fn into hash...\n";
my $header = <IN>;
my %h = parse_header($header);
while (my $line = <IN>) {
	chomp $line;
	my @f = split("\t", $line);
	my $fend = $f[$h{fend}];
	my $chr = $f[$h{chr}];
	my $coord = $f[$h{coord}];
	my $strand = $f[$h{strand}];
	my $frag = $f[$h{frag}];
	my $frag_len = $f[$h{frag_len}];
	my $Probe = $f[$h{Probe}]; #Modified by yousra

	!defined($fends{$fend}) or die "non-unique fend";
	$fends{$fend} = {};
	$fends{$fend}->{chr} = $chr;
	$fends{$fend}->{coord} = $coord;
	$fends{$fend}->{Probe}= $Probe; #Modified by yousra
}
close(IN);
######################################################################################################
# Write mat_probe file modified by yousra
#####################################################################################################

 my $in = $mat_fn_prefix.".mat";
 my $outp0= $mat_fn_prefix.".P0.mat";
 my $outp1= $mat_fn_prefix.".P1.mat";
 my $outp2= $mat_fn_prefix.".P2.mat";
 my $p0=0;
 my $p1=0;
 my $p2=0;
 open(IN,$in) || die;
 open(OUT0,">", $outp0) || die;
 print OUT0 "fend1\tfend2\tcount\tprobe1\tprobe2\ttype\n";
 open(OUT1,">", $outp1) || die;
 print OUT1 "fend1\tfend2\tcount\tprobe1\tprobe2\ttype\n";
 open(OUT2,">", $outp2) || die; 
 print OUT2 "fend1\tfend2\tcount\tprobe1\tprobe2\ttype\n";

 my $head = <IN>;
 my %y = parse_header($head);

 print STDERR "Traversing mattable\n";
 my $count = 0;
 while (my $line = <IN>) {
	++$count;
	print STDERR "$count\n" if ($count % 1000000 == 0);
	chomp $line;
	my @f = split("\t", $line);
	my $fend1 = $f[$y{fend1}];
	my $fend2 = $f[$y{fend2}];
	my $count = $f[$y{count}];

	if ($fends{$fend1}->{Probe}  != $fends{$fend2}->{Probe} ) {
		print OUT1 $fend1,"\t",$fend2,"\t",$count,"\t",$fends{$fend1}->{Probe},"\t",$fends{$fend2}->{Probe},"\tP1","\n";
		$p1+=$count;
	} elsif (($fends{$fend1}->{Probe}  == $fends{$fend2}->{Probe}) &&  $fends{$fend2}->{Probe}==0) {	
		print OUT0 $fend1,"\t",$fend2,"\t",$count,"\t",$fends{$fend1}->{Probe},"\t",$fends{$fend2}->{Probe},"\tP0","\n";
		$p0+=$count;
	} elsif (($fends{$fend1}->{Probe}  == $fends{$fend2}->{Probe}) &&  $fends{$fend2}->{Probe}==1) {
		print OUT2 $fend1,"\t",$fend2,"\t",$count,"\t",$fends{$fend1}->{Probe},"\t",$fends{$fend2}->{Probe},"\tP2","\n";
		$p2+=$count;
	}

}

 my $outstats= $mat_fn_prefix.".read_probe.stats";
 open (OUTSTATS,">", $outstats) || die;
 print OUTSTATS "P0\tP1\tP2","\n";
 print OUTSTATS $p0,"\t",$p1,"\t",$p2,"\n";
my $tot = $p0+$p1+$p2;
 print OUTSTATS 100*($p0/$tot),"\t",100*($p1/$tot),"\t",100*($p2/$tot),"\n";
 close(OUT1);
 close(OUT2);
 close(OUT0);
 close(OUTSTATS);

print "Done!\n";

######################################################################################################
# Subroutines
######################################################################################################

# check if C is between and A,B
sub between
{
	my ($A, $B, $C) = @_;
	$B = $B - $A;
	$C = $C - $A;
	return (($B>$C && $C>0) || ($B<$C && $C<0));
}

sub apprx_lines
{
	my ($fn) = @_;
	my $tmp = "/tmp/".getlogin()."_apprx_lines.tmp";
	system("head -n 100000 $fn > $tmp");
	my $size_head = -s $tmp;
	my $size_all = -s $fn;
	$size_head > 0 or die;
	return (int($size_all/$size_head*100000));
}

sub perc_str
{
	my ($n, $total) = @_;
    return ($n." (".(int(1000 * $n / $total)/10)."%)");
}

sub perc_str2
{
	my ($n, $total) = @_;
    return ((int(1000 * $n / $total)/10)."%");
}

# returns first element above/below value in sorted array
sub binary_search 
{
	my $arr = shift;
	my $value = shift;
	my $above = shift; 

	my $left = 0;
	my $right = $#$arr;
	
	while ($left <= $right) {
		my $mid = ($right + $left) >> 1;
		my $c = $arr->[$mid] <=> $value;
		return $mid if ($c == 0);
		if ($c > 0) {
			$right = $mid - 1;
		} else {
			$left  = $mid + 1;
		}
	}
	$left = -1 if ($left > $#$arr);
	$right = -1 if ($right < 0);
	return (($above eq "+") ? $left : $right);
}

sub parse_header
{
	my ($header) = @_;
	chomp($header);
	my @f = split("\t", $header);
	my %result;
	for (my $i = 0; $i <= $#f; $i++) {
		$result{$f[$i]} = $i;
	}
	return %result;
}
