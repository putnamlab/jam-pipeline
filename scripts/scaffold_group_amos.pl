#!/usr/bin/perl -w
# Read connected components from the link graph (contigs connected, transitively,
# by links of reasonably unique orientation)
# Group them greedily, stopping when the next group would put us over 500k 
# total contig length (biggest component for one genome of interest summed to almost 600k)
# Keep maps of contig to component, contig to group
#      Also save a file of components in each group
#           Tab-separated list of contigs per component
#           One component per line
# Scan through Amos-formatted contigs and write each to file for proper group
#      (too many components, ~ 500,000, to put each in a file)

use strict;
use Carp;

use Getopt::Long;

my $components = "";
my $contigs = "";
my $max     = 1.2e6;
my $help = 0;
die "Options problem\n" unless
  &GetOptions(
							"help|h"	         => \$help,
							"components|com=s" => \$components,
							"contigs|con=s"    => \$contigs,
							"maxlen|max|m=i"   => \$max,  # max total length of contigs in a group
						 );
carp "usage: $0 <args>\n" if $help;

sub openstring {
	my $file = shift;
	return (($file =~ /.gz$/)? "gunzip -c $file |" : $file);
}

open COMP, &openstring($components)
	or die "Cannot open components file $components\n";

my $compid = 'aaaaa';
open GROUP, "> g_$compid.components";
my $groupsize = 0;
my $startgroup = $compid;
my %endgroup   = ();

my @cnum2group = ();
my @cnum2comp  = ();
my $total = 0;

while (<COMP>) {
	chomp;
	my $compsize = 0;
	my @contigs  = split;
	my @cnums    = ();
	for my $contig (@contigs) {
		my ($cnum, $len) = ($contig =~ /^(\d+)\[\d+,(\d+)\]$/);
		die "Cannot parse $contig in group $compid, whole list: $_\n"
			unless $len;
		$compsize += $len;
		$cnum2comp[$cnum] = $compid;
		push @cnums, $cnum;
	}
	if ($groupsize and ($groupsize + $compsize > $max)) {
		# Don't add this component to current group, start new one
		print STDERR "Done with group $startgroup ending with $endgroup{$startgroup}, length $groupsize\n";
		$total += $groupsize;
		close GROUP;
		$groupsize = 0;
		$endgroup{$compid} = $startgroup = $compid;
		open GROUP, "> g_$compid.components";
	}
	$groupsize += $compsize;
	$endgroup{$startgroup} = $compid;
	print GROUP "$compid.$compsize\t$_\n";
	for my $cnum (@cnums) {
		$cnum2group[$cnum] = $startgroup;
	}
	$compid++;
}
print STDERR "Done with group $startgroup ending with $endgroup{$startgroup}, length $groupsize\n";
close GROUP;
$total += $groupsize;
print STDERR "Total length of groups: $total\n";

open CONTIGS, &openstring($contigs);

my $tcontigs = "$$.contigs.temp";
open TCONTIGS, "> $tcontigs" 
	or die "Cannot create temporary file $tcontigs for writing\n";

my $contig = -1;
my $group = "";
while (<CONTIGS>) {
	if (/^##(\d+)/) {
		my $newgroup = $cnum2group[$1] ? $cnum2group[$1] : "";
		undef $cnum2comp[$1];
		unless ($newgroup eq $group) {
			$group = $newgroup;
		}
	}
	print TCONTIGS "$group $_" if $group;
}
close TCONTIGS;
system("sort -s -k 1,1 $tcontigs -o $tcontigs");
open TCONTIGS, "$tcontigs"
	or die "Cannot reopen temporary file $tcontigs for reading\n";

$group = "";

while (<TCONTIGS>) {
	s/^(\S+) //;
	if ($1 ne $group) {
		close GROUPCONTIGS if $group;
		$group = $1;
		open GROUPCONTIGS, "> g_$group.amos"
			or die "Cannot create g_$group.amos for writing\n";
	}
	print GROUPCONTIGS;
}
close TCONTIGS;
# unlink $tcontigs;
close GROUPCONTIGS;

# for $contig (0..$#cnum2comp) {
# 	if (defined $cnum2comp[$contig]) {
#		print STDERR "Missing contig $contig in component $cnum2comp[$contig]\n";
#	}
# }
