#!/usr/bin/perl -w

# This new version assumes that the reads are in sorted order,
# so that we can advance forward-read and reverse-read 
# files in synchrony, even if some are missing.
# (if a reverse is missing, print the current forward and 
# advance that file, or vice versa if a forward is missing).

use strict;
use Carp;

use Getopt::Long;

my $count = 0;
my $help = 0;
my $fwdtag = '.f_';
my $revtag = '.r_';
my $libname = "";
die "Options problem\n" unless
  &GetOptions(
							"fwdtag|f=s"  => \$fwdtag,
							"revtag|r=s"  => \$revtag,
							"libname|n=s" => \$libname,
							"help|h"	    => \$help,
						 );

carp "usage: $0 <args>\n" if $help;

my ($mated, $solo) = (0, 0);
if ($libname) {
	open OUT, ">$libname.mated";
}

# Argument list should be forward-read files

for my $file (@ARGV) {
	my $fwd_file = $file;
	my $rev_file = $fwd_file;
	$rev_file =~ s/$fwdtag/$revtag/;
	print STDERR "FWDS: $fwd_file\n";
	print STDERR "REVS: $rev_file\n";

	unless ($libname) {
		my $out_file = "> $file";
		$out_file =~ s/$fwdtag/_/;
		if ($out_file =~ /\.gz$/) {
			$out_file = "| gzip $out_file";
		}
		open OUT, $out_file;
	}
	# Change input separater to read in whole FASTA records at a time
	$/ = ">";

	if ($fwd_file =~ /^(\S+)\.gz$/) {
		$fwd_file = "gunzip -c $fwd_file |";
	}
	open FWD, $fwd_file;
	if ($rev_file =~ /^(\S+)\.gz$/) {
		$rev_file = "gunzip -c $rev_file |";
	}
	open REV, $rev_file;

	sub parseRead {
	  my $record = shift;
		# print STDERR "+$record+\n";
	  return undef unless $record;
	  my ($num) = ($record =~ /^\s*\S+\.(\d+)/);
		# print STDERR "$num/\n";
	  return $num;
	}
	sub printRead {
		my $record = shift;
		chomp $record;
		return unless ($record and $record =~ /\n[A-Za-z0-9]/); # to skip past empty or summary-only
		$solo++;
		print OUT ">$record";
	}
	sub printReadPair {
		my ($fwd, $rev) = @_;
		chomp $fwd;
		chomp $rev;
		my ($fname, $finfo, $flines) = ($fwd =~ /^\s*(\S+)([^\n]*)\n([[:ascii:]]*)$/);
		my ($rname, $rinfo, $rlines) = ($rev =~ /^\s*(\S+)([^\n]*)\n([[:ascii:]]*)$/);
		if (! $flines or $flines =~ /^#/) {
			&printRead($rev);
		}
		elsif (! $rlines or $rlines =~ /^#/) {
			&printRead($fwd);
		}
		else {
			$fname =~ s/f/m/;
			$rname =~ s/r/m/;
			die "Conflicting names in forward and reverse reads:\n>$fwd\n>$rev\n"
				unless $fname eq $rname;
			$mated++;
			print OUT ">$fname$finfo\n$flines$rlines";
		}
	}

	my $fwd = <FWD>;
	$fwd = <FWD> if $fwd eq '>'; # skip past first delim
	my $rev = <REV>;
	$rev = <REV> if $rev eq '>'; # skip past first delim
	my $fnum = &parseRead($fwd);
	my $rnum = &parseRead($rev);

	READS: while ($fwd or $rev) {
		# If one file is done, run out the other
		if (! $fwd) {
			while ($rev) {
				&printRead($rev);
				$rev = <REV>;
			}
			last READS;
		}
		elsif (! $rev) {
			while ($fwd) {
				&printRead($fwd);
				$fwd = <FWD>;
			}
			last READS;
		}
		# else we have both fwds and revs
		if ($fnum < $rnum) {
			# print forward read & advance forward file
			&printRead($fwd);
			$fnum = &parseRead($fwd = <FWD>);
		}
		elsif ($fnum > $rnum) {
			# print reverse read & advance reverse file
			&printRead($rev);
			$rnum = &parseRead($rev = <REV>);
		}
		else {
			# print together and advance both
			&printReadPair($fwd, $rev);
			$fnum = &parseRead($fwd = <FWD>);
			$rnum = &parseRead($rev = <REV>);
		}
	} # end READS
	unless ($libname) {
		close OUT;
	}
}
print STDERR ("$libname mated: $mated, solo: $solo\n");
