#!/usr/bin/perl -w

use strict;
use Carp;

use Getopt::Long;

my $help = 0;
my $loose = 0;  # if set, join components even on single links
die "Options problem\n" unless
  &GetOptions(
							"help|h"	=> \$help,
							"loose|l" => \$loose,
						 );

carp "usage: $0 <args>\n" if $help;

my %orient = ();
my $prev = "";
my %parent = ();
my %rank   = ();

sub find {
  my $x = shift;
  unless ($parent{$x}) {
    $parent{$x} = $x;
    $rank{$x}   = 0;
    return $x;
  }
  if ($x ne $parent{$x}) {
    # path compression
    $parent{$x} = &find($parent{$x});
  }
  return $parent{$x};
}

sub union {
  my ($x, $y) = @_;
  $x = &find($x);
  $y = &find($y);
  return if ($x eq $y);

  if ($rank{$x} > $rank{$y}) {
    $parent{$y} = $x;
  }
  else {
    $parent{$x} = $y;
    if ($rank{$x} == $rank{$y}) {
      $rank{$y}++;
    }
  }
}

sub process {
  my ($pair) = @_;

  my ($best, $othercount, $bestcount) = ("", 0, 0);
  for my $key (keys %orient) {
    if ($orient{$key} > $bestcount) {
      $othercount += $bestcount;
      $best = $key;
      $bestcount = $orient{$key};
    } 
    else {
      $othercount += $orient{$key};
    }
  }
  if ($loose or (($bestcount - $othercount) >= 2 and $bestcount >= 2*$othercount)) { 
    my ($source, $sink) = split " ", $pair;
    # print "$source\t$sink\t$bestcount\n";
    &union($source, $sink);
  }
}

my %components = ();

sub print_components {
  for my $contig (keys %parent) {
    my $par = &find($contig);
    if (defined($components{$par})) {
      push @{$components{$par}}, $contig;
    }
    else {
      $components{$par} = [$contig];
    }
  }
  for my $par (keys %components) {
    print(join("\t", sort(@{$components{$par}})), "\n");
  }
}

sub dump_orients {
  for my $key (keys %orient) {
    print "o:\t$key\t$orient{$key}\n";
  }
}

sub dump_parents {
  print "\t***\n";
  for my $contig (keys %parent) {
    print "\t$contig\t$parent{$contig}\n";
  }
}

while (<>) {
  my @v = split; 

  next unless ($v[0] lt $v[1]); # need only process double links once
  my $pair = "$v[0] $v[1]";
  if ($pair ne $prev) {
    if ($prev) {
      &process($prev);
    }
    $prev = $pair;
    %orient = ();
  }
  $orient{lc($v[2])}++;
}
&process($prev) unless $loose;

# &dump_parents();
&print_components();
