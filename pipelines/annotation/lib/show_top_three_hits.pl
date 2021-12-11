#!/usr/bin/env perl
use strict;
use warnings;

my $last;
my $last_count;
while(<>) {
	my @row = split;
	if (! defined $last || $row[0] ne $last || $last_count < 2) {
		print;
	}
	if (defined $last && $last eq $row[0]) { $last_count++ }
	else { $last_count = 0 }
	$last = $row[0];
}
