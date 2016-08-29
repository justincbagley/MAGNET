#!/usr/bin/perl

my $firstLine = 1;
while(<>) {
    if ($firstLine == 1) {
	print;
	$firstLine = 0;
	next;
    }
    s/^(.{10})/$1  /;
    print;
}
