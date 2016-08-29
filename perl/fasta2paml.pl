#!/usr/bin/perl

my $usage =  "Usage: $0 [-h] [-v] [-c numChar] [infile]\n" .
    "  -h: help\n" .
    "  -c: long seq names are shortened to 10 char, default: 7 chars from the\n".
    "      beggining is combined with the last 3 chars.  You can change the\n".
    "      behavior by this option.  E.g., -c 10 will shorten the long name\n" .
    "      by taking the first 10 characters of the name.\n".
    " infile should be an aligned fasta, " .
    "STDIN is used if no infile is given\n";

use Getopt::Std;
getopts('hc:') || die "$usage\n";
die "$usage\n" if (defined($opt_h));

my $args = join " ", @ARGV;

if (defined ($opt_c)) {
    $args = "-c $opt_c $args";
}

if (@ARGV) {
    system("fasta2phylip.pl $args | phylip2paml.pl");
} else {  # input is piped stdin
    open OUT, "|fasta2phylip.pl $args | phylip2paml.pl" || die "Can't open fasta2phylip.pl pipe\n";
    while(<>) {
	print OUT $_;
    }
    close OUT;
}

exit;
