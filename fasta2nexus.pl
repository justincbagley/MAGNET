#!/usr/bin/perl

# Converts an aligned fasta (aa or dna) seq file to phylip format

my $usage = "Usage: $0 [-h] [-v] [infile]\n" .
            "  -h  help\n" .
            "  -a  for amino acid\n" .
            " infile should be an aligned fasta, " .
            "STDIN is used if no infile is given\n";

use IO::File;
use POSIX qw(tmpnam);
use Getopt::Std;

getopts('ah') || die "$usage\n";

die "$usage\n" if (defined ($opt_h));
die "$usage\n" if (@ARGV > 1);

@ARGV = ('-') unless @ARGV;

if (defined ($opt_a)) {
    $type="protein";
} else {
    $type="nucleotide";
}

my $file=shift;
ReadInFASTA($file);

do { $name = tmpnam() }
until $fh = IO::File->new($name, O_RDWR|O_CREAT|O_EXCL);

do { $outname = tmpnam() }
until $outfh = IO::File->new($outname, O_RDWR|O_CREAT|O_EXCL);
close $outfh;

END { unlink($name) or die "Couldn't unlink $name : $!" }
END { unlink($outname) or die "Couldn't unlink $name : $!" }

# find the maximum sequence length
my $maxLen = 0;
foreach my $i (@seqDat) {
    my $len = CharLen($i);
    $maxLen = $len if ($maxLen < $len);
}

# make a temporari TEXT file which PAUP can read.
foreach my $i (0..$#seqName) {
    print $fh "$seqName[$i]\t$seqDat[$i]";
    my $len = CharLen($seqDat[$i]);
    if ($len < $maxLen) {
	warn "WARN $len bases in $seqName[$i], < $maxLen bases: " .
	    "\'?\' (missing) appended to fill up\n";
	for ($j = 0; $j < $maxLen - $len; $j++) {
	    print $fh "?";
	}
    }
    print $fh "\n";
}
close($fh);

# convert to NEXUS and store in a temp file
open (PAUP, "|paup -n > /dev/null 2>&1");
print PAUP "ToNEXUS format=text fromfile=$name datatype=$type tofile=$outname replace=yes; quit";
close (PAUP);

# print out
open (INFILE, "<$outname");
while(<INFILE>) {
    if ($type eq "nucleotide") {
	s/datatype=nucleotide/datatype=dna/;
    }
    print;
}

# takes an arg; name of a file from which data are read
# Then read in the data and store them in @seqName and @seqDat
sub ReadInFASTA {
    my $infile = shift;
    my @line;
    my $i = -1;
    my @result = ();

    open (INFILE, "<$infile") || die "Can't open $infile\n";

    while (<INFILE>) {
        chomp;
        if (/^>/) {  # name line in fasta format
            $i++;
            s/^>\s*//; s/^\s+//; s/\s+$//;
            #@line = split (/\|/);     # note it takes only the name before |
            #$line[0] =~ s/\s+$//;
            $seqName[$i] = $_;
            $seqDat[$i] = "";
        } else {
            s/^\s+//; s/\s+$//;
	    s/\s+//g;
            next if (/^$/);            # skip empty line
            s/[uU]/T/g;                  # change U to T
            $seqDat[$i] = $seqDat[$i] . uc($_);
        }
    }
    close(INFILE);
}



sub CharLen {  # returns number of characters in a string
    my $string = shift;
    return scalar(split (//, $string));
}
