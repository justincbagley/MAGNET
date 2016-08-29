#!/usr/bin/perl

# Converts an aligned fasta (aa or dna) seq file to SITES format

my $usage = 
    "Usage: $0 [-hnv] [-b 1|2|3] [-g grouping [-i groupIdentifier]] [-s siteList] [-f siteListFile] [-o offset] [infile]\n" .
    "  -h: help\n" .
    "  -b: beginning of the reading frame (1,2,3), default:1\n".
    "  -g: list of SampleSizes of each group (e.g., 5,3,10), default: no grouping\n".
    "  -i: group Identifier (e.g. 'sp1,sp2,sp3')\n".
    "  -s: specify the site lists of exons, do not put any spaces between elements\n" .
    "     example siteList : 1-4,8,90-\n" .
    "       Default: all sites are considered as non-coding\n".
    "  -f: read site list from a file\n".
    "  -v: The specified sites are actually non-coding (rather than exons)\n".
    "  -o: offset, subtract this value from each listed sites (default = 0)\n".
    "  -n  print name conversion in STDERR\n\n" .
    " Convert FASTA file to Jody Hey's SITES format (similar to phylip format).\n".
    " infile should be an aligned fasta, STDIN is used if no infile is given\n".
    " For -f option, you can use spaces, tab, comma, or new-line delimited\n".
    " site numbers. You can use the range specifiers (see -s), but be\n".
    " CAREFUL not to use spaces around '-'; e.g. 1-4 is ok, but 1 - 4 is not.\n".
    " You can also add comments. For each line, any characters\n".
    " after the 1st # is considered as comments and ignored.\n".
    " Example of site list file:\n\n" .
    "# This file can be given to -f\n".
    "  1-4,     8  # COMMENT: first 5 bases\n".
    "90-\n".
    "# end of file\n";

my $lineWidth = 70;

use IO::File;
use Getopt::Std;
getopts('hvb:g:i:s:f:no:') || die "$usage\n";

die "$usage\n" if (defined ($opt_h));

my $totNumChar = 10;  # number of characters allowed for name in phylip
my $numFrontChar = 7; # When the name is too long, this amount of characters
                      # are used from the beginning of the name, and the rest
                      # are from the end of the name.

my $tmpFH = IO::File->new_tmpfile || die "Unable to make new temp file: $!";

my $firstLine = 1;
my $maxLen = 0;
my $numTaxa = 0;
my $name;

while(<>){
    chomp;
    s/^\s+//; s/\s$//;
    next if (/^$/);

    if (s/^>\s*//) {

	if ($firstLine == 0) {
	    if ($seqLen != $maxLen && $maxLen != 0) {
		warn "WARN: The $numTaxa-th species ($name) have ",
		     "different seq length\n";
		warn "Previous Seq Len: $maxLen, This Seq Len: $seqLen\n";
	    }
	    print $tmpFH "\n";    # end of the previous sequence
	} else {
	    $firstLine = 0;
	}

	$maxLen = $seqLen if ($seqLen > $maxLen); $seqLen = 0;
	$numTaxa ++;

	$name = $_;
	if (CharLen($_) <=10) {
	    printf $tmpFH "%-10s", $_;
	} else  {
	    $shortName = TrimName($_);
	    print STDERR "$_ => $shortName\n" if (defined ($opt_n));
	    printf $tmpFH "%10s", $shortName;
	}
    } else {
	my ($corrected, $notAllowedChar) = CheckAllowedChar($_);
	if ($notAllowedChar ne "")  {
	    warn "Not allowed char converted to N in $name: $notAllowedChar\n";
	}
	# unknown char become 'N'
	$seqLen += CharLen ($corrected);
	print $tmpFH $corrected;
    }
}

print $tmpFH "\n";

# Print comments
print "## converted from FASTA by fasta2sites.pl\n";
# numSeqs and length
print "$numTaxa   $maxLen\n";

## frame (1,2,3) and number of noncoding regions (default to 1)
if (defined($opt_b)) {
    if ($opt_b !~ /^[123]$/) {
	die "ERROR: -b takes either 1, 2, or 3\n";
    }
    print "$opt_b "
} else {
    print "1 ";  #default frame starts from 1
}

my @index = CreateSITESIndex($maxLen);
print scalar(@index), "\n";
foreach my $i (@index) {
    $i =~ s/-/ /;
    print "$i\n";
}

## number of seq groups (default to 1)
my @grpSizes = ($numTaxa);
my $grpNames = ();
if (defined($opt_g)) {
    @grpSizes = split(/,/, $opt_g);
    my $allSamples = Sum(@grpSizes);
    if ($numTaxa != $allSamples) {
	die "ERROR: There are $numTaxa, but your grouping (-g) contains $allSamples.".
	    "  Please check -g.\n";
    }
}

if (defined($opt_i)) { ## grp names <=12 char
    @grpNames = split(/,/,$opt_i);
    my $errFlag = 0;
    foreach my $name (@grpNames) {
	if (CharLen($name) > 12) {
	    warn "ERROR: $name is longer than 12 char\n";
	    $errFlag = 1;
	}
    }
    die "Execution aborted\n" if ($errFlag);
}

my $numNames = scalar(@grpNames);
if (@grpNames < @grpSizes) {
    for (my $i = $numNames + 1; $i <= @grpSizes; $i++) {
	push  @grpNames, "grp$i";
    }
}

print scalar(@grpSizes), "\n";
if (@grpSizes > 1) {
    for (my $i = 0; $i < @grpSizes; $i++) {
	print "$grpNames[$i] $grpSizes[$i]\n";
    }
}

## print the sequence data
seek ($tmpFH, 0, 0) || die "seek: $!";
my $line;
while (defined ($line = $tmpFH->getline())) {
    chomp ($line);
    $missingBases = $maxLen - (CharLen($line) - $totNumChar);
    while ($missingBases > 0) {
	$line = $line . "-";
	$missingBases--;
    }

    print "$line\n";
#    current SITES doesn't like broken lines (contradicting the documentation)
#    for (my $pos=0 ; $pos < length ($line) ;  $pos += $lineWidth) {
#	print substr($line, $pos, $lineWidth), "\n";
#    }
}

sub CharLen {  # returns number of characters in a string
    my $string = shift;
    return scalar(split (//, $string));
}

sub TrimName { # trim a long name
    my $name = shift;
    my @charArray = split (//, $name);
    return join ('', splice (@charArray, 0, $numFrontChar),
		 splice (@charArray, - ($totNumChar - $numFrontChar)));
}

sub CheckAllowedChar {
    my $seqString = shift;
    my $notAllowed = $seqString;

    $notAllowed =~ s/[ATGCNatgcn\.\-\*]//g;
    $seqString =~ s/[^ATGCNatgcn\.\-\*]/N/g;
    return ($seqString, $notAllowed);
}

sub MkSelIndex {
    my ($max, $siteList) = @_;
    $siteList =~ s/^\s+//;
    $siteList =~ s/\s+$//;

    my @sites = split(/\s*,\s*/, $siteList);

    my @result = ();
    foreach my $item (@sites) {
	if ($item =~ /^(\d+)-(\d+)$/) {
	    die "ERROR: 1st number is larger than 2nd in $item\n" if ($1 > $2);
	    $beginPos = $1 - 1;
	    $endPos = $2 - 1;
	} elsif ($item =~ /^-(\d+)$/) {
	    $beginPos = 0;
	    $endPos = $1 - 1;
	} elsif ($item =~ /^(\d+)-$/) {
	    $beginPos = $1 - 1;
	    $endPos = $max-1;
	} elsif ($item =~ /^(\d+)$/) {
	    $beginPos = $1 - 1;
	    $endPos = $1 - 1;
	} else {
	    die "$siteList given as the list of sites.  " . 
		"Make sure it is comma delimitted, and each element is " .
		    " one of the forms: 23-26, 29, -10, 40-\n";  
	}
	
	if (defined ($opt_o)) {  # deal with offset
	    $beginPos -= $opt_o;
	    $endPos -= $opt_o;
	}
	$beginPos = 0 if ($beginPos < 0);
	next if ($endPos < 0);  # ignore the sites < 0
	push (@result, $beginPos..$endPos);
    }
    return(@result);
}

sub MkIndexFromFile {
    my ($max, $file) = @_;
    open (IN,"<$file") || die "Can't open the file $file\n";
    my @result=();
    while(<IN>) {
	chomp();
	s/#.*$//;  # remove comments
	
	s/^\s+//; s/\s+$//;
	next if (/^$/);
	s/,/\t/g;  # convert commans to tab
	s/\s+/\t/g;
	unless(/^[\d\t-]+$/) {
	    warn "This line contains non-numeric, skipped:\n$_\n";
	    next;
	}
	my @line = split;
	push @result, @line;
    }
    my $siteString = join ",", @result;
    @result = MkSelIndex($max, $siteString);
    return @result;
}

sub CreateSITESIndex {
    my $maxLen = shift;

    if (defined ($opt_s) && defined ($opt_f)) {
	die "ERROR: -s and -f can't be used at the same time\n$usage\n";
    }

    my @index = ();
	
    if (defined($opt_s)) {
	@index = MkSelIndex($maxLen, $opt_s);
    } elsif (defined($opt_f)) {
	@index = MkIndexFromFile($maxLen, $opt_f);
    } else {
	@index = ();  # default is that there is no exon.
    }

    # note that, SITES requires index of non-coding region.
    # So need to make the complements without $opt_v
    if (! defined ($opt_v)) {
	my @allSites = 0..($maxLen - 1) ;
	@index = InANotInB (\@allSites, \@index);
    }

    my @tmpIndex = map {$_ + 1} @index;  # convert to 1-offset

    @index = SiteIndexToRange(@tmpIndex);
    return @index;
}

sub SiteIndexToRange {
    my @index = @_;

    return () if (@index == 0);

    @index = sort {$a <=> $b} (@index);
    @index = ExtractUnique (@index);

    $index[$#index+1] = $index[$#index] + 2;  # this force the final index
                                              # to be printed out.
    my $begin = shift @index;
    my $prev = $begin;
    my @results = ();
    for (my $i = 0; $i < @index; $i++) {
	if ($index[$i] == $prev + 1) {
	    $prev = $index[$i];
	    next;
	}
	# now put the range into the results, and reset
	push @results, "$begin-$prev";
	$begin = $prev = $index[$i];
    }

    return (@results);
}

sub InANotInB {
    my ($aRef, $bRef) =@_;
    my %seen = ();
    my @aonly =();

    foreach my $item (@$bRef) { $seen{$item} = 1};
    foreach my $item (@$aRef) {
	push (@aonly, $item) unless $seen{$item};
    }
    return (@aonly);
}

sub ExtractUnique {
    my %seen=();
    my @unique = ();

    foreach my $item (@_) {
        push (@unique, $item) unless $seen{$item}++;
    }
    return @unique;
}

sub Sum {
    $result = 0;
    foreach my $i (@_) {
	$result += $i;
    }
    return $result;
}
