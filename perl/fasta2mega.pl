#!/usr/bin/perl

my $usage="Usage: $0 [-h] [-f paml|fasta] [inputFile]\n" .
  "  -f: input format. paml, or fasta.  If no flag, fasta is default\n" .
  " STDIN is used as the input if no fastaFile is given\n";

my $sep = "\t";  # if you use tab in the sequence name, change this to
                 # other characters such as ","
my $debug = 0;

use Getopt::Std;
getopts('hf:') || die "$usage\n";
die "$usage\n" if (defined($opt_h));

@ARGV = ('-') unless @ARGV; # take STDIN when no arg.

## read in seq data
my $seqFile = shift @ARGV;
my @dat = ();

if (defined($opt_f) && lc($opt_f) eq 'paml') {
    @dat = ReadInPAML($seqFile);
} else {
    @dat = ReadInFASTA($seqFile);
}

PrintMega(@dat);

#@dat = AdjustSeqLength(@dat);    # attach '?' for shorter sequences
#@dat = RemoveGapOnlySites(@dat); # @partitionSize is modified in here

exit(0);


# takes an arg; name of a file from which data are read Then read in
# the data and make an array.  Each element of this array corresponds
# to a sequence, name tab data.
sub ReadInFASTA {
    my $infile = shift;
    my @line;
    my $i = -1;
    my @result = ();
    my @seqName = ();
    my @seqDat = ();

    open (INFILE, "<$infile") || die "Can't open $infile\n";

    while (<INFILE>) {
        chomp;
        if (/^>/) {  # name line in fasta format
            $i++;
            s/^>\s*//; s/^\s+//; s/\s+$//;
            $seqName[$i] = $_;
            $seqDat[$i] = "";
        } else {
            s/^\s+//; s/\s+$//;
	    s/\s+//g;                  # get rid of any spaces
            next if (/^$/);            # skip empty line
            s/[uU]/T/g;                  # change U to T
            $seqDat[$i] = $seqDat[$i] . uc($_);
        }

	# checking no occurence of internal separator $sep.
	die ("ERROR: \"$sep\" is an internal separator.  Line $. of " .
	     "the input FASTA file contains this charcter. Make sure this " . 
	     "separator character is not used in your data file or modify " .
	     "variable \$sep in this script to some other character.\n")
	    if (/$sep/);

    }
    close(INFILE);

    foreach my $i (0..$#seqName) {
	$result[$i] = $seqName[$i] . $sep . $seqDat[$i];
    }
    return (@result);
}

sub ReadInPAML {
    my $infile = shift;
    my @line;

    open (INFILE, "<$infile") || die "Can't open $infile\n";

    my $header = <INFILE>;
    $header =~ s/^\s+//; $header =~ s/\s+$//;
    my @sizeInfo = split /\s+/, $header;
    my $numSeq = $sizeInfo[0];
    my $lenSeq = $sizeInfo[1];
    
    $multiPartitions = 0;
    if (@sizeInfo != 2) {
	if (@sizeInfo == 3 && $sizeInfo[2] eq "G") {
	    $multiPartitions = 1;
	} else {
	    die "Is this PAML file?\n";
	}
    }
    
    if ($multiPartitions == 1) { # extract partition infos
	$header = <INFILE>;
	my @line = split /\s+/, $header;
	my $sig = shift @line;
	my $numPart = shift @line;
	die "ERROR: 2nd line of input paml file doesn't appear to conatins " .
	    "the paritition info\n" if ($sig ne 'G' || $numPart != @line);

	my $sumPartSize = Sum(@line);
	if (($rUnit eq 'codon' && $lenSeq == $sumPartSize * 3) ||
	    ($rUnit eq 'raw' && $lenSeq == $sumPartSize)) {
	    @partitionSize = @line;
	} elsif ($rUnit eq 'codon' && $lenSeq == $sumPartSize) {
	    warn "WARN: It appears that partition sizes (2nd line in " .
		"in-file) are specified by number of nucleotides.  " .
		    "Adjusted to codon counts\n";
	    @partitionSize = map {$_/3} @line;
	} elsif ($rUnit eq 'raw'  && $lenSeq == $sumPartSize * 3) {
	    warn "WARN: It appears that partition sized (2nd line in " .
		"in-file) are specified by number of codons.  Adjusted to ".
		    "nucleotides counts\n";
	    @partitionSize = map {$_* 3} @line;
	} elsif (! defined ($opt_p)) {
	    die "Partition sizes (2nd line in in-file) are not correct?\n";
	}
    }

    my $cnt=0;
    my @result = ();
    my $i = -1;
    while (<INFILE>) {
	chomp;
	s/\s+$//; s/^\s+//;
	next if (/^$/);
	@line = split /\s+/;
	die "ERROR: line $. doesn't have 2 columns\n$_\n" if (@line != 2);
	
	$i++;
	my $seqName = $line[0];
	$line[1] =~ s/[uU]/T/g;                  # change U to T
	my $seqDat = uc($line[1]);
	my $curSeqLen = CharLen($seqDat);
	if ($curSeqLen != $lenSeq) {
	    warn "Line $.: Header says $lenSeq bases, but $seqName has" .
		" $curSeqLen bases\n";
	}
	push @result, "$seqName\t$seqDat";
    }
    if ($i + 1 != $numSeq) {
	warn "The header says $numSeq sequences, but only ", $i + 1,
	" sequences are found in the input file\n";
    }

    close(INFILE);

    return @result;
}

sub PrintMega{

    my @seqDat =  GetSeqDat(@_);
    my @seqName =  GetSeqName(@_);

    print "#MEGA\n";
    print "!title $seqFile;\n";
    print "!format DataType=DNA indel=- identical=. missing=?;\n";
#    print "!Gene=junk Property=Coding CodonStart=1;\n";
    for my $i (0..$#seqDat) {
	print "#$seqName[$i]\n$seqDat[$i]\n";
    }
}

sub Ck2ArrayRefs {
    my ($vectRef1, $vectRef2) = @_;
    unless (@_ == 2 && ref($vectRef1) eq 'ARRAY' && ref($vectRef2) eq 'ARRAY'){
	die "args to Ck2ArrayRef should be ARRAY REF, ARRAY REF\n";
    }
    my $len = @$vectRef1;
    if ($len != scalar(@$vectRef2) && $len < 1) {
	return -1;
    }
    return ($len);
}

# note this function take only \@seqDat (no names)
# 2nd argument is reference to an array of index (0-offset)
# Note that global $rUnit (unit of resampling) influence the behavior
sub SelectSites {
    my ($seqDatRef, $indexRef) = @_;
    my @result = ();

    for my $seq (@$seqDatRef) {
	my @line = ($rUnit eq 'codon') ? MkTripletArray($seq) : 
	    split (//, $seq);
	@line = SelectElements (\@line, $indexRef);
	my $temp = join ("", @line);
	push @result, $temp;
    }
    return (@result);
}

# for a given string, it separates into triplets, and return
# the resulting array.
# if the last element is less than a triplet, it will be removed
sub MkTripletArray {
    my $seq = shift;
    $seq =~ s/\s+//g;
    $seq =~ s/(.{3})/$1 /g;
    $seq =~ s/\s+$//;
    my @result = split(/ /, $seq);
    pop @result unless ($result[$#result] =~ /.{3}/);
    return @result;
}

# Two array refs as the argument
# 2nd array contains the index of elements to be selected (0-offset)
sub SelectElements {
    my ($arrayRef, $indexRef) = @_;
    unless (@_ == 2 && ref($arrayRef) eq 'ARRAY' && ref($indexRef) eq 'ARRAY'){
	die "args to SelectElements() should be ARRAY REF, ARRAY REF\n";
    }

    my @result = ();
    foreach my $posi (@$indexRef) {
	push @result, $$arrayRef[$posi];
    }
    return @result;
}

sub GetSeqDat {
    my @data = @_;
    my @line;
    my @result = ();

    foreach my $i (@data) {
	@line = split (/$sep/, $i);
	push @result, $line[1];
    }

    return (@result)
}

sub GetSeqName {
    my @data = @_;
    my @line;
    my @result = ();

    foreach my $i (@data) {
	@line = split (/$sep/, $i);
	push @result, $line[0];
    }
    return (@result)
}

sub MaxSeqLen {
    my @data = GetSeqDat(@_);
    my $maxLen = 0;
    foreach $i (@data) {
	my $len = CharLen($i);
	$maxLen = $len if ($len > $maxLen);
    }
    return ($maxLen);
}

# take std seq data (name\tseq), and attach "?" for the shorter sequences
sub AdjustSeqLength {
    my @data = @_;
    my @seqDat = GetSeqDat(@_);
    my @seqName = GetSeqName(@_);
    my $maxLen = MaxSeqLen(@_);
    
    foreach $i (0 .. $#seqDat) {
	my $thisLen = CharLen ($seqDat[$i]);
	if ($thisLen == $maxLen)  {
	    ; # do nothing
	} elsif ($thisLen < $maxLen) {
	    my $diff = $maxLen - $thisLen;
	    warn "WARN: $seqName[$i] shorter.  " .
		"$diff '?' (missing character) were added at the end\n";
	    for ($j=0; $j < $diff; $j++) {
		$data[$i] = $data[$i] . "?";
	    }
	} else {
	    die "ERROR: the length of sequence $seqName[$i] is $thisLen, " .
		"longer than \$maxLen = $maxLen.  Weird!!";
	}
    }
    return (@data);
}

sub RemoveGapOnlySites {
    my @seqDat = GetSeqDat(@_);
    my @seqName = GetSeqName(@_);
    my $maxLen = MaxSeqLen(@_);
    my @gapSites = ();
    my @notGapSites = ();
    my ($posi, $seqNumber);
    my @seqMat = ();

    # make 2 dimensional matrix
    foreach $seqNumber (0..$#seqDat) {
	my @tmpArray = split(//, $seqDat[$seqNumber]);
	# Check the length
	if (@tmpArray != $maxLen)  {
	    die "ERROR: the sequence $seqName[$i] is not same length " .
		"as \$maxLen = $maxLen.  Weird!!";
	}
	push @seqMat, [ @tmpArray ];
    }

    # now identify the all gap sites
    for $posi (0 .. ($maxLen-1)) {
	my $gap = 1;
	for $seqNumber (0 .. $#seqMat){
	    if ($seqMat[$seqNumber][$posi] !~ /^[-\?]$/) {
		$gap = 0;
		last;
	    }
	}
	if ($gap == 1) {  # all sequences have a gap at these sites
	    push (@gapSites, $posi+1); # now unit-offset
	} else {          # there are some non-gap character at these sites
	    push (@notGapSites, $posi);
	}
    }

    # select sites and make results
    my @result = ();
    for $seqNumber (0 .. $#seqMat) {
	my @thisSeq = SelectElements($seqMat[$seqNumber], \@notGapSites);
	my $line = $seqName[$seqNumber] . $sep . (join("", @thisSeq));
	push (@result, $line);
    }

    if (@gapSites > 0) {
	warn ("Following sites consist of all gaps, removed from analysis\n");
	print STDERR join(" ", @gapSites);
	print STDERR "\n";

	# decrements @partitionSize
	AdjustPartitionSize(@gapSites);
    }
    return (@result);
}

# count the number of characters in a string
sub CharLen {
    my $string = shift;
    my @charString = split (//, $string);
    return scalar(@charString);
}

# this function take two scalars and return the larger value
sub Larger {
    my ($a, $b) = @_;

    return (($a > $b) ? $a : $b);
}

sub Sum {
    my @a = @_;
    my $sum = 0;
    foreach my $i (@a) {
	$sum += $i;
    }
    return $sum;
}

sub sortByColumn {
# numerical sort by a column, return an array
#    sortbyColumn ($col_num, $order, @record)
# @record is an array with each element representing a space delimited record
# example
#    ("473 p1 S     0:06 -bash", "541 p2 SW    0:00 ps-a", ....)
# $col_num -- the column by which the record is sorted by (left-most col is 0)
# $order can be "a" (ascending) or "d" (descending),
# sort column can be hyphnated numbers (e.g. 10-4-150)

    local $col_num = shift(@_);
    local $order = shift(@_);
    local @record = @_ ;
    local ($sortCol);
    
    ## test if the sort column is hyphnated or plain number
    local $sortMethod = "number";
    foreach $sortCol (@record) {
	if ( (split(/\s+/,$sortCol))[$col_num] =~ /\d+-\d+/ ) {
	    $sortMethod = "hyphnated";
	    last ;
	}
    }

    return sort $sortMethod @record;

## two sub-sub-routines
    sub number {
	# $col_num, $order are the given arguments
	# the left-most column is 0 
	local $first = (split(/\s+/, $a))[$col_num];
	local $second = (split(/\s+/, $b))[$col_num];
# argument errors not well trapped here
	($first,$second) = ($second, $first) if ($order eq "d");
	
	return ($first <=> $second);
    }

#probably I don't need the "sub number"
    sub hyphnated {
	# $col_num, $order are the given arguments
	local ($each_number, $cmp_result, @temp_swap);

	## separte the hyphnated numbers and put them in the following arrays
        local @first = split(/-/, ((split(/\s+/, $a))[$col_num]));
	local @second = split(/-/, ((split(/\s+/, $b))[$col_num]));

	## ascending (default) or descending order
	if ($order eq "d") {
	    @temp_swap = @first;
	    @first = @second;
	    @second = @temp_swap;
	}
	
	## comparison of array elements
	for ($each_number = 0; $each_number <=
	     (($#first < $#second) ? $#first : $#second) ; $each_number++) {
	    $cmp_result = ($first[$each_number] <=> $second[$each_number]);
	    last if ($cmp_result);
	}

	## if the size of two arrays differ
	if ( ($cmp_result == 0) && ($#first != $#second) ) {
	    return (($#first < $#second) ? -1 : 1);
	} else {
	    return $cmp_result;
	}
    }
}
