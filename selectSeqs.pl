#!/usr/bin/perl

my $usage="\nSelect specified sequences from a fasta file and print them out\n" .
    "Usage: $0 [-hvp] -f seqNamesFile [fastaFile]\n" .
    "   or  $0 [-hv] -m 'pattern' [fastaFile]\n" .
    "   or  $0 [-hv] -n 'i,j,k,...' [fastaFile]\n" .
    "   or  $0 [-hv] -l 'i-j,...' [fastaFile]\n".
    " -f: specify the filename which contains a list of sequence names\n" .
    "     Each line contains a name. Comments can be added with \"#\".\n" .
    " -p: 'p'attern.  With -f, each line contains perl type regular expression.\n" .
    "     If -p, -v & -f are specified simultaneously, the output includes\n".
    "     sequences which does NOT much any of the regular expressions.\n".
    "     NOTE: '#' should not be a part of the regular expression.\n". 
    " -m 'pattern': extract the sequences whose name contains the pattern. \n".
    "               You can use perl type regular expressions\n".
    " -n 'range': take list of integers.  e.g. -n '-3,5,7-9, 12-' will\n".
    "             select sequences 1,2,3,5,7,8,9,12,13,...\n".
    "             This option can be used to reorder or doubling up some seqs\n".
    "             e.g., -n '3,1,2,2' will bring the 3rd sequence to the front,\n".
    "             and the 2nd sequence get doubled\n".
    "  -l 'size-range': select the sequences based on sequence length.\n".
    "                   only ATGC and degenerate characters (including N) are\n".
    "                   counted, but '-', '?' etc are not counted.\n".
    "                   It takes a list of integers. e.g. -l '500-' will select\n".
    "                   sequences longer than 499bp. -l '100-200,300-400' will\n".
    "                   select sequences in the given ranges (ends are included).\n".
    " -v: Select the sequences which are NOT specified or NOT matching the " .
       "pattern\n".
    "\nIf name of input file (fastaFile) is not given, STDIN is used.\n" .
    "But -f, -m, or -n is mandatory.\n" .
    "Note -vs is not tested well.\n";

# Changelog
# Version 120524
#  - -l option is added
#
# Version 110909
#  - -p option for -f is added
#
# Version 090716
#  - -n option get added.


my $sep = "\t";
my $lineWidth = 70;   # used to break the long sequences into lines.

our ($opt_h, $opt_f, $opt_l, $opt_m, $opt_n, $opt_v, $opt_p);
use Getopt::Std;
getopts('hf:l:m:n:vp') || die "$usage\n";
die "$usage\n" if (defined($opt_h));

if (defined($opt_f) + defined($opt_m) + defined($opt_n) + defined($opt_l) != 1) {
    die "$usage\nERROR: Please specify ONLY 1 of -f, -m, -n, or -l options.\n\n";
}

if (defined ($opt_f)) {
    @selectedSeqs = ReadSeqNameFile($opt_f);
}

die "ERROR: -p can be used only with -f.\n$usage\n" 
    if (defined($opt_p) && (defined($opt_m) || defined($opt_n)));

@ARGV = ('-') unless @ARGV; # take STDIN when no arg.
my $seqFile = shift @ARGV;
my @dat = ReadInFASTA($seqFile);

if (defined ($opt_f)) {
    if (defined($opt_p)) {
	@dat = SelectSeqsByPattern(\@dat, \@selectedSeqs);
    } else { # exact match
	if (defined ($opt_v)) {
	    # this part isn't tested well, so check it.
	    @dat = SelectSeqsNotSpecified(\@dat, \@selectedSeqs);
	} else {
	    @dat = SelectSeqs(\@dat, \@selectedSeqs);
	}
    }
} elsif (defined ($opt_m)) {
    my @tmpPattern = ($opt_m);
    @dat = SelectSeqsByPattern(\@dat, \@tmpPattern);
} elsif (defined ($opt_n)) {
    @dat = SelectSeqsByNumber(\@dat, $opt_n);
} elsif (defined ($opt_l)) {
    @dat = SelectSeqsByLength(\@dat, $opt_l);
}

PrintFASTA(@dat);

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

sub ReadSeqNameFile {
    my $file = shift;
    my @result = ();
    open(INFILE, "<$file") || die "Can't open $file\n";

    while (<INFILE>) {
        chomp;
        s/#.*$//;    # remove comments (#)
        s/^\s+//; s/\s+$//;
        next if (/^$/);
        push @result, $_;
    }
    close(INFILE);
    return @result;
}

sub SelectSeqsNotSpecified {
    my ($seqARef, $nameARef) = @_;
    my @seqName = GetSeqName(@$seqARef);
    my @seqDat = GetSeqDat(@$seqARef);

    my %nameH = ();
    foreach my $n (@$nameARef) {
	$nameH{$n} = 1;
    }
    
    my @result = ();
    for(my $i = 0; $i < @seqName; $i++) {
	next if (defined($nameH{$seqName[$i]}));
	push @result, $seqName[$i] . $sep . $seqDat[$i];
    }
    return (@result);
}

sub SelectSeqs {
    my ($seqARef, $nameARef) = @_;
    my @seqName = GetSeqName(@$seqARef);
    my @seqDat = GetSeqDat(@$seqARef);

    # make a hash table
    my %seqHash = ();
    for my $i (0..$#seqName) {
	if (exists($seqHash{$seqName[$i]})) {
	    die "ERROR: In fasta file, there are more than 1 entry " .
		"which has the name $seqName[$i]\n";
	} else {
	    $seqHash{$seqName[$i]} = $seqDat[$i];
	}
    }

    # select the specified seqs
    foreach my $name (@$nameARef) {
	if (exists($seqHash{$name})) {
	    my $tmp = $name . $sep . $seqHash{$name};
	    push @result, $tmp;
	} else {
	    warn "WARN: $name didn't occur in the input file\n";
	}
    }
    return @result;
}

sub SelectSeqsByLength {
    my ($seqARef, $rangeTxt) = @_;

    my $numSeqs = scalar(@$seqARef);

    my @seqName = GetSeqName(@$seqARef);
    my @seqDat = GetSeqDat(@$seqARef);

    my @lenArr = ();
    foreach my $thisSeq (@seqDat) {
	push @lenArr, SeqLenDegen($thisSeq);
    }

    my $maxLen = Max(@lenArr);

    my @allowedRange = MkSelIndex($maxLen, $rangeTxt); # this is 0-offset
    map {$_ ++} @allowedRange; # 1-offset now

    if (defined($opt_v)) {
	my @allIndex = 1..$maxLen;
	my @complement = InANotInB(\@allIndex, \@allowedRange);
	@allowedRange = @complement;
    }

    @result = ();
    foreach my $i (0..$#seqDat) {
	if (MemberQ($lenArr[$i], \@allowedRange)) {
	    push @result, $seqName[$i] . $sep . $seqDat[$i];
	    warn "INFO: $seqName[$i] = $lenArr[$i] bp. selected\n";
	} else {
	    #warn "INFO: $seqName[$i] = $lenArr[$i] bp. ignored\n";
	}
    }
    
    return @result;
}

# IUPAC
sub SeqLenDegen {
    my $seq = shift;

    if ($aminoAcidSeq) {
	$seq =~ s/[^ACDEFGHIKLMNPQRSTVWYZ]//gi;
    } else {
	$seq =~ s/[^ACGTURYMKSWBDHVN]//gi;
    }

    return length($seq);
}

sub SelectSeqsByNumber {
    my ($seqARef, $rangeTxt) = @_;

    my $numSeqs = scalar(@$seqARef);

    my @index = MkSelIndex($numSeqs, $rangeTxt);

    if (defined($opt_v)) {
	my @allIndex = 0..($numSeqs-1);
	my @complementIndex = InANotInB(\@allIndex, \@index);
	@index = @complementIndex;
    }

    my @result = ();

    foreach my $i (@index) {
	push @result, $$seqARef[$i];
    }
    return @result;
}

# This is from selectSites.pl
# returns 0 offset index
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
	push (@result, $beginPos..$endPos);
    }
    return(@result);
}

sub SelectSeqsByPattern {
    my ($seqARef, $patternARef) = @_;
    my @seqName = GetSeqName(@$seqARef);
    my @seqDat = GetSeqDat(@$seqARef);

    my @result = ();
    my @matchedNames = ();
    foreach my $pattern (@$patternARef) {
	for(my $i = 0; $i < @seqName; $i++) {
	    if ( $seqName[$i] =~ /$pattern/ ) {
		push @matchedNames, $seqName[$i];
	    }
	}
    }

    @matchedNames = ExtractUnique(@matchedNames);

    if (defined($opt_v)) {
	my @complement = InANotInB(\@seqName, \@matchedNames);
	@matchedNames = @complement;
    }

    foreach my $n (@matchedNames) {
# This has a problem if seqName contains reqular expression char
#	my @indexArr = grep {$seqName[$_] =~ /^$n$/} 0..$#seqName;
	my @indexArr = grep {$seqName[$_] eq $n} 0..$#seqName;
	if (@indexArr != 1) {
	    die "ERROR: ". scalar(@indexArr) . " sequences with name: $n\n";
	}
	push @result, $n . $sep . $seqDat[$indexArr[0]];
    }
    
    if (defined($opt_v)) {
	print STDERR "INFO: ". scalar(@result) ." names didn't match the following patterns:\n"; 
    } else {
	print STDERR "INFO: ". scalar(@result) ." names matched the following patterns:\n"; 
    }    
    foreach my $pattern (@$patternARef) {
	print STDERR "  PATTERN: '$pattern'\n";
    }

    print STDERR join("\n", @matchedNames) . "\n\n";
    
    return @result;
}

sub PrintFASTA {
    my @seqName = GetSeqName(@_);
    my @seqDat = GetSeqDat(@_);
    for my $i (0..$#seqDat) {
#	print ">$seqName[$i]\n$seqDat[$i]\n";
	print ">$seqName[$i]\n";
	my $seq = $seqDat[$i];
	for (my $pos=0 ; $pos < length ($seq) ;  $pos += $lineWidth) {
	    print substr($seq, $pos, $lineWidth), "\n";
	}
	
    }
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

sub Max {
    my $max = shift;
    foreach $item (@_) {
        if ($item > $max) {
            $max = $item;
        }
    }
    return $max;
}

sub MemberQ {
    my ($x, $arrRef) = @_;
    foreach my $item (@$arrRef) {
        if ($x eq $item) {
            return 1;
        }
    }
    return 0;
}
