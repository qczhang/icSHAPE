#! /usr/bin/perl
#
# collapse fastq file to remove PCR duplicates
#   history 20131113
#   version 0.1
#   copyright @cliff zhang
#
use strict;
use warnings;
use File::Basename;
use Getopt::Std;
use FindBin;

my $splitFastqBin = "$FindBin::Bin/../bin/splitFastq";
my $readCollapseBin = "$FindBin::Bin/../bin/readCollapse";

use vars qw ($opt_h $opt_V $opt_D $opt_U $opt_1 $opt_2 $opt_o $opt_p $opt_q $opt_a $opt_f );
&getopts('hVDU:1:2:o:p:q:af:');

my $usage = <<_EOH_;
## --------------------------------------
collapse fastq file to remove PCR duplicates

Command:
$0 -1 fastq_PE_reads_1 -2 fastq_PE_reads_2 -U fastq_SE_reads

# what it is:
 -U     single ends read
 -1     paired ends read 1
 -2     paired ends read 2

# more options:
 -o     single ends read output 
 -p     PE read output 1
 -q     PE read output 2

 -f     unique fasta after collapse

 -a     input is fasta file ( not yet done )

_EOH_
;

&main();

sub main {
    my %parameters = &init();

    my $inFile = $parameters{input1}; my $outFile = $parameters{output1};
    if ( not $parameters{isPairEnds} ) {
        print STDERR "Collapsing file $inFile...\n\t", `date`;
        if ( -e $parameters{output1} ) 
            { print STDERR "Warning! $parameters{output1} exisits, will be overwritten.\n"; print STDERR `/bin/rm $parameters{output1}`; }
        if  ( ( defined $parameters{FASTA} ) and ( -e $parameters{FASTA} ) ) 
            { print STDERR "Warning! $parameters{FASTA} exisits, will be overwritten.\n"; print STDERR `/bin/rm $parameters{FASTA}`; }
    }
    else {
        if ( -e $parameters{output1} ) 
            { print STDERR "Warning! $parameters{output1} exisits, will be overwritten.\n"; print STDERR `/bin/rm $parameters{output1}`; }
        if ( -e $parameters{output2} ) 
            { print STDERR "Warning! $parameters{output2} exisits, will be overwritten.\n"; print STDERR `/bin/rm $parameters{output2}`; }
        if  ( ( defined $parameters{FASTA} ) and ( -e $parameters{FASTA} ) ) 
            { print STDERR "Warning! $parameters{FASTA} exisits, will be overwritten.\n"; print STDERR `/bin/rm $parameters{FASTA}`; }
        $inFile = "/tmp/tmp" . $$ . ".fastq"; $outFile = "/tmp/tmpOut" . $$ . ".fastq";
        _mergePairEndReads ( $parameters{input1}, $parameters{input2}, $inFile );
    }
 
    my @file2Collapse = ( $inFile );
    my ( $cPos, $cLen ) = _estimateSplit ( $inFile );
    if ( $cLen > 0 ) {
        print STDERR "File $inFile too large, will be splitted...\n\t", `date`;
        shift @file2Collapse;
        my $tmpOutDir = "/tmp";
        print STDERR "$splitFastqBin $inFile $tmpOutDir $cPos $cLen new\n\t", `date`;
        my $filesSplited = `$splitFastqBin $inFile $tmpOutDir $cPos $cLen new`;
        if ( $filesSplited =~ /^ERROR/i ) { die ( "Error! splitting file failed. quiting...\n\t" . `date` ); }       ## cleaning may be needed
        else { @file2Collapse = _getSplitedFile ( $filesSplited ); }
    }

    my $totalReads = 0;  my $uniqReads = 0;  my @outputFastas = ();
    foreach my $file ( @file2Collapse ) {
        my $collapseResults = "";
        if ( defined $parameters{FASTA} ) {
            print STDERR "$readCollapseBin $file $outFile append $parameters{FASTA}\n\t", `date`;
            $collapseResults = `$readCollapseBin $file $outFile append $parameters{FASTA}`;
        }
        else {
            print STDERR "$readCollapseBin $file $outFile append\n\t", `date`;
            $collapseResults = `$readCollapseBin $file $outFile append`;
        }
        my ( $volReads, $volUniqReads, $volUniqRatio, $volFasta ) = _parseCollapseOutput ( $collapseResults );
        $totalReads += $volReads;  $uniqReads += $volUniqReads;
        push @outputFastas, $volFasta;
    }
    my $uniqRatio = sprintf ("%.2f", $uniqReads/$totalReads);

    if ( $parameters{isPairEnds} ) { 
        _splitPairEndReads ( $outFile, $parameters{output1}, $parameters{output2} );
        print STDERR `/bin/rm -f $inFile`;
        print STDERR `/bin/rm -f $outFile`;
    }

    if ( not $parameters{isPairEnds} ) { print STDERR "Collapsing file $inFile finished.\n\t", `date`; }
    else { print STDERR "Collapsing file $parameters{input1} and $parameters{input2} finished.\n\t", `date`; }

    print "Read collapse successful! Total count: $totalReads, unique count: $uniqReads, unique ratio: $uniqRatio.\n";
    1;
}

sub init {
    my %parameters = ();

    die $usage if ( $opt_h || ( ( not $opt_U ) && ( ( not $opt_1 ) || ( not $opt_2 ) ) ) );
    $opt_V = 0 if ( not defined $opt_V );
    $opt_D = 0 if ( not defined $opt_D );

    my $pwd = `pwd`;  chomp $pwd;
    if ( defined $opt_U ) {
        $parameters{input1} = $opt_U;
        my ($fileName, $fileDir, $fileSuffix) = fileparse($parameters{input1}, qr/\.[^.]*/);
        if ( defined $opt_o ) { $parameters{output1} = $opt_o; }
        else { $parameters{output1} = $pwd . "/" . $fileName . ".collapsed" . $fileSuffix; }
        $parameters{isPairEnds} = 0;
    }
    elsif ( defined $opt_1 && defined $opt_2 )  {
        $parameters{input1} = $opt_1;
        my ($fileName1, $fileDir1, $fileSuffix1) = fileparse($parameters{input1}, qr/\.[^.]*/);
        $parameters{input2} = $opt_2;
        my ($fileName2, $fileDir2, $fileSuffix2) = fileparse($parameters{input2}, qr/\.[^.]*/);

        if ( defined $opt_p ) { $parameters{output1} = $opt_p; }
        else { $parameters{output1} = $pwd . "/" . $fileName1 . ".collapsed" . $fileSuffix1; }
        if ( defined $opt_q ) { $parameters{output2} = $opt_q; }
        else { $parameters{output2} = $pwd . "/" . $fileName2 . ".collapsed" . $fileSuffix2; }

        die "Error! Paired ends reads should have different names.\n" if ( $fileName1 eq $fileName2 );
        $parameters{isPairEnds} = 1;
    }
    else { die $usage; }

    if ( defined $opt_a ) { $parameters{isFasta} = 1; }
    if ( defined $opt_f ) { $parameters{FASTA} = $opt_f; }

    return ( %parameters );
}

sub _mergePairEndReads {
    my ( $readFile1, $readFile2, $peFile ) = @_;

    ## should test whether they are of the same length
    print STDERR "merge two PE fastq files...\n\t", `date`;
    system ( "paste $readFile1 $readFile2 > $peFile" ); 
    
    1;
}

sub _splitPairEndReads {
    my ( $peFile, $readFile1, $readFile2 ) = @_;

    print STDERR "split into two PE fastq files...\n\t", `date`;
    open ( PE, $peFile ); open ( R1, ">$readFile1" ); open ( R2, ">$readFile2" );
    while ( my $line = <PE> ) {
        my ( $r1, $r2 ) = split ( /\t/, $line );
        print R1 $r1, "\n"; print R2 $r2;
    }
    close PE; close R1; close R2; 
    1;
}

sub _estimateSplit
{
    my $file = shift;
    my $fileSize = -s $file;

    if ( not -e $file ) { die "File $file does not exists!\n"; }
    my $pos = 1;
    my $len = int ( _log4 ( $fileSize / 1000000000 ) );
    open my $fileHandler , "<", $file;
    my $line = <$fileHandler>; $line =<$fileHandler>;
    close $fileHandler;
    my $readLen = length ( $line );
    if ( $readLen <= 10 ) { $len = 0; }
    else { 
        $pos = int ( $readLen / 4 + 10 ); 
        $pos = 1 if ( $pos < 1 );
        while ( ( $pos + $len ) > $readLen ) {
            $pos = $pos - $len;
            if ( $pos < 1 ) {
                $pos = 1;
                last;
            }
        }
    }

    return ( $pos, $len );
}

sub _getSplitedFile
{
    my $splitResults = shift;

    my @data = split ( /\s/, $splitResults );
    my @filesSplited = @data[5..$#data];

    return @filesSplited;
}

sub _parseCollapseOutput
{
    my $collapseResults = shift;

    chomp $collapseResults;
    my ( $volReads ) = ( $collapseResults =~ /Total count: (\d+),/ );
    my ( $volUniqReads ) = ( $collapseResults =~ /unique count: (\d+),/ );
    my ( $uniqRatio ) = ( $collapseResults =~ /unique ratio: ([\.\d]+)./ );
    my ( $fasta ) = ( $collapseResults =~ /fasta file: (.+)$/ );

    return ( $volReads, $volUniqReads, $uniqRatio, $fasta );
}


sub _log4
{
    my $num = shift;
    return log ($num) / log (4);
}
