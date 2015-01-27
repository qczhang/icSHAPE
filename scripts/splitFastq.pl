#! /usr/bin/perl
#
# split fastq files into volumns or libraries
# 
#
use strict;
use warnings;
use File::Basename;
use Getopt::Std;
use FindBin;

my $splitFastqbin = "$FindBin::Bin/../bin/splitByBarcode";

use vars qw ($opt_h $opt_V $opt_D $opt_U $opt_1 $opt_2 $opt_l $opt_b $opt_s $opt_d );
&getopts('hVD1:2:U:l:b:s:d:');

my $usage = <<_EOH_;
## --------------------------------------
split fastq file

Command:
$0 -1 fastq_PE_reads_1 -2 fastq_PE_reads_2 -U fastq_SE_reads -l split_by_library -b barcode_position_and_length

# what it is:
 -1     paired ends read 1
 -2     paired ends read 2
 -U     single ends read
 -b     barcode position and length, separated by semicolon. e.g. 5:4
 -l     specify a string to denote different libraries. Use the format BARCODE1:LIB_NAME1::BARCODE2:LIB_NAME2::...

# more options:
 -s     simple statistics if split by library
 -d     output directory

_EOH_
;

&main();

sub main {
    my %parameters = &init();

    my $outDir = $parameters{outDirectory};
    if ( not -e $outDir ) { mkdir $outDir or die "Error! Create $outDir failed.\n"; }

    my @files1 = (); my @files2 = ();
    my @fs1 = split ( /:/, $parameters{input1} );
    foreach my $f ( @fs1 ) { foreach ( glob ( $f ) ) { push @files1, $_; } }
    if ( $parameters{isPairEnds} ) {
        my @fs2 = split ( /:/, $parameters{input2} );
        foreach my $f ( @fs2 ) { foreach ( glob ( $f ) ) { push @files2, $_; } }
        die "Error! input different number of files for PE mode\n" if ( scalar ( @files1 ) != scalar ( @files2 ) );
    }

    for ( my $fileIdx = 0; $fileIdx < scalar ( @files1 ); $fileIdx++ ) {
        my $file1 = $files1[$fileIdx];
        if ( not $parameters{isPairEnds} )  {
            # note that append will lump reads into files in the outputdirectory with names like LIB_DMSO1.fastq, LIB_DMSO2.fastq...
            print STDERR "$splitFastqbin $file1 $parameters{outDirectory} $parameters{BCPOS} $parameters{BCLENGTH} append $parameters{libraryCode}\n\t", `date`;
            print STDERR `$splitFastqbin $file1 $parameters{outDirectory} $parameters{BCPOS} $parameters{BCLENGTH} append $parameters{libraryCode}`;
            if ( $? != 0 ) { die "Error! splitting file $file1 not successful!\n"; } ## check return status
        }
        else {
            my $file2 = $files2[$fileIdx];
            my $inFastqFile = "/tmp/tmp$$.fastq";
            _mergePairEndReads ( $file1, $file2, $inFastqFile );
            print STDERR "$splitFastqbin $inFastqFile $parameters{outDirectory} $parameters{BCPOS} $parameters{BCLENGTH} append $parameters{libraryCode}\n\t", `date`;
            print STDERR `$splitFastqbin $inFastqFile $parameters{outDirectory} $parameters{BCPOS} $parameters{BCLENGTH} append $parameters{libraryCode}`;
            if ( $? != 0 ) { print STDERR `/bin/rm -f $inFastqFile`; die "Error! splitting file $inFastqFile not successful!\n"; }
            print STDERR `/bin/rm -f $inFastqFile`;
        }
    }

    if ( $parameters{isPairEnds} )  {
        foreach my $lib ( keys %{$parameters{library}} ) {
            my $outFile = $parameters{outDirectory} . "/" . $lib . ".fastq"; ##print $outFile, "\n";
            my $outFile1 = $parameters{outDirectory} . "/r1." . $lib . ".fastq";
            my $outFile2 = $parameters{outDirectory} . "/r2." . $lib . ".fastq";
            _splitPairEndReads ( $outFile, $outFile1, $outFile2 );
            print STDERR `/bin/rm -f $outFile`;
        }
    }

    my $statFile = $parameters{outDirectory} . "/splitFastq.stat";
    if ( defined $parameters{BCCOUNT} ) {
        print STDERR `/bin/mv $statFile $parameters{BCCOUNT}`;
        $statFile = $parameters{BCCOUNT};
    }
    print STDERR "Spliting successful! Check spliting statistics at $statFile\n";

    1;
}


sub init 
{
    my %parameters = ();

    die $usage if ( $opt_h );
    $opt_V = 0 if ( not defined $opt_V );
    $opt_D = 0 if ( not defined $opt_D );

    if ( defined $opt_U ) {
        my ($inFileName, $inFileDir, $inFileSuffix) = fileparse($opt_U, qr/\.[^.]*/);
        die "Error! input file must be in fastq format!\n" if ( $inFileSuffix ne ".fastq" );
        $parameters{input1} = $opt_U;
        $parameters{isPairEnds} = 0;
    }
    elsif ( defined $opt_1 && defined $opt_2 )  {
        my ($inFileName, $inFileDir, $inFileSuffix) = fileparse($opt_1, qr/\.[^.]*/);
        die "Error! input file must be in fastq format!\n" if ( $inFileSuffix ne ".fastq" );
        ($inFileName, $inFileDir, $inFileSuffix) = fileparse($opt_2, qr/\.[^.]*/);
        die "Error! input file must be in fastq format!\n" if ( $inFileSuffix ne ".fastq" );

        $parameters{input1} = $opt_1;
        $parameters{input2} = $opt_2;
        $parameters{isPairEnds} = 1;
    }
    else { die $usage; }

    if ( defined $opt_d ) { $parameters{outDirectory} = $opt_d; }
    else { $parameters{outDirectory} = `pwd`; chomp $parameters{outDirectory}; }

    if ( defined $opt_l ) {
        die "Error! library barcode position and length must be specified in split_by_library mode!\n" if ( not defined $opt_b ); 
        $parameters{libraryCode} = $opt_l;

        my %library = ();
        my @lib_code = split ( /::/, $opt_l );
        foreach my $lc ( @lib_code ) {
            my ( $code, $lib ) = split ( /:/, $lc );
            $library{$lib} = $code;

            my $outFile = $parameters{outDirectory}. "/" . $lib . ".fastq";
            if ( -e $outFile ) {
                print STDERR "Warning! $outFile exists...will be erased.\n";
                print STDERR `/bin/rm $outFile`;
            }
        }
        $parameters{library} = \%library;

        my $unmatchedFile = $parameters{outDirectory} . "/unmatched.fastq";
        if ( -e $unmatchedFile ) {
            print STDERR "Warning! $unmatchedFile exists...will be erased.\n";
            print STDERR `/bin/rm $unmatchedFile`;
        }

        my ( $barcodePos, $barcodeLen ) = split ( /:/, $opt_b );
        $parameters{BCPOS} = $barcodePos;
        $parameters{BCLENGTH} = $barcodeLen;
    }

    if ( defined $opt_s )  { $parameters{BCCOUNT} = $opt_s; }

    return ( %parameters );

    1;
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

    ## should test whether they are of the same length
    print STDERR "split into two PE fastq files...\n\t", `date`;
    open ( PE, $peFile );
    open ( R1, ">$readFile1" );
    open ( R2, ">$readFile2" );
    while ( my $line = <PE> ) {
        my ( $r1, $r2 ) = split ( /\t/, $line );
        print R1 $r1, "\n";
        print R2 $r2;
    }
    close PE;
    close R1;
    close R2;
    
    1;
}

