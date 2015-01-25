#! /usr/bin/perl
#
# trimming
#   history 20150101
#   version 0.1
#   copyright @cliff zhang
#
use strict;
use warnings;
use File::Basename;
use Getopt::Std;
use FindBin;

my $trimmomatic = "$FindBin::Bin/../bin/trimmomatic-0.30.jar";
my $simpletrim = "$FindBin::Bin/../bin/simpleTrim";

use vars qw ($opt_h $opt_V $opt_D $opt_U $opt_1 $opt_2 $opt_o $opt_p $opt_q $opt_a $opt_l $opt_t $opt_c );
&getopts('hVDU:1:2:o:p:q:a:l:t:c:');

my $usage = <<_EOH_;
## --------------------------------------
trimming fastq file to remove possible adapter contamination and also barcode regions

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

 -a     adapter sequences

 -l     crop of leading nucleotides
 -t     crop of tailing nucleotides
 -c     fastq coding

_EOH_
;

&main();

sub main {
    my %parameters = &init();

    my $inFile = $parameters{input1}; my $outFile = $parameters{output1};
    my $cmd = "$parameters{java} -jar $trimmomatic ";
    if ( not $parameters{isPairEnds} )  {
        print STDERR "Trimming file SE reads $inFile...\n\t", `date`;
        if ( -e $parameters{output1} ) { print STDERR "Warning! $parameters{output1} exisits, will be overwritten.\n"; print STDERR `/bin/rm $parameters{output1}`; }
        $cmd .= "SE -$parameters{coding} $trimlog $parameters{input1} $parameters{output1} $parameters{output1}.unpaired $options";
    }
    else {
        print STDERR "Trimming file PE reads $parameters{input1} and $parameters{input2}...\n\t", `date`;
        if ( -e $parameters{output1} ) { print STDERR "Warning! $parameters{output1} exisits, will be overwritten.\n"; print STDERR `/bin/rm $parameters{output1}`; }
        if ( -e $parameters{output2} ) { print STDERR "Warning! $parameters{output2} exisits, will be overwritten.\n"; print STDERR `/bin/rm $parameters{output2}`; }
        $inFile = "/tmp/tmp" . $$ . ".fastq"; $outFile = "/tmp/tmpOut" . $$ . ".fastq";
        _mergePairEndReads ( $parameters{input1}, $parameters{input2}, $inFile );

        $cmd = "PE -$parameters{coding} $trimlog $parameters{input1} $parameters{input2} $parameters{output1} $parameters{output1}.unpaired $parameters{output2}.unpaired $options";
    }

    system ( "$cmd" );
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
        else { $parameters{output1} = $pwd . "/" . $fileName . ".trimmed" . $fileSuffix; }
        $parameters{isPairEnds} = 0;
    }
    elsif ( defined $opt_1 && defined $opt_2 )  {
        $parameters{input1} = $opt_1;
        my ($fileName1, $fileDir1, $fileSuffix1) = fileparse($parameters{input1}, qr/\.[^.]*/);
        $parameters{input2} = $opt_2;
        my ($fileName2, $fileDir2, $fileSuffix2) = fileparse($parameters{input2}, qr/\.[^.]*/);

        if ( defined $opt_p ) { $parameters{output1} = $opt_p; }
        else { $parameters{output1} = $pwd . "/" . $fileName1 . ".trimmed" . $fileSuffix1; }
        if ( defined $opt_q ) { $parameters{output2} = $opt_q; }
        else { $parameters{output2} = $pwd . "/" . $fileName2 . ".trimmed" . $fileSuffix2; }

        die "Error! Paired ends reads should have different names.\n" if ( $fileName1 eq $fileName2 );
        $parameters{isPairEnds} = 1;
    }
    else { die $usage; }

    if ( defined $opt_a ) { $parameters{adapter} = $opt_a; }
    if ( defined $opt_l ) { $parameters{leading} = $opt_l; }
    if ( defined $opt_t ) { $parameters{tailing} = $opt_t; }
    if ( defined $opt_c ) { $parameters{coding} = $opt_c; }

    return ( %parameters );
}
