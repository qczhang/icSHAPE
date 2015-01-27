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
my $simpleTrim = "$FindBin::Bin/../bin/simpleTrim";

use vars qw ($opt_h $opt_V $opt_D $opt_U $opt_1 $opt_2 $opt_o $opt_p $opt_q $opt_a $opt_l $opt_t $opt_c $opt_m );
&getopts('hVDU:1:2:o:p:q:a:l:t:c:m:');

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

 -c     fastq coding
 -a     adapter sequences
 -l     crop of leading nucleotides
 -t     crop of tailing nucleotides
 -m     min length

_EOH_
;

&main();

sub main {
    my %parameters = &init();

    my $inFile1 = $parameters{input1};
    if ( not $parameters{isPairEnds} )  {
        print STDERR "Trimming file SE reads $inFile1...\n\t", `date`;
        if ( -e $parameters{output1} ) 
            { print STDERR "Warning! $parameters{output1} exisits, will be overwritten.\n"; print STDERR `/bin/rm $parameters{output1}`; }

        if ( ( defined $parameters{leading} ) and ( $parameters{leading} > 0 ) ) {
            print STDERR "$simpleTrim -U $inFile1 -o $parameters{output1}.tmp -d $parameters{leading} -m $parameters{minLength}\n";
            print STDERR `$simpleTrim -U $inFile1 -o $parameters{output1}.tmp -d $parameters{leading} -m $parameters{minLength}`;
            if ( $? ) { die "Error in removing leading nucleotides for file $inFile1!\n"; }

            $inFile1 = $parameters{output1} . ".tmp";
        }

        if ( defined $parameters{adapter} ) {
            print STDERR "java -mx512m -jar $trimmomatic SE -$parameters{coding} -trimlog $parameters{trimlog} $inFile1 $parameters{output1} ILLUMINACLIP:$parameters{adapter}:2:30:10 TRAILING:20 MINLEN:$parameters{minLength}\n";
            print STDERR `java -mx512m -jar $trimmomatic SE -$parameters{coding} -trimlog $parameters{trimlog} $inFile1 $parameters{output1} ILLUMINACLIP:$parameters{adapter}:2:30:10 TRAILING:20 MINLEN:$parameters{minLength}`;
            if ( $? ) { die "Error in running trimmomatic for file $inFile1!\n"; }
        }
    }
    else {
        my $inFile2 = $parameters{input2};
        print STDERR "Trimming file PE reads $parameters{input1} and $parameters{input2}...\n\t", `date`;
        if ( -e $parameters{output1} ) 
        { print STDERR "Warning! $parameters{output1} exisits, will be overwritten.\n"; print STDERR `/bin/rm $parameters{output1}`; }
        if ( -e $parameters{output2} ) 
            { print STDERR "Warning! $parameters{output2} exisits, will be overwritten.\n"; print STDERR `/bin/rm $parameters{output2}`; }

        if ( ( defined $parameters{leading} ) and ( $parameters{leading} > 0 ) ) {
            print STDERR "$simpleTrim -1 $parameters{input1} -2 $parameters{input2} -p $parameters{output1}.tmp -q $parameters{output2}.tmp -d $parameters{leading} -m $parameters{minLength}\n";
            print STDERR `$simpleTrim -1 $parameters{input1} -2 $parameters{input2} -p $parameters{output1}.tmp -q $parameters{output2}.tmp -d $parameters{leading} -m $parameters{minLength}`;
            if ( $? ) { die "Error in removing leading nucleotides for file $inFile1 and $inFile2!\n"; }

            $inFile1 = $parameters{output1} . ".tmp";
            $inFile2 = $parameters{output2} . ".tmp";
        }

        if ( defined $parameters{adapter} ) {
            print STDERR "java -mx512m -jar $trimmomatic PE -$parameters{coding} -trimlog  $parameters{trimlog} $inFile1 $inFile2 $parameters{output1} $parameters{output1}.unpaired $parameters{output2} $parameters{output2}.unpaired ILLUMINACLIP:$parameters{adapter}:2:30:10 TRAILING:20 MINLEN:$parameters{minLength}\n";
            print STDERR `java -mx512m -jar $trimmomatic PE -$parameters{coding} -trimlog  $parameters{trimlog} $inFile1 $inFile2 $parameters{output1} $parameters{output1}.unpaired $parameters{output2} $parameters{output2}.unpaired ILLUMINACLIP:$parameters{adapter}:2:30:10 TRAILING:20 MINLEN:$parameters{minLength}`;
            if ( $? ) { die "Error in running trimmomatic for file $inFile1 and $inFile2!\n"; }
        }
    }

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
        $parameters{trimlog} = $parameters{output1}. ".trimlog";
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
        $parameters{trimlog} = $parameters{output1}. ".trimlog";
        $parameters{isPairEnds} = 1;
    }
    else { die $usage; }

    if ( defined $opt_c ) { $parameters{coding} = $opt_c; }
    else { $parameters{coding} = "phred33"; }
    if ( defined $opt_a ) { $parameters{adapter} = $opt_a; }
    if ( defined $opt_l ) { $parameters{leading} = $opt_l; }
    if ( defined $opt_t ) { $parameters{tailing} = $opt_t; }
    if ( defined $opt_m ) { $parameters{minLength} = $opt_m; }

    return ( %parameters );
}
