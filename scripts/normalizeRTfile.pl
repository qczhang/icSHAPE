#! /usr/bin/perl
#
#
use strict;
use warnings;
use File::Basename;
use Getopt::Std;

use vars qw ($opt_h $opt_V $opt_D $opt_i $opt_o $opt_m $opt_d $opt_l $opt_f );
&getopts('hVDi:o:r:m:d:l:f:');

my $usage = <<_EOH_;
## --------------------------------------
Normalize base density file

Command:
$0 -i input_baseDensity_file -o output_normalized_baseDensity_file

# what it is:
 -i     input_baseDensity_file
 -o     output_normalized_baseDensity_file

# more options:
 -d     head to skip
 -l     tail to skip
 -f     scaling form (the benchmarked value will be scaled to this)

 -m     normalize method

_EOH_
;

&main();

sub main {
    my %parameters = &init();
    &normalizeBaseDensity ( $parameters{input}, $parameters{output}, $parameters{normalizeMethod}, $parameters{headToSkip}, $parameters{tailToSkip}, $parameters{scalingForm} );
}

sub init {
    my %parameters = ();

    die $usage if ( $opt_h || ( not defined $opt_i ) || ( not defined $opt_o ) );
    $opt_V = 0 if ( not defined $opt_V );
    $opt_D = 0 if ( not defined $opt_D );

    $parameters{input} = $opt_i;
    $parameters{output} = $opt_o;

    if ( defined $opt_m ) { $parameters{normalizeMethod} = $opt_m; }
    else { $parameters{normalizeMethod} = "median:decile2"; }
    if ( defined $opt_d ) { $parameters{headToSkip} = $opt_d; }
    else { $parameters{headToSkip} = 0; }
    if ( defined $opt_l ) { $parameters{tailToSkip} = $opt_l; }
    else { $parameters{tailToSkip} = 0; }
    if ( defined $opt_f ) { $parameters{scalingForm} = $opt_f; }
    else { $parameters{scalingForm} = 100; }

    return ( %parameters );
}


sub normalizeBaseDensity {
    my ( $baseDensityFile, $outputFile, $normalizeMethod, $headToSkip, $tailToSkip, $scalingForm ) = @_;
    print STDERR "Normalize base density from file $baseDensityFile...\n\t", `date`;

    $headToSkip++; 
    open ( IN, $baseDensityFile );
    open ( OUT, ">$outputFile" ) or die "Error opening $outputFile for writing transcript base density.\n";
    print OUT "#transcript\tlength\ttype\tbase_frequency, start from position 1.\n";

    my $transcript = "";  my $len = 0;  my $rpkm = 0; my @baseDensity = ();
    while ( my $line = <IN> )  {
        next if ( $line =~ /^#/ );

        my $scalingFactor = 1;
        chomp $line;
        ## remember that the first element in @baseDensity is base zero
        #  in RT stop file, it means no modified base in the read
        #  in background file, it means null, so should always be zero
        ( $transcript, $len, $rpkm, @baseDensity ) = split ( /\t/, $line );
        my $trimed_last = $len - $tailToSkip;      ## again, to avoid base zero
        while ( $trimed_last < ( $headToSkip + 40 ) ) {
            print STDERR "Warning! Transcript $transcript too short.\n";

            $headToSkip = int ( $headToSkip/2 );
            $tailToSkip = int ( $tailToSkip/2 );
            $trimed_last = $len - $tailToSkip;
            last if ( ( $headToSkip == 0 ) && ( $tailToSkip == 0 ) );
        }

        $scalingFactor = &calcScalingFactor ( \@baseDensity, $headToSkip, $trimed_last, $normalizeMethod );

        if ( $scalingFactor > 1 ) {
            $scalingFactor = ( $scalingFactor / $scalingForm );

            print OUT $transcript, "\t", $len, "\tbaseDensity\t", $rpkm, "\t", $scalingFactor;
            if ( $normalizeMethod =~ /log/i ) {
                for ( my $idx = 1; $idx <= $len; $idx++ )  {
                    print OUT "\t", sprintf ( "%.3f", &log2 ( $baseDensity[$idx]/$scalingFactor + 1 ) );
                }
            }
            else {
                for ( my $idx = 1; $idx <= $len; $idx++ )  {
                    print OUT "\t", sprintf ( "%.3f", ( $baseDensity[$idx]/$scalingFactor ) );
                }
            }
            print OUT "\n";
        }

        $line = <IN>;
        chomp $line;
        ( $transcript, $len, $rpkm, @baseDensity ) = split ( /\t/, $line );
        $scalingFactor = &calcScalingFactor ( \@baseDensity, $headToSkip, $trimed_last, $normalizeMethod );

        if ( $scalingFactor > 1 ) {
            $scalingFactor = ( $scalingFactor / $scalingForm );

            print OUT $transcript, "\t", $len, "\tRTstop\t", $rpkm, "\t", $scalingFactor;
            if ( $normalizeMethod =~ /log/i ) {
                for ( my $idx = 1; $idx <= $len; $idx++ )  {
                    print OUT "\t", sprintf ( "%.3f", &log2 ( $baseDensity[$idx]/$scalingFactor + 1 ) );
                }
            }
            else {
                for ( my $idx = 1; $idx <= $len; $idx++ )  {
                    print OUT "\t", sprintf ( "%.3f", ( $baseDensity[$idx]/$scalingFactor ) );
                }
            }
            print OUT "\n";
        }
    }
    close IN;
    close OUT;

    1;
}

sub calcScalingFactor  {
    my ( $ref_array, $start, $end, $normalizeMethod ) = @_;

    my @sortedTrimed_array = sort { $b <=> $a } ( @{$ref_array}[$start..$end] );

    my $len = scalar(@sortedTrimed_array); 
    my @rangeOfSelection = ();
    if ( $normalizeMethod =~ /smart/i ) {
        ## smartly skip those non-zero elements
        for ( my $idx = 0; $idx < $len; $idx++ )  {
            last if ( $sortedTrimed_array[$idx] <= 0 );
            push @rangeOfSelection, $sortedTrimed_array[$idx];
        }
    }
    else {
        my $selectStart = 0; 
        my $selectEnd = $len - 1;
        if ( $normalizeMethod =~ /upper/i ) {
            ## get the upper half
            $selectEnd = int ( $len / 2 ) - 1;
        }
        elsif ( $normalizeMethod =~ /quartile(\d+)/i ) {
            ## get a specified quarter
            $selectStart = int ( ($1-1)*$len / 4 );
            $selectEnd = int ( $1*$len / 4 ) - 1;
        }
        elsif ( $normalizeMethod =~ /decile(\d+)/i ) {
            ## get a specified decile
            $selectStart = int ( ($1-1)*$len / 10 );
            $selectEnd = int ( $1*$len / 10 ) - 1;
        }
        elsif ( $normalizeMethod =~ /vigintile(\d+)/i ) {
            ## get a specified decile
            $selectStart = int ( ($1-1)*$len / 20 );
            $selectEnd = int ( $1*$len / 20 ) - 1;
        }

        @rangeOfSelection = @sortedTrimed_array[$selectStart..$selectEnd];
    }
    $len = scalar(@rangeOfSelection); 

    my $scalingFactor = 1;
    if ( $normalizeMethod =~ /median/i ) {
        my $median = 1;
        if ( $len % 2 == 0 ) {  $median = ($rangeOfSelection[$len/2-1] + $rangeOfSelection[$len/2]) /2;  }
        else {  $median = $rangeOfSelection[($len-1)/2];  }
        $scalingFactor = $median;
    }
    elsif ( $normalizeMethod =~ /mean/i ) {
        my $mean = eval ( join "+", @rangeOfSelection );
        $mean /= scalar (@rangeOfSelection);
        $scalingFactor = $mean;
    }
    elsif ( $normalizeMethod =~ /peak/i ) {
        $scalingFactor = $rangeOfSelection[0];
    }

    return $scalingFactor;
}

sub log2 {
    my $num = shift;
    return log($num)/log(2);
}


