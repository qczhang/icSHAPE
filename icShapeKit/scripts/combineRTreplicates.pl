#! /usr/bin/perl
#
#
use strict;
use warnings;
use File::Basename;
use Getopt::Std;

use vars qw ($opt_h $opt_V $opt_D $opt_i $opt_o $opt_f  );
&getopts('hVDi:o:f:');

my $usage = <<_EOH_;
## --------------------------------------
Calculate average signals

Command:
 $0 -i input_signal_files -o output_combined_signal_file 
# what it is:
 -i     input signal files
 -o     output combined signal files

# more options:
 -f     input format

_EOH_
;


&main();

sub main 
{
    my %parameters = &init();

    my %trans_len = ();
    my %trans_rpkm = ();
    my %trans_scalingFactor = ();
    my %trans_rtstop = ();
    my %trans_bd = ();
    my %trans_enrich = ();

    my @inputFiles = split ( /:/, $parameters{input} );
    foreach my $inputFile ( @inputFiles ) {
        combineSignals ( $inputFile, \%trans_len, \%trans_rpkm, \%trans_scalingFactor, \%trans_rtstop, \%trans_bd, \%trans_enrich, 
            inputFormat => $parameters{inputFormat} );
    }

    outputSignals ( $parameters{output}, \%trans_len, \%trans_rpkm, \%trans_scalingFactor, \%trans_rtstop, \%trans_bd, \%trans_enrich, 
        inputFormat => $parameters{inputFormat} );

    1;
}


sub init 
{
    my %parameters = ();

    die $usage if ( $opt_h || ( not defined $opt_i ) || ( not defined $opt_o ) );
    $opt_V = 0 if ( not defined $opt_V );
    $opt_D = 0 if ( not defined $opt_D );

    $parameters{input} = $opt_i;
    $parameters{output} = $opt_o;

    if ( defined $opt_f )  { $parameters{inputFormat} = $opt_f; }
    else { $parameters{inputFormat} = "bd";  }

    return ( %parameters );
}

sub combineSignals 
{
    my $signalFile = shift;
    my $ref_trans_len = shift;
    my $ref_trans_rpkm = shift;
    my $ref_trans_scalingFactor = shift;
    my $ref_trans_rtstop = shift;
    my $ref_trans_bd = shift;
    my $ref_trans_enrich = shift;
    my %parameters = @_;

    print STDERR "read signal from $signalFile\n", `date`;
    my ( $id, $idOld, $type, $len, $rpkm, $scalingFactor, $score, @rtstops, @baseDensities, @data ) = ( "", "", "", 0, 0, 0, 0, (), () );
    my $line = "";  my $lineOld = ""; my $lineCount = 0;
    open ( SG, $signalFile );
    while ( $line = <SG> ) {
        next if ( $line =~ /^#/ );
        $lineCount++;
        print STDERR "  line: $lineCount\n" if ( $lineCount % 10000 == 0 );

        chomp $line;
        ( $id, @data ) = split ( /\t/, $line );

        if ( $id eq $idOld ) {
            if ( $parameters{inputFormat} =~ /normalized/ ) {
                if ( $line =~ /baseDensity/ ) {
                    if ( $lineOld =~ /RTstop/ ) {
                        ( $id, $len, $type, $rpkm, $scalingFactor, @baseDensities ) = split ( /\t/, $line );
                        ( $id, $len, $type, $rpkm, $scalingFactor, @rtstops ) = split ( /\t/, $lineOld );
                    }
                }
                elsif ( $line =~ /RTstop/ ) {
                    if ( $lineOld =~ /baseDensity/ ) {
                        ( $id, $len, $type, $rpkm, $scalingFactor, @baseDensities ) = split ( /\t/, $lineOld );
                        ( $id, $len, $type, $rpkm, $scalingFactor, @rtstops ) = split ( /\t/, $line );
                    }
                }
            }
            else {
                ( $id, $len, $rpkm, @baseDensities ) = split ( /\t/, $lineOld );
                ( $id, $len, $rpkm, @rtstops ) = split ( /\t/, $line );
            }

            if ( not defined $ref_trans_len->{$id} ) {
                $ref_trans_len->{$id} = $len;
                $ref_trans_rpkm->{$id} = $rpkm;
                $ref_trans_scalingFactor->{$id} = $scalingFactor if ( $parameters{inputFormat} =~ /normalized/ );
                $ref_trans_rtstop->{$id} = [ @rtstops ];
                $ref_trans_bd->{$id} = [ @baseDensities ];
            }
            else {
                die "Inconsistent transcript $id!\n" if ( $len != $ref_trans_len->{$id} );
                $ref_trans_rpkm->{$id} .= "," . $rpkm;
                $ref_trans_scalingFactor->{$id} .= "," . $scalingFactor if ( $parameters{inputFormat} =~ /normalized/ );
                if ( $parameters{inputFormat} =~ /normalized/ ) {
                    $ref_trans_scalingFactor->{$id} .= "," . $scalingFactor;
                    for ( my $idxPos = 0; $idxPos < $len; $idxPos++ ) {
                        $ref_trans_rtstop->{$id}[$idxPos] += $rtstops[$idxPos];
                        $ref_trans_bd->{$id}[$idxPos] += $baseDensities[$idxPos];
                    }
                }
                else {
                    for ( my $idxPos = 0; $idxPos <= $len; $idxPos++ ) {
                        $ref_trans_rtstop->{$id}[$idxPos] += $rtstops[$idxPos];
                        $ref_trans_bd->{$id}[$idxPos] += $baseDensities[$idxPos];
                    }
                }
            }
        }
        else {
            $idOld = $id;
            $lineOld = $line;
        }
    }
    close SG;

    1;
}

sub outputSignals 
{
    my $signalFile = shift;
    my $ref_trans_len = shift;
    my $ref_trans_rpkm = shift;
    my $ref_trans_scalingFactor = shift;
    my $ref_trans_rtstop = shift;
    my $ref_trans_bd = shift;
    my $ref_trans_enrich = shift;
    my %parameters = @_;

    print STDERR "output signal to $signalFile\n", `date`;
    open ( SG, ">$signalFile" );
    foreach my $id ( keys %{$ref_trans_len} ) {
        if ( $parameters{inputFormat} =~ /normalized/ ) {
            print SG $id, "\t", $ref_trans_len->{$id}, "\t", $ref_trans_rpkm->{$id}, "\t", $ref_trans_scalingFactor->{$id};
            for ( my $idxPos = 0; $idxPos < $ref_trans_len->{$id}; $idxPos++ ) { print SG "\t", $ref_trans_bd->{$id}[$idxPos]; }
            print SG "\n";
            print SG $id, "\t", $ref_trans_len->{$id}, "\t", $ref_trans_rpkm->{$id}, "\t", $ref_trans_scalingFactor->{$id};
            for ( my $idxPos = 0; $idxPos < $ref_trans_len->{$id}; $idxPos++ ) { print SG "\t", $ref_trans_rtstop->{$id}[$idxPos]; }
            print SG "\n";
        }
        else {
            print SG $id, "\t", $ref_trans_len->{$id}, "\t", $ref_trans_rpkm->{$id};
            for ( my $idxPos = 0; $idxPos <= $ref_trans_len->{$id}; $idxPos++ ) { print SG "\t", $ref_trans_bd->{$id}[$idxPos]; }
            print SG "\n";
            print SG $id, "\t", $ref_trans_len->{$id}, "\t", $ref_trans_rpkm->{$id};
            for ( my $idxPos = 0; $idxPos <= $ref_trans_len->{$id}; $idxPos++ ) { print SG "\t", $ref_trans_rtstop->{$id}[$idxPos]; }
            print SG "\n";
        }
    }
    close SG;

    1;
}

