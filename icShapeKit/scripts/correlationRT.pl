#! /usr/bin/perl
#
#
use strict;
use warnings;
use File::Basename;
use Getopt::Std;
use Data::Dumper;

use lib qw "/home/qczhang/lib/perllib";
use Statistics::Basic qw(:all);
use Locale::Schedule::Simple qw( &waitForFile &finishTag );

use vars qw ($opt_h $opt_V $opt_D $opt_1 $opt_2 $opt_s $opt_f $opt_c $opt_w $opt_t $opt_b $opt_T );
&getopts('hVD1:2:o:s:f:c:w:t:b:T:');

my $usage = <<_EOH_;
## --------------------------------------
Calculate average signals

Command:
 $0 -1 input_signal_file1 -2 input_signal_file2 
# what it is:
 -1     input signal file 1
 -2     input signal file 2

# more options:
 -f     input format
 -s     selection strategy
 -w     smooth window size
 -c     correlation window size
 -t     threshold to select window for correlation calculation
 -T     threshold of coverage

 -b     background base density

_EOH_
;


&main();

sub main 
{
    my %parameters = &init();

    my ( $ref_trans_len1, $ref_trans_rpkm1, $ref_trans_rtstop1 ) = &readSignals ( $parameters{input1}, 
        inputFormat => $parameters{inputFormat},
        coverageThreshold => $parameters{coverageThreshold});
    my ( $ref_trans_len2, $ref_trans_rpkm2, $ref_trans_rtstop2 ) = &readSignals ( $parameters{input2}, 
        inputFormat => $parameters{inputFormat},
        coverageThreshold => $parameters{coverageThreshold});

    my $ref_trans_background = undef;
    if ( ( defined $parameters{selectionStrategy} ) and ( $parameters{selectionStrategy} eq "useBackground" ) ) {
        ( $ref_trans_len1, $ref_trans_rpkm1, $ref_trans_background ) 
            = &readSignals ( $parameters{background}, readFormat => $parameters{selectionStrategy}, inputFormat => "bd" );
    }

    foreach my $trans ( keys %{$ref_trans_rtstop1} ) {
        next if ( not defined $ref_trans_rtstop2->{$trans} );
        next if ( ( not defined $ref_trans_len1->{$trans} ) or ( not defined $ref_trans_len2->{$trans} ) );
        if ( $ref_trans_len1->{$trans} != $ref_trans_len2->{$trans} ) { print "Error! skip transcript of different length in two replicates.\n"; next; }

        my $rpkm = $ref_trans_rpkm1->{$trans};
        if ( ( not defined $parameters{selectionStrategy} ) or ( $parameters{selectionStrategy} ne "useBackground" ) )  
            { $rpkm = ( $ref_trans_rpkm1->{$trans} + $ref_trans_rpkm2->{$trans} ) / 2; }
        my $ref_smoothed1 = $ref_trans_rtstop1->{$trans}; 
        my $ref_smoothed2 = $ref_trans_rtstop2->{$trans};
        if ( defined $parameters{smoothWindow} ) {
            $ref_smoothed1 = smooth ( $ref_trans_len1->{$trans}, $ref_trans_rtstop1->{$trans} );
            $ref_smoothed2 = smooth ( $ref_trans_len2->{$trans}, $ref_trans_rtstop2->{$trans} );
        }

        if ( $parameters{selectionStrategy} and ( $parameters{selectionStrategy} =~ /background/i ) ) {
            next if ( not defined $ref_trans_background->{$trans} );
            &calcCorrelation ( $trans, $ref_trans_len1->{$trans}, $rpkm, $ref_smoothed1, $ref_smoothed2, 
                windowLen => $parameters{correlationWindow}, 
                winThreshold => $parameters{windowThreshold}, 
                useBackground => $ref_trans_background->{$trans} );
        }
        else {
            &calcCorrelation ( $trans, $ref_trans_len1->{$trans}, $rpkm, $ref_smoothed1, $ref_smoothed2,
                windowLen => $parameters{correlationWindow}, 
                winThreshold => $parameters{windowThreshold} );
        }
    }

    1;
}


sub init 
{
    my %parameters = ();

    die $usage if ( $opt_h || ( not defined $opt_1 ) || ( not defined $opt_2 ) );
    $opt_V = 0 if ( not defined $opt_V );
    $opt_D = 0 if ( not defined $opt_D );

    $parameters{input1} = $opt_1;
    $parameters{input2} = $opt_2;

    if ( defined $opt_f )  { $parameters{inputFormat} = $opt_f; }
    else { $parameters{inputFormat} = "bd";  }

    if ( defined $opt_s )  { $parameters{selectionStrategy} = $opt_s; }
    if ( defined $opt_w )  { $parameters{smoothWindow} = $opt_w; }
    if ( defined $opt_c )  { $parameters{correlationWindow} = $opt_c; } 
    else { $parameters{correlationWindow} = "all";  }
    if ( defined $opt_t )  { $parameters{windowThreshold} = $opt_t; }
    else { $parameters{windowThreshold} = 0;  }
    if ( defined $opt_T )  { $parameters{coverageThreshold} = $opt_T; }
    else { $parameters{coverageThreshold} = 0;  }

    if ( defined $opt_b )  { $parameters{background} = $opt_b; }

    if ( $parameters{inputFormat} eq "enrich" ) {
        die "Must specify background file for enrichment score!\n" if ( ( $parameters{selectionStrategy} !~ /background/i ) or ( not defined $parameters{background} ) );
    }

    print STDERR Dumper ( \%parameters ) if ( $opt_D );
    return ( %parameters );
}

sub readSignals 
{
    my $signalFile = shift;
    my %parameters = @_;

    print STDERR "read signal from $signalFile\n", `date`;
    my %trans_len = (); my %trans_rpkm = (); my %trans_rtstop = ();
    my ( $id, $len, $rpkm, $scalingFactors, $score, $fgrtstop, $pos0, @baseDensities ) = ( "", 0, 0, 0, 0, 0, 0, () );
    my @scores = ();  my @fgrtstops = ();
    my $line = "";  my $line2 = ""; my $lineCount = 0;
    open ( SG, $signalFile );
    while ( $line = <SG> ) {
        next if ( $line =~ /^#/ );
        $lineCount++;
        print STDERR "  line: $lineCount\n" if ( $lineCount % 10000 == 0 );

        if ( $parameters{inputFormat} eq "enrich" ) {
            ( $id, $len, $rpkm, $scalingFactors, @baseDensities ) = split ( /\t/, $line );
            @scores = ();  @fgrtstops = ();
            my ( $scalingFactor ) = split ( /,/, $scalingFactors );
            for ( my $idxPos = 0; $idxPos < scalar(@baseDensities); $idxPos++ ) {
                ( $score, $fgrtstop ) = split ( /,/, $baseDensities[$idxPos] );
                #$score = 0 if ( $score eq "NULL" );
                push @scores, $score;
                push @fgrtstops, $fgrtstop;
            }
            my $coverage = eval ( join ( "+", @fgrtstops ) ) / scalar ( @fgrtstops );
            if ( $coverage > $parameters{coverageThreshold} ) {
                $trans_len{$id} = $len;
                my @rpkms = split ( /,/, $rpkm );
                $trans_rpkm{$id} = eval ( join ( "+", @rpkms ) ) / scalar ( @rpkms );
                $trans_rtstop{$id} = [ @scores ];
            }
        }
        elsif ( $parameters{inputFormat} eq "bd" ) {
            $line2 = <SG>;
            if ( $parameters{readFormat} and ( $parameters{readFormat} eq "useBackground" ) ) {
                chomp $line;
                ( $id, $len, $rpkm, $pos0, @baseDensities ) = split ( /\t/, $line );
                $trans_len{$id} = $len;
                my @rpkms = split ( /,/, $rpkm );
                $trans_rpkm{$id} = eval ( join ( "+", @rpkms ) ) / scalar ( @rpkms );
                $trans_rtstop{$id} = [ @baseDensities ];
            }
            else {
                # by default it will read the RT stop line
                chomp $line2;
                ( $id, $len, $rpkm, $pos0, @baseDensities ) = split ( /\t/, $line2 );
                my $coverage = eval ( join ( "+", @baseDensities ) ) / scalar ( @baseDensities );
                if ( $coverage > $parameters{coverageThreshold} ) {
                    $trans_len{$id} = $len;
                    $trans_rpkm{$id} = $rpkm;
                    $trans_rtstop{$id} = [ @baseDensities ];
                }
            }
        }
    }
    close SG;

    return ( \%trans_len, \%trans_rpkm, \%trans_rtstop );
}

sub calcCorrelation 
{
    my ( $trans, $len, $rpkm, $ref_rtstop1, $ref_rtstop2, %parameters ) = @_;
    my $winLen = $parameters{windowLen};
    my $winThreshold = $parameters{winThreshold};

    my @window_rep1 = (); my @window_rep2 = ();
    if ( $winLen eq "all" ) {
        for ( my $idxPos = 0; $idxPos < $len; $idxPos++ ) {
            if ( ( $ref_rtstop1->[$idxPos] ne "NULL" ) and ( $ref_rtstop2->[$idxPos] ne "NULL" ) ) {
                if ( defined $parameters{useBackground} ) {
                    if ( $parameters{useBackground}->[$idxPos] ne "NULL" ) {
                        if ( $parameters{useBackground}->[$idxPos] >= $winThreshold ) {
                            push @window_rep1, $ref_rtstop1->[$idxPos];
                            push @window_rep2, $ref_rtstop2->[$idxPos];
                        }
                    }
                }
                else {
                    if ( ( $ref_rtstop1->[$idxPos] >= $winThreshold ) and ( $ref_rtstop2->[$idxPos] >= $winThreshold ) ) {
                        push @window_rep1, $ref_rtstop1->[$idxPos];
                        push @window_rep2, $ref_rtstop2->[$idxPos];
                    }
                }
            }
        }
        $winLen = scalar ( @window_rep1 );
        my $corr = correlation ( \@window_rep1, \@window_rep2 );
        print $trans, "\t", $len, "\t", $winLen, "\t", $rpkm, "\t", $corr, "\n";
    }
    else {
        my $winTotal = $winLen * $winThreshold;
        my @bg_window_rep1 = (); my @bg_window_rep2 = ();
        my $winPos = 0; my $totalCov1 = 0; my $totalCov2 = 0;
        for ( my $idxPos = 0; $idxPos < $len; $idxPos++ ) {
            if ( $winPos == $winLen )  {
                if ( defined $parameters{useBackground} ) {
                    $totalCov1 = eval ( join "+", @bg_window_rep1 );
                    $totalCov2 = $winTotal + 1;
                }
                else {
                    $totalCov1 = eval ( join "+", @window_rep1 );
                    $totalCov2 = eval ( join "+", @window_rep2 );
                }

                if ( ( $totalCov1 >= $winTotal ) and ( $totalCov2 >= $winTotal ) ) {
                    ## this means average coverage higher than 10
                    my $corr = correlation ( \@window_rep1, \@window_rep2 );
                    print $trans, "\t", $len, "\t", $rpkm, "\t", $corr, "\n";
                }

                $winPos = 0;
            }

            if ( ( $ref_rtstop1->[$idxPos] ne "NULL" ) and ( $ref_rtstop2->[$idxPos] ne "NULL" ) ) {
                if ( defined $parameters{useBackground} ) {
                    if ( $parameters{useBackground}->[$idxPos] ne "NULL" ) {
                        $window_rep1[$winPos] = $ref_rtstop1->[$idxPos];
                        $window_rep2[$winPos] = $ref_rtstop2->[$idxPos];
                        $bg_window_rep1[$winPos] = $parameters{useBackground}->[$idxPos];

                        $winPos++;
                    }
                }
                else {
                    $window_rep1[$winPos] = $ref_rtstop1->[$idxPos];
                    $window_rep2[$winPos] = $ref_rtstop2->[$idxPos];

                    $winPos++;
                }
            }
        }
    }

    1;
}

sub smooth
{
    ## so far only for windows of 5 points 
    my $len = shift;
    my $ref_rtstop = shift;

    my @smoothed = ();
    $smoothed[0] = ( $ref_rtstop->[0] + $ref_rtstop->[1] + $ref_rtstop->[2] ) / 3;
    $smoothed[1] = ( $ref_rtstop->[0] + $ref_rtstop->[1] + $ref_rtstop->[2] + $ref_rtstop->[3] ) / 4;
    for ( my $idxPos = 2; $idxPos < $len-2; $idxPos++ ) {
        $smoothed[$idxPos] = ( $ref_rtstop->[$idxPos-2] + $ref_rtstop->[$idxPos-1] + $ref_rtstop->[$idxPos] + $ref_rtstop->[$idxPos+1] + $ref_rtstop->[$idxPos+2] ) / 5;
    }
    $smoothed[$len-2] = ( $ref_rtstop->[$len-4] + $ref_rtstop->[$len-3] + $ref_rtstop->[$len-2] + $ref_rtstop->[$len-1] ) / 4;
    $smoothed[$len-1] = ( $ref_rtstop->[$len-3] + $ref_rtstop->[$len-2] + $ref_rtstop->[$len-1] ) / 3;
    
    return \@smoothed;
}
