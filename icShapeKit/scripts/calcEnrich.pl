#! /usr/bin/perl
#
# generate uniq alignment sam file
#   1, uniq hit are defined as the hit fragment starting position and the length
#
use strict;
use warnings;
use File::Basename;
use Getopt::Std;
use Data::Dumper;

use lib "/home/qczhang/lib/perllib";
use Locale::Schedule::Simple qw( &waitForFile &finishTag );

use vars qw ($opt_h $opt_V $opt_D $opt_o $opt_f $opt_b $opt_g $opt_d $opt_l $opt_e $opt_w $opt_x $opt_y $opt_A $opt_C $opt_R $opt_T $opt_S );
&getopts('hVDo:f:b:g:d:l:e:w:x:y:A:C:R:T:S:');

my $usage = <<_EOH_;
## --------------------------------------
Calculate enrichment file using RT stop as foreground and base density as background

Command:
 $0 -o output_enrichment_file -f foreground_file  -b background_file
# what it is:
 -o     output enrichment file
 -f     RT stop file (must be normalized)
 -b     base density file (must be normalized)

# more options:
 -e     enrich method, subtraction or dividing
 -w     winsor method, defined as, e.g. factor5:scaling10
 -x     subtraction factor
 -y     dividing factor

 -d     number of leading nucleotides to crop (enrichment values will be set to 1 in normal and zero in log space)
 -l     number of tailing nucleotides to crop (enrichment values will be set to 1 in normal and zero in log space)
 -g     input signals are in log rather than in normal space 

 -C     use a tag file as the agent in checking the availability of the input file
 -R     use a tag file to mark the output file ready for use
 -A     use a tag to check whether the splited fastq file is ready to use
 -T     time out
 -S     every time sleep

_EOH_
;


&main();

sub main {
    my %parameters = &init();

    my $bgStatus = waitForFile ( $parameters{backgroundFile}, $parameters{waitForTag}, $parameters{sleepTime}, $parameters{timeout} );
    die "background file $parameters{backgroundFile} not ready!\n" if ( $bgStatus != 1 );
    my $fgStatus = waitForFile ( $parameters{foregroundFile}, $parameters{waitForTag}, $parameters{sleepTime}, $parameters{timeout} );
    die "foreground file $parameters{foregroundFile} not ready!\n" if ( $fgStatus != 1 );

    my ( $ref_bg_len, $ref_bg_rpkm, $ref_bg_scalingFactor_bd, $ref_bg_baseDensity, $ref_bg_avgRTstop, $ref_bg_scalingFactor_rt, $ref_bg_RTstop ) = &readSignal ( $parameters{backgroundFile} );
    my ( $ref_fg_len, $ref_fg_rpkm, $ref_fg_scalingFactor_bd, $ref_fg_baseDensity, $ref_fg_avgRTstop, $ref_fg_scalingFactor_rt, $ref_fg_RTstop ) = &readSignal ( $parameters{foregroundFile} );
    my $ref_enrichment = &calcEnrichment ( $ref_bg_len, $ref_bg_baseDensity, $ref_bg_RTstop, $ref_fg_len, $ref_fg_baseDensity, $ref_fg_RTstop, $parameters{enrichMethod}, \%parameters );

    &winsorization ( $ref_fg_len, $ref_enrichment, $parameters{headToSkip}, $parameters{tailToSkip}, $parameters{winsorMethod} ) if ( defined $parameters{winsorMethod} );
    &outputEnrichment ( $parameters{enrichmentFile}, $ref_fg_len, $ref_enrichment, $ref_fg_rpkm, $ref_fg_scalingFactor_rt, $ref_fg_RTstop, $ref_bg_rpkm,  $ref_bg_scalingFactor_bd,$ref_bg_scalingFactor_rt, $ref_bg_baseDensity, $ref_bg_RTstop );

    finishTag ( $parameters{enrichmentFile}, $parameters{waitForTag} ) if ( defined $parameters{waitForTag} );
    finishTag ( $parameters{readyForUse} ) if ( defined $parameters{readyForUse} );

    1;
}


sub init {
    my %parameters = ();

    die $usage if ( $opt_h || ( not defined $opt_o ) || ( not defined $opt_f ) || ( not defined $opt_b ) );
    $opt_V = 0 if ( not defined $opt_V );
    $opt_D = 0 if ( not defined $opt_D );

    $parameters{enrichmentFile} = $opt_o;
    $parameters{backgroundFile} = $opt_b;
    $parameters{foregroundFile} = $opt_f;

    if ( defined $opt_e ) {
        $parameters{enrichMethod} = $opt_e;
    }
    else {
        $parameters{enrichMethod} = "complex";
    }

    if ( defined $opt_w ) {
        $parameters{winsorMethod} = $opt_w;
    }

    if ( defined $opt_d ) {
        #$parameters{headToSkip} = $opt_d;
        $parameters{headToSkip} = 0;
    }
    else {
        $parameters{headToSkip} = 0;
    }
    if ( defined $opt_l ) {
        #$parameters{tailToSkip} = $opt_l;
        $parameters{tailToSkip} = 0;
    }
    else {
        $parameters{tailToSkip} = 0;
    }

    if ( defined $opt_x ) {
        $parameters{subFactor} = $opt_x;
    }
    if ( defined $opt_y ) {
        $parameters{divFactor} = $opt_y;
    }

    if ( defined $opt_g ) {
        $parameters{logspace} = 1;
    }
    else {
        $parameters{logspace} = 0;
    }

    if ( defined $opt_A ) {
        $parameters{waitForTag} = $opt_A;
    }

    if ( defined $opt_C ) {
        $parameters{checkForAvailability} = $opt_C;
    }

    if ( defined $opt_R ) {
        $parameters{readyForUse} = $opt_R;
    }

    if ( defined $opt_T ) {
        $parameters{timeout} = $opt_T;
    }

    if ( defined $opt_S ) {
        $parameters{sleepTime} = $opt_S;
    }
    else {
        $parameters{sleepTime} = 600;
    }

    print Dumper \%parameters if ( $opt_D );
    return ( %parameters );
}


sub readSignal {
    my $signalFile = shift;

    my %trans_len = ();  my %trans_rpkm = ();  my %trans_scalingFactor_bd = ();
    my %trans_baseDensity = ();
    my %trans_avgRTstop = ();
    my %trans_scalingFactor_rt = ();
    my %trans_RTstop = ();

    my $lineCount = 0;
    open ( SG, $signalFile );
    while ( my $line = <SG> )  {
        next if ( $line =~ /^#/ );
        chomp $line;
        my ( $id, $len, $type, $rpkm, $scalingFactor, @baseDensities ) = split ( /\t/, $line );
        if ( $type eq "baseDensity" ) {
            if ( ( defined $trans_len{$id} ) and ( $trans_len{$id} != $len ) ) {
                print STDERR "Warning! $id of inconsistent transcript length!\n";
            }
            $trans_len{$id} = $len;
            $trans_rpkm{$id} = $rpkm;
            $trans_scalingFactor_bd{$id} = $scalingFactor;
            $trans_baseDensity{$id} = \@baseDensities;
        }
        elsif ( $type eq "RTstop" ) {
            if ( ( defined $trans_len{$id} ) and ( $trans_len{$id} != $len ) ) {
                print STDERR "Warning! $id of inconsistent transcript length!\n";
            }
            $trans_len{$id} = $len;
            $trans_avgRTstop{$id} = $rpkm;
            $trans_scalingFactor_rt{$id} = $scalingFactor;
            $trans_RTstop{$id} = \@baseDensities;
        }
    }
    close SG;

    return ( \%trans_len, \%trans_rpkm, \%trans_scalingFactor_bd, \%trans_baseDensity, \%trans_avgRTstop, \%trans_scalingFactor_rt, \%trans_RTstop );
}

sub calcEnrichment {
    my ( $ref_bg_len, $ref_bg_baseDensity, $ref_bg_RTstop, $ref_fg_len, $ref_fg_baseDensity, $ref_fg_RTstop, $enrichMethod, $ref_parameters ) = @_;

    my $skipHead = $ref_parameters->{headToSkip};
    my $skipTail = $ref_parameters->{tailToSkip};
    my $subFactor = ( defined $ref_parameters->{subFactor} ) ? $ref_parameters->{subFactor} : 1;
    my $divFactor = ( defined $ref_parameters->{divFactor} ) ? $ref_parameters->{divFactor} : 1;
    my %trans_enrichment = ();
    foreach my $trans ( keys %{$ref_fg_RTstop} ) {
        my @siganlEnrichment = ();
        if ( ( not defined $ref_bg_len->{$trans} ) or ( not defined $ref_bg_baseDensity->{$trans} ) or ( not defined $ref_bg_RTstop->{$trans} ) ) {
            print STDERR "Error! transcript $trans not defined in background file.\n";
            next;
        }
        elsif ( not defined $ref_fg_len->{$trans} ) {
            print STDERR "Error! transcript $trans not defined in foreground file.\n";
            next;
        }
        elsif ( $ref_bg_len->{$trans} != $ref_fg_len->{$trans} ) {
            print STDERR "Error! transcript $trans is of different length in background and foreground files.\n";
            next;
        }

        if ( $ref_parameters->{logspace} or ( $enrichMethod eq "subtraction" ) )  {
            #next if ( not defined $ref_bg_RTstop->{$trans} );
            for ( my $idx = 0; $idx < $skipHead; $idx++ )  {
                $siganlEnrichment[$idx] = "NULL";
            }
            for ( my $idx = $skipHead; $idx < ( $ref_fg_len->{$trans} - $skipTail ); $idx++ )  {
                $siganlEnrichment[$idx] = sprintf ( "%.3f", ( $ref_fg_RTstop->{$trans}[$idx] - $subFactor * $ref_bg_RTstop->{$trans}[$idx] ) );
            }
            for ( my $idx = ( $ref_fg_len->{$trans} - $skipTail ); $idx < $ref_fg_len->{$trans}; $idx++ )  {
                $siganlEnrichment[$idx] = "NULL";
            }
        }
        elsif ( $enrichMethod eq "dividing" )  {
            #next if ( not defined $ref_bg_baseDensity->{$trans} );
            for ( my $idx = 0; $idx < $skipHead; $idx++ )  {
                $siganlEnrichment[$idx] = "NULL";
            }
            for ( my $idx = $skipHead; $idx < ( $ref_fg_len->{$trans} - $skipTail ); $idx++ )  {
                if ( $ref_bg_baseDensity->{$trans}[$idx] > 0 ) {
                    $siganlEnrichment[$idx] = sprintf ( "%.3f", ( $divFactor * ( $ref_fg_baseDensity->{$trans}[$idx] ) / $ref_bg_baseDensity->{$trans}[$idx] ) );
                }
                else {
                    $siganlEnrichment[$idx] = "NULL";
                }
            }
            for ( my $idx = ( $ref_fg_len->{$trans} - $skipTail ); $idx < $ref_fg_len->{$trans}; $idx++ )  {
                $siganlEnrichment[$idx] = "NULL";
            }
        }
        else  {
            #next if ( ( not defined $ref_bg_RTstop->{$trans} ) or ( not defined $ref_bg_baseDensity->{$trans} ) );
            for ( my $idx = 0; $idx < $skipHead; $idx++ )  {
                $siganlEnrichment[$idx] = "NULL";
            }
            for ( my $idx = $skipHead; $idx < ( $ref_fg_len->{$trans} - $skipTail ); $idx++ )  {
                if ( $ref_bg_baseDensity->{$trans}[$idx] > 0 ) {
                    $siganlEnrichment[$idx] = sprintf ( "%.3f", ( $divFactor * ( $ref_fg_RTstop->{$trans}[$idx] - $subFactor * $ref_bg_RTstop->{$trans}[$idx] ) / $ref_bg_baseDensity->{$trans}[$idx] ) );
                }
                else {
                    $siganlEnrichment[$idx] = "NULL";
                }
            }
            for ( my $idx = ( $ref_fg_len->{$trans} - $skipTail ); $idx < $ref_fg_len->{$trans}; $idx++ )  {
                $siganlEnrichment[$idx] = "NULL";
            }
        }

        $trans_enrichment{$trans} = \@siganlEnrichment;
    }

    return \%trans_enrichment;
}

sub winsorization{
    my ( $ref_transcript_len, $ref_transcript_signal, $headToSkip, $tailToSkip, $winsorMethod ) = @_;
    print STDERR "Winsorization...\n\t", `date`;

    my $scaleAbove0 = 1;
    my ( $winsorFactor ) = ( $winsorMethod =~ /factor(\d+)/ );
    $winsorFactor = 5 if ( not defined $winsorFactor );
    my ( $winsorScaling ) = ( $winsorMethod =~ /scaling(\d+)/ );
    $headToSkip++; 
    foreach my $trans ( keys %{$ref_transcript_signal} )  {
        my $trimed_last = $ref_transcript_len->{$trans} - $tailToSkip - 1;      ## again, to avoid base zero
        print STDERR "$trans\t$ref_transcript_len->{$trans}\t$headToSkip\t$trimed_last\t" if ( $opt_V );
        my ( $winsorUpper, $winsorLower ) = &winsorWindow ( $ref_transcript_signal->{$trans}, $headToSkip, $trimed_last, $winsorFactor );
        print STDERR "\t$winsorUpper\t$winsorLower\n" if ( $opt_V );

        if ( ( $scaleAbove0 ) and ( $winsorLower < 0 ) ) { $winsorLower = 0; }
        if ( $winsorUpper == $winsorLower ) {
            print STDERR "Not enough resolution! Skip transcript $trans\n";
            next;
        }
        if ( defined $winsorScaling ) {
            for ( my $idx = 0; $idx < $ref_transcript_len->{$trans}; $idx++ )  {
                if ( $ref_transcript_signal->{$trans}[$idx] ne "NULL" ) {
                    if ( $ref_transcript_signal->{$trans}[$idx] > $winsorUpper ) {
                        $ref_transcript_signal->{$trans}[$idx] = $winsorUpper; 
                    }
                    elsif ( $ref_transcript_signal->{$trans}[$idx] < $winsorLower ) {
                        $ref_transcript_signal->{$trans}[$idx] = $winsorLower; 
                    }

                    $ref_transcript_signal->{$trans}[$idx] = sprintf ( "%.3f", ( ( $ref_transcript_signal->{$trans}[$idx] - $winsorLower ) *$winsorScaling/( $winsorUpper - $winsorLower )) );
                }
            }
        }
        else {
            for ( my $idx = 0; $idx < $ref_transcript_len->{$trans}; $idx++ )  {
                if ( $ref_transcript_signal->{$trans}[$idx] ne "NULL" ) {
                    if ( $ref_transcript_signal->{$trans}[$idx] > $winsorUpper ) {
                        $ref_transcript_signal->{$trans}[$idx] = $winsorUpper; 
                    }
                    elsif ( $ref_transcript_signal->{$trans}[$idx] < $winsorLower ) {
                        $ref_transcript_signal->{$trans}[$idx] = $winsorLower; 
                    }
                }
            }
        }
    }

    1;
}

sub winsorWindow  {
    my ( $ref_array, $start, $end, $winsorFactor ) = @_;
    $winsorFactor /= 100;
    ## winsorFactor is defined in pentage 
    die ( "Error! Winsor factor must be less than 0.5.\n" ) if ( $winsorFactor >= 0.5 );

    #print STDERR "$start\t$end\n";
    my @sortedTrimed_array = ();
    for ( my $idx = $start; $idx <= $end; $idx++ ) {
        push @sortedTrimed_array, $ref_array->[$idx] if ( $ref_array->[$idx] ne "NULL" );
    }
    @sortedTrimed_array = sort { $b <=> $a } @sortedTrimed_array;

    my $len = scalar(@sortedTrimed_array); 
    my $winsorLen = int ( $winsorFactor * $len );
    print STDERR "\t", $len, "\t", $winsorLen if ( $opt_V );

    return ( $sortedTrimed_array[$winsorLen], $sortedTrimed_array[$len-$winsorLen-1] );
}

sub outputEnrichment  {
    my ( $outputFile, $ref_transcript_len, $ref_transcript_enrichment, $ref_fg_rpkm, $ref_fg_scalingFactor_rt, $ref_fg_RTstop, $ref_bg_rpkm, $ref_bg_scalingFactor_bd, $ref_bg_scalingFactor_rt, $ref_bg_baseDensity, $ref_bg_RTstop ) = @_;
    print STDERR "Output enrichment scores to file $outputFile...\n\t", `date`;

    open ( OUT, ">$outputFile" ) or die "Error opening $outputFile for writing transcript base density.\n";
    print OUT "#transcript\tlength\tenrichment score, start from position 1.\n";
    foreach my $trans ( keys %{$ref_transcript_enrichment} )  {
        print OUT $trans, "\t", $ref_transcript_len->{$trans}, "\t", $ref_fg_rpkm->{$trans}, ",", $ref_bg_rpkm->{$trans}, "\t", $ref_fg_scalingFactor_rt->{$trans}, ",", $ref_bg_scalingFactor_bd->{$trans}, ",", $ref_bg_scalingFactor_rt->{$trans};
        for ( my $idx = 0; $idx < $ref_transcript_len->{$trans}; $idx++ ) {
            print OUT "\t", $ref_transcript_enrichment->{$trans}[$idx];
            print OUT ",", sprintf ( "%.3f", ( $ref_fg_scalingFactor_rt->{$trans}*$ref_fg_RTstop->{$trans}[$idx] ) );
            print OUT ",", sprintf ( "%.3f", ( $ref_bg_scalingFactor_rt->{$trans}*$ref_bg_RTstop->{$trans}[$idx] )  );
            print OUT ",", sprintf ( "%.3f", ( $ref_bg_scalingFactor_bd->{$trans}*$ref_bg_baseDensity->{$trans}[$idx] )  );
        }
        print OUT "\n";
    }
    close OUT;

    1;
}
