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

use vars qw ($opt_h $opt_V $opt_D $opt_i $opt_o $opt_t $opt_T $opt_s $opt_e );
&getopts('hVDi:o:t:T:s:e:');

my $usage = <<_EOH_;
## --------------------------------------
Calculate average signals

Command:
 $0 -i input_signal_file 
# what it is:
 -i     input signal file 

# more options:
 -t     threshold to select window for correlation calculation
 -T     threshold of coverage
 -s     head to skip
 -e     end to skip

_EOH_
;


&main();

sub main 
{
    my %parameters = &init();

    my ( $ref_trans_len1, $ref_trans_rpkm1, $ref_trans_rtstop1 ) = &readSignals ( $parameters{input1}, 
        outputFile => $parameters{output},
        coverageThreshold => $parameters{coverageThreshold},
        winThreshold => $parameters{windowThreshold},
        headToSkip => $parameters{headToSkip},
        tailToSkip => $parameters{tailToSkip}
    );

    1;
}


sub init 
{
    my %parameters = ();

    die $usage if ( $opt_h || ( not defined $opt_i ) );
    $opt_V = 0 if ( not defined $opt_V );
    $opt_D = 0 if ( not defined $opt_D );

    $parameters{input1} = $opt_i;
    $parameters{output} = $opt_o;

    if ( defined $opt_t )  { $parameters{windowThreshold} = $opt_t; }
    else { $parameters{windowThreshold} = 0;  }
    if ( defined $opt_T )  { $parameters{coverageThreshold} = $opt_T; }
    else { $parameters{coverageThreshold} = 0;  }
    if ( defined $opt_s )  { $parameters{headToSkip} = $opt_s; }
    else { $parameters{headToSkip} = 0;  }
    if ( defined $opt_e )  { $parameters{tailToSkip} = $opt_e; }
    else { $parameters{tailToSkip} = 0;  }

    print STDERR Dumper ( \%parameters ) if ( $opt_D );
    return ( %parameters );
}

sub readSignals 
{
    my $signalFile = shift;
    my %parameters = @_;

    my $winThreshold = $parameters{winThreshold};
    my $coverageThreshold = $parameters{coverageThreshold};
    my $headToSkip = $parameters{headToSkip};
    my $tailToSkip = $parameters{tailToSkip};

    if ( defined $parameters{outputFile} ) { open ( OUT, ">$parameters{outputFile}" ) or die "Error opening $parameters{outputFile} for writing!\n"; }
    print STDERR "read signal from $signalFile\n", `date`;
    my $line = "";  my $lineCount = 0;
    open ( SG, $signalFile );
    while ( $line = <SG> ) {
        next if ( $line =~ /^#/ );
        $lineCount++;
        print STDERR "  line: $lineCount\n" if ( $lineCount % 10000 == 0 );

        my ( $id, $len, $rpkmString, $scalingFactors, @baseDensities ) = split ( /\t/, $line );
        my @scores = ();  my @fgrtstops = (); my @bgBDs = ();
        #my ( $scalingFactor ) = split ( /,/, $scalingFactors );
        for ( my $idxPos = 0; $idxPos < scalar(@baseDensities); $idxPos++ ) {
            my ( $score, $fgrtstop, $bgrtstop, $bgBD ) = split ( /,/, $baseDensities[$idxPos] );
            push @scores, $score;
            push @fgrtstops, $fgrtstop if ( $fgrtstop ne "NULL" ); 
            push @bgBDs, $bgBD if ( $bgBD ne "NULL" );
        }

        my $coverage = eval ( join ( "+", @fgrtstops ) ) / scalar ( @fgrtstops );
        if ( $coverage > $coverageThreshold ) {
            my @rpkms = split ( /,/, $rpkmString );
            my $rpkm = eval ( join ( "+", @rpkms ) ) / scalar ( @rpkms );
            if ( defined $parameters{outputFile} ) {  print OUT $id, "\t", $len, "\t", $rpkm;  }
            else {  print $id, "\t", $len, "\t", $rpkm;  }

            for ( my $idxPos = 0; $idxPos < $len; $idxPos++ ) {
                my $selected = 0;
                
                if ( ( $idxPos > $headToSkip ) and ( $idxPos < ( $len - $tailToSkip ) ) )  {
                    if ( ( $bgBDs[$idxPos] ne "NULL" ) and ( $bgBDs[$idxPos] >= $winThreshold ) )  { $selected = 1; }
                }

                if ( defined $parameters{outputFile} ) {  
                    if ( $selected ) { print OUT "\t", $scores[$idxPos]; } else { print OUT "\tNULL";  }
                }
                else {
                    if ( $selected ) { print "\t", $scores[$idxPos]; } else { print "\tNULL";  }
                }
            }
            if ( defined $parameters{outputFile} ) {  print OUT "\n";  }
            else {  print "\n";  }
        }
    }
    close SG;

    1;
}


