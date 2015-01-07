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

use vars qw ($opt_h $opt_V $opt_D $opt_i $opt_o $opt_r $opt_p $opt_c $opt_A $opt_C $opt_R $opt_T $opt_S );
&getopts('hVDi:o:r:pc:A:C:R:T:S:');

my $usage = <<_EOH_;
## --------------------------------------
calculate expression values for all transcripts in FPKM from SAM files 
  (attention: now only for PE reads)

Command:
$0 -i input_sam_file -o output_baseDensity_file -r transcript_rpkm_file

# what it is:
 -i     input_sam_file
 -o     output_baseDensity_file
 -r     rpkm_file

# more options:
 -p     calculate RT stops
 -c     cutoff of RPKM

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

    my $RPKMStatus = waitForFile ( $parameters{rpkmFile}, $parameters{waitForTag}, $parameters{sleepTime}, $parameters{timeout} );
    die "Reference RPKM file not ready!\n" if ( $RPKMStatus != 1 );

    my ( $ref_transcript_len, $ref_transcript_FPKM ) = &readRPKM ( $parameters{rpkmFile}, $parameters{minload} );

    my %baseDensity = ();
    my %RTstop = ();
    foreach my $tran ( keys %{$ref_transcript_len} ) {
        for ( my $idx = 0; $idx <= $ref_transcript_len->{$tran}; $idx++ ) {
            $baseDensity{$tran}[$idx] = 0;
            $RTstop{$tran}[$idx] = 0;
        }
    }

    my @samFiles = split ( /:/, $parameters{input} );
    foreach my $samFile ( @samFiles ) {
        my $samStatus = waitForFile ( $samFile );
        die "Input SAM file not ready!\n" if ( $samStatus != 1 );

        &calcBaseDensity ( $samFile, \%baseDensity, \%RTstop, $ref_transcript_len, $ref_transcript_FPKM );
    }

    &output_baseDensity ( $parameters{output}, $ref_transcript_len, $ref_transcript_FPKM, \%baseDensity, \%RTstop );
    finishTag ( $parameters{output}, $parameters{waitForTag} ) if ( defined $parameters{waitForTag} );
    finishTag ( $parameters{readyForUse} ) if ( defined $parameters{readyForUse} );
}

sub init {
    print STDERR "Initializing...\n\t", `date`;

    my %parameters = ();

    die $usage if ( $opt_h || ( not defined $opt_i ) || ( not defined $opt_o ) );
    $opt_V = 0 if ( not defined $opt_V );
    $opt_D = 0 if ( not defined $opt_D );

    $parameters{input} = $opt_i;
    $parameters{output} = $opt_o;
    $parameters{rpkmFile} = $opt_r;

    if ( defined $opt_p ) {
        $parameters{isRTstop} = 1;
    }
    else {
        $parameters{isRTstop} = 0;
    }

    if ( defined $opt_c ) {
        $parameters{minload} = $opt_c;
    }
    else {
        $parameters{minload} = 0;
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

# deprecated
sub _waitForRPKM  {
    my $ref_parameters = shift;
    print STDERR "Checking the availability of $ref_parameters->{rpkmFile}...\n\t", `date`;

    my $RPKMFileStatus = 1;
    if ( $ref_parameters->{waitForTag} ) {
        my $sleep = $ref_parameters->{sleepTime};
        my $tryTag1 = $ref_parameters->{rpkmFile} . "." . $ref_parameters->{waitForTag};
        my $tryTag2 = $ref_parameters->{rpkmFile};  $tryTag2 =~ s/rpkm/$ref_parameters->{waitForTag}/;

        while ( ( not -e $tryTag1 ) && ( not -e $tryTag2 ) ) {
            print STDERR "RPKM file $ref_parameters->{rpkmFile} is not ready yet, waiting...\n\t", `date`;

            if ( defined $ref_parameters->{timeout} )  {
                if ( $ref_parameters->{timeout} <= 0 ) {
                    print STDERR "ERROR! Time out.\n\t", `date`;
                    $RPKMFileStatus = -1;
                    last;
                }
                elsif ( $ref_parameters->{timeout} < $sleep ) { 
                    $sleep = $ref_parameters->{timeout};
                    $ref_parameters->{timeout} = -1;
                }
                else  {
                    $ref_parameters->{timeout} -= $sleep;
                }
            }

            sleep $sleep;
        }
    }
    else {
        if ( not -e $ref_parameters->{rpkmFile} ) {
            $RPKMFileStatus = -1;
        }
        elsif ( -z $ref_parameters->{rpkmFile} ) {
            $RPKMFileStatus = 0;
        }
    }

    return $RPKMFileStatus;
}


sub readRPKM  {
    my $rpkmFile = shift;
    my $minload = shift;
    print STDERR "Read transcript abundance information from file $rpkmFile...\n\t", `date`;

    my %transcript_len = ();
    my %transcript_RPKM = ();
    open ( IN, $rpkmFile );
    while ( my $line = <IN> )  {
        next if ( $line =~ /^#/ );

        chomp $line;
        my ( $transcript, $len, $uniqRead, $multiRead, $rpkm ) = split ( /\t/, $line );
        next if ( $minload && ( $rpkm < $minload ) );
        $transcript_len{$transcript} = $len;
        $transcript_RPKM{$transcript} = $rpkm;
    }
    close IN;

    return ( \%transcript_len, \%transcript_RPKM );
}


sub calcBaseDensity {
    my ( $inputSamFile, $ref_baseDensity, $ref_RTstop, $ref_transcript_len, $ref_transcript_FPKM ) = @_;
    print STDERR "Calculate base density from file $inputSamFile...\n\t", `date`;

    my %hitID_start = ();  my %hitID_end = (); 
    my $readID = "";  my $hitID = "";  my $hitCount = 0;  my $totalCount = 0;
    open ( SAM, $inputSamFile );
    while ( my $line = <SAM> ) {
        next if ( $line =~ /^@/ );

        $totalCount++;
        print STDERR $totalCount, "\n" if ( $totalCount % 1000000 == 0 );

        my ( $read, $tag, $hit, $pos, $mapq, $cigar, $rnext, $pnext, $tlen, $seq, $qual ) = split ( /\t/, $line );
        if ( ( $tag == 99 ) or ( $tag == 355 ) ) {
            next if ( $rnext ne "=" );
            ## PE
        }
        elsif ( ( $tag == 0 ) or ( $tag == 256 ) ) {
            $tlen = length ( $seq );
            next if ( not $tlen );
            ## SE
        }
        else {
            next;
        }

        if ( $read ne $readID ) {
            if ( $readID ) {
                foreach my $hitID ( keys %hitID_start ) {
                    next if ( not defined $ref_transcript_len->{$hitID} );

                    for ( my $idx = $hitID_start{$hitID}; $idx < $hitID_end{$hitID}; $idx++ )
                    #{   $ref_baseDensity->{$hitID}[$idx] = ( $ref_baseDensity->{$hitID}[$idx] + $ref_transcript_FPKM->{$hitID}/$hitCount );  }
                    #$ref_RTstop->{$hitID}[$hitID_start{$hitID}-1] = ( $ref_RTstop->{$hitID}[$hitID_start{$hitID}-1] + $ref_transcript_FPKM->{$hitID}/$hitCount );
                    {   $ref_baseDensity->{$hitID}[$idx] = ( $ref_baseDensity->{$hitID}[$idx] + 1/$hitCount );  }
                    $ref_RTstop->{$hitID}[$hitID_start{$hitID}-1] = ( $ref_RTstop->{$hitID}[$hitID_start{$hitID}-1] + 1/$hitCount );
                }
            }
            $hitCount = 0;
            %hitID_start = ();
            %hitID_end = ();
        }

        $readID = $read;
        $hitCount++;
        #$hitCount = $hitCount + $ref_transcript_FPKM->{$hit};
        if ( defined $ref_transcript_len->{$hit} ) {
            $hitID_start{$hit} = $pos;
            $hitID_end{$hit} = $pos + $tlen;
            $hitID_end{$hit} = $ref_transcript_len->{$hit} + 1 if ( $ref_transcript_len->{$hit} < $hitID_end{$hit} );
        }
    }
    close SAM;

    if ( $readID ) {
        foreach my $hitID ( keys %hitID_start ) {
            for ( my $idx = $hitID_start{$hitID}; $idx < $hitID_end{$hitID}; $idx++ )
            #{   $ref_baseDensity->{$hitID}[$idx] = ( $ref_baseDensity->{$hitID}[$idx] + $ref_transcript_FPKM->{$hitID}/$hitCount );  }
            #$ref_RTstop->{$hitID}[$hitID_start{$hitID}-1] = ( $ref_RTstop->{$hitID}[$hitID_start{$hitID}-1] + $ref_transcript_FPKM->{$hitID}/$hitCount );
            {   $ref_baseDensity->{$hitID}[$idx] = ( $ref_baseDensity->{$hitID}[$idx] + 1/$hitCount );  }
            $ref_RTstop->{$hitID}[$hitID_start{$hitID}-1] = ( $ref_RTstop->{$hitID}[$hitID_start{$hitID}-1] + 1/$hitCount );
        }
    }

    1;
}

sub output_baseDensity  {
    my ( $outputFile, $ref_transcript_len, $ref_transcript_rpkm, $ref_baseDensity, $ref_RTstop ) = @_;
    print STDERR "Output base density to file $outputFile...\n\t", `date`;
    
    open ( OUT, ">$outputFile" ) or die "Error opening $outputFile for writing transcript base density.\n";
    print OUT "#transcript\tbase frequency, start from position 0.\n";
    foreach my $seq ( keys %{$ref_baseDensity} )  {
        #my $totalBaseDensity = eval ( join ( "+", @{$ref_baseDensity->{$seq}} ) );
        #my $avgBaseDensity = sprintf ( "%.4f", $totalBaseDensity / $ref_transcript_len->{$seq} );
        #print OUT $seq, "\t", $ref_transcript_len->{$seq}, "\t", $avgBaseDensity;
        print OUT $seq, "\t", $ref_transcript_len->{$seq}, "\t", $ref_transcript_rpkm->{$seq};
        for ( my $idx = 0; $idx <= $ref_transcript_len->{$seq}; $idx++ ) {
            print OUT "\t", sprintf ( "%.3f", $ref_baseDensity->{$seq}[$idx] );
        }
        print OUT "\n";

        #my $totalRTstop = eval ( join ( "+", @{$ref_RTstop->{$seq}} ) );
        #my $avgRTstop = sprintf ( "%.4f", $totalRTstop / $ref_transcript_len->{$seq} );
        #print OUT $seq, "\t", $ref_transcript_len->{$seq}, "\t", $avgRTstop;
        print OUT $seq, "\t", $ref_transcript_len->{$seq}, "\t", $ref_transcript_rpkm->{$seq};
        for ( my $idx = 0; $idx <= $ref_transcript_len->{$seq}; $idx++ ) {
            print OUT "\t", sprintf ( "%.3f", $ref_RTstop->{$seq}[$idx] );
        }
        print OUT "\n";
    }
    close OUT;

    1;
}

# deprecated
sub _finishTag  {
    my $ref_parameters = shift;

    if ( $ref_parameters->{waitForTag} ) {
        my $finishTag = $ref_parameters->{output} . "." . $ref_parameters->{waitForTag};
        print STDERR `touch $finishTag`;
    }

    1;
}
