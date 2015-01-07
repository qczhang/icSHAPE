#! /usr/bin/perl
#
#
use strict;
use warnings;
use File::Basename;
use Getopt::Std;

use vars qw ($opt_h $opt_V $opt_D $opt_i $opt_o $opt_r $opt_c );
&getopts('hVDi:o:r:c:');

my $usage = <<_EOH_;
## --------------------------------------
calculate RT stops from sam file

Command:
$0 -i input_sam_file -o output_RTstop_file -r transcript_rpkm_file

# what it is:
 -i     input sam file
 -o     output RTstop file (with base density information)
 -r     rpkm file

# more options:
 -c     cutoff of RPKM

_EOH_
    ;


&main();

sub main {
    my %parameters = &init();

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
    foreach my $samFile ( @samFiles ) { &calcBaseDensity ( $samFile, \%baseDensity, \%RTstop, $ref_transcript_len, $ref_transcript_FPKM ); }

    &output_baseDensity ( $parameters{output}, $ref_transcript_len, $ref_transcript_FPKM, \%baseDensity, \%RTstop );
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

    if ( defined $opt_c ) { $parameters{minload} = $opt_c; }
    else { $parameters{minload} = 0; }

    return ( %parameters );
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
        if ( ( $tag == 99 ) or ( $tag == 355 ) ) { next if ( $rnext ne "=" ); } ## PE
        elsif ( ( $tag == 0 ) or ( $tag == 256 ) ) { $tlen = length ( $seq ); next if ( not $tlen ); } ## SE
        else { next; }

        if ( $read ne $readID ) {
            if ( $readID ) {
                foreach my $hitID ( keys %hitID_start ) {
                    next if ( not defined $ref_transcript_len->{$hitID} );

                    for ( my $idx = $hitID_start{$hitID}; $idx < $hitID_end{$hitID}; $idx++ )
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
        print OUT $seq, "\t", $ref_transcript_len->{$seq}, "\t", $ref_transcript_rpkm->{$seq};
        for ( my $idx = 0; $idx <= $ref_transcript_len->{$seq}; $idx++ ) { print OUT "\t", sprintf ( "%.3f", $ref_baseDensity->{$seq}[$idx] ); }
        print OUT "\n";

        print OUT $seq, "\t", $ref_transcript_len->{$seq}, "\t", $ref_transcript_rpkm->{$seq};
        for ( my $idx = 0; $idx <= $ref_transcript_len->{$seq}; $idx++ ) { print OUT "\t", sprintf ( "%.3f", $ref_RTstop->{$seq}[$idx] ); }
        print OUT "\n";
    }
    close OUT;

    1;
}
