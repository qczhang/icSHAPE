#! /usr/bin/perl
#
use strict;
use warnings;

use lib "$ENV{ICSHAPE}/module";
use Getopt::Std;
use icSHAPEutil qw( &readGTF &getBioType &get5primeLen &get3primeLen &getExonLen );

use vars qw ($opt_h $opt_V $opt_D $opt_i $opt_o $opt_t $opt_p $opt_g $opt_f $opt_c );
&getopts('hVDi:t:o:p:g:f:c:');

my $usage = <<_EOH_;
## --------------------------------------
meta gene analysis

Command:
$0 -i transcript_profile -t heatmap_file -o metagene_profile 

# what it is:
 -i     profile for every transcript, in the format of:
        transcript_ID  [more info] position1_data  position2_data  ...

 -o     profile of meta gene
        if running in mode -p(position file), metagene is calculated in the range of -50..50; if in the mode of -g(gene annotations), metagene is calculated for all protein coding genes in the range of -50..TSS..+100,-100..TES..+50

 -t     heatmap of all transcripts

# more options:
 -p     position file (in the format of:
        transcript_ID   position
        ...

 -g     gtf annotation file (annotations from ensembl are suggested)
 -f     gtf file format (default: ensembl)

 -c     specify the first column in the input transcript profile. Default is 2, which means the profile data for that transcript starts from column 2. But sometimes, you may include other information in starting columns and your real profile starts from later columns - tell the script using this option.

_EOH_
;

&main ();

sub main {
    my %parameters = &init();

    my $transcriptFile = $parameters{transcriptFile};
    my $heatmapFile = $parameters{heatmapFile};
    my $metageneFile = $parameters{metageneFile};
    my $profileCol = $parameters{profileCol};

    my $ref_position = undef;
    my $ref_annotation = undef;
    if ( defined $parameters{positionFile} ) {
        $ref_position = readPosition ( $parameters{positionFile} );
    }
    elsif ( defined $parameters{gtfFile} ) {
        $ref_annotation = readGTF ( $parameters{gtfFile}, 
            attribute => "transcript_id", 
            source => $parameters{gtfSource}, 
            verbose => 1, 
            skip => "Blat" ); 
        &smallfix ( $ref_annotation );
    }

    my @avgNum = ( ); my @avgProfile = ( );
    if ( defined $parameters{positionFile} ) { @avgNum = newArray ( 100 ); @avgProfile = newArray ( 100 ); }
    elsif ( defined $parameters{gtfFile} ) { @avgNum = newArray ( 300 ); @avgProfile = newArray ( 300 ); }

    my $lineCount = 0;
    open ( TR, $transcriptFile );
    while ( my $line = <TR> ) {
        next if ( $line =~ /^#/ );
        $lineCount++; #last if ( $lineCount > 10 );

        chomp $line;
        my ( $transID, @data ) = split ( /\s/, $line );
        my @profile = @data[$profileCol..$#data];
        my $biotypeString = getBioType ( $transID, $ref_annotation );

        if ( defined $parameters{positionFile} ) {
            if ( defined $ref_position->{$transID} ) { &positionProfile ( $transID, \@profile, $heatmapFile, \@avgNum, \@avgProfile, $ref_position->{$transID} ); }
        }
        elsif ( defined $parameters{gtfFile} ) {
            if ( $biotypeString eq "protein coding" ) {
                if ( ( defined $ref_annotation->{$transID}{start_codon} ) and ( defined $ref_annotation->{$transID}{stop_codon} ) ) {
                    &geneProfile ( $transID, \@profile, $heatmapFile, \@avgNum, \@avgProfile, $ref_annotation );
                }
                else { print STDERR "Warning! protein coding transcript $transID has no start and stop codon annotation...skipped!\n"; }
            }
        }
    }
    close TR;

    &printAverage ( $metageneFile, \@avgNum, \@avgProfile );
}


## ------------------------------------
sub init
{
    my %parameters = ();
    die $usage if ( $opt_h || ( not $opt_i ) || ( not $opt_o ) || ( not $opt_t ) );

    $opt_V = 0 if ( not defined $opt_V );
    $opt_D = 0 if ( not defined $opt_D );
    my $pwd = `pwd`;  chomp $pwd;

    $parameters{transcriptFile} = $opt_i;
    $parameters{heatmapFile} = $opt_t;
    $parameters{metageneFile} = $opt_o;

    if ( defined $opt_p ) { $parameters{positionFile} = $opt_p; }
    if ( defined $opt_g ) { $parameters{gtfFile} = $opt_g; }
    if ( defined $opt_f ) { $parameters{gtfSource} = $opt_f; }
    else { $parameters{gtfSource} = "ensembl"; }
    if ( defined $opt_c ) { 
        die "Column out of range!\n" if ( $opt_c < 2 );
        $parameters{profileCol} = $opt_c - 2; 
    }
    else { $parameters{profileCol} = 0; }

    return ( %parameters );
}

sub newArray 
{
    my $len = shift;
    my %parameters = @_;

    my @array = ();
    for ( my $idx = 0; $idx < $len; $idx++ ) { $array[$idx] = 0; }
    return @array;
}

sub positionProfile 
{
    my $transID = shift; my $ref_profile = shift; 
    my $heatmapFile = shift;
    my $ref_avgNum = shift;  my $ref_avgProfile = shift;
    my $ref_position = shift; 

    if ( defined $heatmapFile ) {
        open ( OUT, ">>$heatmapFile" );
        print OUT $transID;
    }
    foreach my $featurePos ( keys %{$ref_position} ) {
        print STDERR "Warning! feature position out of range!\n" if ( $featurePos < 0 );

        print OUT "\t", $featurePos if ( defined $heatmapFile );
        for ( my $pos = 0; $pos <= 100; $pos++ ) {
            my $relPos = $pos + $featurePos - 50;
            if ( $relPos >= 0 ) {
                if ( $ref_profile->[$relPos] ne "NULL" ) {
                    $ref_avgNum->[$pos]++; 
                    $ref_avgProfile->[$pos] += $ref_profile->[$relPos];
                }
                print OUT "\t", $ref_profile->[$relPos] if ( defined $heatmapFile );
            }
            else { print OUT "\tNA" if ( defined $heatmapFile ); }

            print OUT "\n" if ( defined $heatmapFile );
        }
    }

    if ( defined $heatmapFile ) { close OUT;  }
}

sub geneProfile 
{
    my $transID = shift; my $ref_profile = shift; 
    my $heatmapFile = shift;
    my $ref_avgNum = shift;  my $ref_avgProfile = shift;
    my $ref_annotation = shift; 

    my $exonLen = &getExonLen ( $ref_annotation->{$transID} );
    my $fivePrimeLen = &get5primeLen ( $ref_annotation->{$transID} );
    my $threePrimeLen = &get3primeLen ( $ref_annotation->{$transID} );
    my $CDSlen = $exonLen - $fivePrimeLen - $threePrimeLen;
    my $TSSpos = $fivePrimeLen;  my $TESpos = $fivePrimeLen + $CDSlen;
    ## check strand ## note that position are 0-indexed
    if ( defined $heatmapFile ) {
        open ( OUT, ">>$heatmapFile" );
        print OUT $transID;
    }

    if ( $fivePrimeLen < 50 ) {
        print STDERR "$transID 5 prime length $fivePrimeLen less than 50!\n" if ( $opt_V );
        if ( defined $heatmapFile ) { for ( my $posIdx = 0; $posIdx < ( 50 - $fivePrimeLen ); $posIdx++ ) { print OUT "\tNULL"; } }
        for ( my $posIdx = 0; $posIdx < $TSSpos; $posIdx++ ) {
            my $value = "NULL";
            if ( defined $ref_profile->[$posIdx] ) { $value = $ref_profile->[$posIdx]; }
            else { print STDERR "Warning! Inconsistent annotation? Position $posIdx of $transID not defined!\n"; }
            if ( defined $heatmapFile ) { print OUT "\t", $value; }
            if ( $value ne "NULL" ) {
                $ref_avgNum->[50-$fivePrimeLen+$posIdx]++;
                $ref_avgProfile->[50-$fivePrimeLen+$posIdx] += $value;
            }
        }
    }
    else {
        for ( my $posIdx = ( $TSSpos - 50 ); $posIdx < $TSSpos; $posIdx++ ) {
            my $value = "NULL";
            if ( defined $ref_profile->[$posIdx] ) { $value = $ref_profile->[$posIdx]; }
            else { print STDERR "Warning! Inconsistent annotation? Position $posIdx of $transID not defined!\n"; }
            if ( defined $heatmapFile ) { print OUT "\t", $value; }
            if ( $value ne "NULL" ) {
                $ref_avgNum->[50-$fivePrimeLen+$posIdx]++;
                $ref_avgProfile->[50-$fivePrimeLen+$posIdx] += $value;
            }
        }
    }

    if ( $CDSlen < 100 ) {
        print STDERR "$transID CDS length $CDSlen less than 50!\n" if ( $opt_V );
        for ( my $posIdx = $TSSpos; $posIdx < ( $TSSpos + $CDSlen ); $posIdx++ ) {
            my $value = "NULL";
            if ( defined $ref_profile->[$posIdx] ) { $value = $ref_profile->[$posIdx]; }
            else { print STDERR "Warning! Inconsistent annotation? Position $posIdx of $transID not defined!\n"; }
            if ( defined $heatmapFile ) { print OUT "\t", $value; }
            if ( $value ne "NULL" ) {
                $ref_avgNum->[50-$fivePrimeLen+$posIdx]++;
                $ref_avgProfile->[50-$fivePrimeLen+$posIdx] += $value;
            }
        }
        if ( defined $heatmapFile ) {
            for ( my $posIdx = 0; $posIdx < ( 100 - $CDSlen ); $posIdx++ ) { print OUT "\tNULL"; }
            for ( my $posIdx = 0; $posIdx < ( 100 - $CDSlen ); $posIdx++ ) { print OUT "\tNULL"; }
        }
        for ( my $posIdx = ( $TESpos - 100 ); $posIdx < $TESpos; $posIdx++ ) {
            my $value = "NULL";
            if ( defined $ref_profile->[$posIdx] ) { $value = $ref_profile->[$posIdx]; }
            else { print STDERR "Warning! Inconsistent annotation? Position $posIdx of $transID not defined!\n"; }
            if ( defined $heatmapFile ) { print OUT "\t", $value; }
            if ( $value ne "NULL" ) {
                $ref_avgNum->[250-$fivePrimeLen-$CDSlen+$posIdx]++;
                $ref_avgProfile->[250-$fivePrimeLen-$CDSlen+$posIdx] += $value;
            }
        }
    }
    else {
        for ( my $posIdx = $TSSpos; $posIdx < ( $TSSpos + 100 ); $posIdx++ ) {
            my $value = "NULL";
            if ( defined $ref_profile->[$posIdx] ) { $value = $ref_profile->[$posIdx]; }
            else { print STDERR "Warning! Inconsistent annotation? Position $posIdx of $transID not defined!\n"; }
            if ( defined $heatmapFile ) { print OUT "\t", $value; }
            if ( $value ne "NULL" ) {
                $ref_avgNum->[50-$fivePrimeLen+$posIdx]++;
                $ref_avgProfile->[50-$fivePrimeLen+$posIdx] += $value;
            }
        }
        for ( my $posIdx = ( $TESpos - 100 ); $posIdx < $TESpos; $posIdx++ ) {
            my $value = "NULL";
            if ( defined $ref_profile->[$posIdx] ) { $value = $ref_profile->[$posIdx]; }
            else { print STDERR "Warning! Inconsistent annotation? Position $posIdx of $transID not defined!\n"; }
            if ( defined $heatmapFile ) { print OUT "\t", $value; }
            if ( $value ne "NULL" ) {
                $ref_avgNum->[250-$fivePrimeLen-$CDSlen+$posIdx]++;
                $ref_avgProfile->[250-$fivePrimeLen-$CDSlen+$posIdx] += $value;
            }
        }
    }

    if ( $threePrimeLen < 50 ) {
        print STDERR "$transID 3 prime length $threePrimeLen less than 50!\n" if ( $opt_V );
        for ( my $posIdx = $TESpos; $posIdx < ( $TESpos + $threePrimeLen ); $posIdx++ ) {
            my $value = "NULL";
            if ( defined $ref_profile->[$posIdx] ) { $value = $ref_profile->[$posIdx]; }
            else { print STDERR "Warning! Inconsistent annotation? Position $posIdx of $transID not defined!\n"; }
            if ( defined $heatmapFile ) { print OUT "\t", $value; }
            if ( $value ne "NULL" ) {
                $ref_avgNum->[250-$fivePrimeLen-$CDSlen+$posIdx]++;
                $ref_avgProfile->[250-$fivePrimeLen-$CDSlen+$posIdx] += $value;
            }
        }
        if ( defined $heatmapFile ) {
            for ( my $posIdx = 0; $posIdx < ( 50 - $threePrimeLen ); $posIdx++ ) { print OUT "\tNULL"; }
        }
    }
    else {
        for ( my $posIdx = $TESpos; $posIdx < ( $TESpos + 50 ); $posIdx++ ) {
            my $value = "NULL";
            if ( defined $ref_profile->[$posIdx] ) { $value = $ref_profile->[$posIdx]; }
            else { print STDERR "Warning! Inconsistent annotation? Position $posIdx of $transID not defined!\n"; }
            if ( defined $heatmapFile ) { print OUT "\t", $value; }
            if ( $value ne "NULL" ) {
                $ref_avgNum->[250-$fivePrimeLen-$CDSlen+$posIdx]++;
                $ref_avgProfile->[250-$fivePrimeLen-$CDSlen+$posIdx] += $value;
            }
        }
    }

    if ( defined $heatmapFile ) { print OUT "\n"; close OUT;  }
}

sub printAverage 
{
    my $outFile = shift;
    my $ref_numArray = shift;
    my $ref_avgArray = shift;

    open ( OUT, ">$outFile" ) or die ( "Could not open file $outFile for reading!\n" );
    print OUT "count";
    for ( my $idx = 0; $idx < scalar ( @{$ref_numArray} ); $idx++ ) {
        print OUT "\t", $ref_numArray->[$idx];
    }
    print OUT "\n";

    print OUT "average";
    for ( my $idx = 0; $idx < scalar ( @{$ref_numArray} ); $idx++ ) {
        if ( $ref_numArray->[$idx] ) {
            $ref_avgArray->[$idx] /= $ref_numArray->[$idx];
        }
        else {
            $ref_avgArray->[$idx] = "NULL";
        }
        print OUT "\t", $ref_avgArray->[$idx];
    }
    print OUT "\n";
    close OUT;
}

sub readPosition
{
    my $file = shift;
    my %parameters = @_;

    open ( PS, $file ) or die ( "Cannot open file $file of position information for reading!\n" );
    print STDERR "read in position file $file\n\t", `date`;
    my %feature_pos = ();
    while ( my $line = <PS> ) {
        next if ( $line =~ /^#/ );
        chomp $line;
        my @data = split ( /\t/, $line );

        my $transID = $data[0];
        my $pos = $data[1];
        $feature_pos{$transID}{$pos} = 1;
    }
    close PS;

    return ( \%feature_pos );
}

sub smallfix
{
    my $ref_annotation = shift;
    $ref_annotation->{MM5_8SrRNA}{biotype} = [ "rRNA gene" ]; $ref_annotation->{MM18SrRNA}{biotype} = [ "rRNA gene" ]; $ref_annotation->{MM28SrRNA}{biotype} = [ "rRNA gene" ];
    $ref_annotation->{"lcl|gene=robert:MM18SrRNA"}{biotype} = [ "rRNA gene" ]; $ref_annotation->{"lcl|gene=robert:MM28SrRNA"}{biotype} = [ "rRNA gene" ];
}
