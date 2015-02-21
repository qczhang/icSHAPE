#!/usr/bin/perl
use strict;
use warnings;

use lib "$ENV{ICSHAPE}/module";
use Getopt::Std;
use icSHAPEutil qw( &readGTF_ensembl_new );

use lib "/home/qczhang/lib/perllib/";
use Data::Dumper;

use vars qw ($opt_h $opt_V $opt_D $opt_i $opt_o $opt_a $opt_f $opt_b $opt_s );
&getopts('hVDi:a:o:f:b:s');

my $usage = <<_EOH_;
## --------------------------------------
convert genomic coordinates to transcript coordinates

Command:
$0 -i genome_coordinates -o transcript_coordinates -a annotation_file

# what it is:
 -i     genomic coordinates
 -o     transcript coordinates
 -a     annotation file (GTF format)

 -f     annotation source (e.g., ensembl, gencode, ...)
 -b     bin width 
 -s     skip nonoverlapped regions

_EOH_
;

&main ();

sub main {
    my %parameters = &init();

    my $genomeFile = $parameters{genomeFile};
    my $transcriptFile = $parameters{transcriptFile};
    my $annotationFile = $parameters{annotationFile};

    my $ref_position = readGenomePosition ( $genomeFile, format => "bed" );
    my $ref_annotation = readGTF_ensembl_new ( $parameters{annotationFile} );

    my $ref_bin = binize ( $ref_annotation->{gene_info}, $ref_annotation->{chr_size} );

    &convert ( $ref_position, $ref_annotation, $ref_bin, 1000000 );
    &printPosition ( $ref_position, $transcriptFile, skipNonOverlap => $parameters{skipNonOverlap} );

    1;
}


## ------------------------------------
sub init
{
    my %parameters = ();
    die $usage if ( $opt_h || ( not $opt_i ) || ( not $opt_o ) || ( not $opt_a ) );

    $opt_V = 0 if ( not defined $opt_V );
    $opt_D = 0 if ( not defined $opt_D );
    my $pwd = `pwd`;  chomp $pwd;

    $parameters{genomeFile} = $opt_i;
    $parameters{transcriptFile} = $opt_o;
    $parameters{annotationFile} = $opt_a;

    if ( defined $opt_f ) { $parameters{annotationSource} = $opt_f; }
    else { $parameters{annotationSource} = "ensembl"; }
    if ( defined $opt_b ) { $parameters{bw} = $opt_b; }
    else { $parameters{bw} = 1; }
    if ( defined $opt_s ) { $parameters{skipNonOverlap} = 1; }
    else { $parameters{skipNonOverlap} = 0; }

    return ( %parameters );
}

sub readGenomePosition
{
    my $input_file = shift;
    my %parameters = @_;

    open (IN, "<$input_file") or die ( "Cannot open $input_file for reading!\n" );
    print STDERR "Read positions from file $input_file.\n\t", `date`;
    my @absPosition = ();
    my $count = 0;
    while (my $line = <IN> )  {
        next if ( $line =~ /^#/ );

        chomp $line;
        my ( $chr, $start, $end, $strand, @data ) = split ( /\t/, $line );
        $chr =~ s/^chr//;

        my $format = ( defined $parameters{format} ) ? "$parameters{format}" : "bed";
        if ( $format eq "bed" ) { 
            $absPosition[$count]{chr} = $chr; $absPosition[$count]{strand} = $strand;
            $absPosition[$count]{start} = $start; $absPosition[$count]{end} = $end;
        }
        elsif ( $format eq "simple" ) {
            $absPosition[$count]{chr} = $chr; $absPosition[$count]{strand} = "+";
            $absPosition[$count]{start} = $start; $absPosition[$count]{end} = $end;
            $count++;
            $absPosition[$count]{chr} = $chr; $absPosition[$count]{strand} = "-";
            $absPosition[$count]{start} = $start; $absPosition[$count]{end} = $end;
        }
        elsif ( $format eq "single" ) {
            $absPosition[$count]{chr} = $chr; $absPosition[$count]{strand} = "+";
            $absPosition[$count]{start} = $start; $absPosition[$count]{end} = $start+1;
            $count++;
            $absPosition[$count]{chr} = $chr; $absPosition[$count]{strand} = "-";
            $absPosition[$count]{start} = $start; $absPosition[$count]{end} = $start+1;
        }
        $count++;
    }
    close IN;
    print STDERR "$count positions in total.\n\t", `date`;

    return \@absPosition;
}

sub binize
{
    my $ref_featurePos = shift;
    my $ref_chr_size = shift;
    my %parameters = @_;

    print STDERR "Binize genome to speed up searching.\n\t", `date`;
    my $bw = ( defined $parameters{bw} ) ? $parameters{bw} : 1000000;
    my %bin = ();
    foreach my $chr ( keys %{$ref_chr_size} ) {
        my $count = int ( $ref_chr_size->{$chr} / $bw ) + 1;
        for ( my $idx = 0; $idx < $count; $idx++ ) { $bin{$chr}{"+"}[$idx] = ();  $bin{$chr}{"-"}[$idx] = ();  }
    }

    foreach my $featureID ( keys %{$ref_featurePos} ) {
        my $chr = $ref_featurePos->{$featureID}{chr};
        my $strand = $ref_featurePos->{$featureID}{strand};
        my $start = int ( $ref_featurePos->{$featureID}{start} / $bw );
        my $end = int ( $ref_featurePos->{$featureID}{end} / $bw );
        for ( my $idx = $start; $idx <= $end; $idx++ ) { push @{$bin{$chr}{$strand}[$idx]}, $featureID; }
    }

    return \%bin;
}

sub convert
{
    my $ref_positions = shift;
    my $ref_annotation = shift;
    my $ref_bin = shift;
    my $bw = shift;

    print STDERR "Checking overlapped transcripts in input positions...\n\t", `date`;
    for ( my $idx = 0; $idx < scalar ( @{$ref_positions} ); $idx++ ) {
        print STDERR "position $idx\n\t", `date` if ( $idx and ( $idx % 10000 == 0 ) ); 

        my @overlapRegion = ();
        my $chr = $ref_positions->[$idx]{chr};
        my $strand = $ref_positions->[$idx]{strand};
        my $start = $ref_positions->[$idx]{start};
        my $end = $ref_positions->[$idx]{end}; 
        for ( my $idxBin = int ( $start / $bw ); $idxBin <= int ( $end / $bw ); $idxBin++ ) {
            ## get genes in the bin
            foreach my $gene ( @{$ref_bin->{$chr}{$strand}[$idxBin]} ) {
                # get transcripts for each gene
                next if ( ( $end < $ref_annotation->{gene_info}{$gene}{start} ) or ( $start > $ref_annotation->{gene_info}{$gene}{end} ) );
                foreach my $transID ( @{$ref_annotation->{gene_info}{$gene}{transcript}} ) {
                    next if ( ( $end < $ref_annotation->{transcript_info}{$transID}{start} ) or ( $start > $ref_annotation->{transcript_info}{$transID}{end} ) );
                    &overlapTrans ( \@overlapRegion, $ref_annotation, $transID, $strand, $start, $end );
                    $ref_positions->[$idx]{overlap} = \@overlapRegion;
                }
            }
        }
    }

    1;
}

sub printPosition
{
    my $ref_positions = shift;
    my $outFile = shift;
    my %parameters = @_;

    open ( OUT, ">$outFile" ) or die "Cannot open $outFile to write!\n";
    print STDERR "Output overlapped transcripts to $outFile.\n\t", `date`;
    for ( my $idx = 0; $idx < scalar ( @{$ref_positions} ); $idx++ ) {
        if ( $ref_positions->[$idx]{overlap} ) {
            foreach my $overlap ( @{$ref_positions->[$idx]{overlap}} ) 
                { print OUT join ( "\t", $ref_positions->[$idx]{chr}, $ref_positions->[$idx]{start}, $ref_positions->[$idx]{end}, $ref_positions->[$idx]{strand}, $overlap), "\n"; }
        }
        elsif ( not $parameters{skipNonOverlap} ) 
            { print OUT join ( "\t", $ref_positions->[$idx]{chr}, $ref_positions->[$idx]{start}, $ref_positions->[$idx]{end}, $ref_positions->[$idx]{strand}, ".", ".", ".", ".", "."), "\n"; }
    }

    1;
}

sub overlapTrans
{
    my $ref_overlapRegion = shift;
    my $ref_annotation = shift; my $transID = shift;
    my $strand = shift; my $absStart = shift; my $absEnd = shift; 

    my $transcript = $ref_annotation->{transcript_info}{$transID};
    my $numExon = scalar ( keys %{$transcript->{exon}} );

    my $relExonStart = 0;  my $startPosiInExon = 0; my $endPosiInExon = 0; my $absStartInExon = 0; my $absEndInExon = 0; 
    for ( my $idxExon = 1; $idxExon <= $numExon; $idxExon++ ) {
        my $idxExonWithStrand = ( $strand eq "+" ) ? $idxExon : ( $numExon - $idxExon + 1 );
        my $exonID = $ref_annotation->{transcript_info}{$transID}{exon}{$idxExon};
        my $exonStart = $ref_annotation->{exon_info}{$exonID}{start};
        my $exonEnd = $ref_annotation->{exon_info}{$exonID}{end};
        my $exonLength = $exonEnd - $exonStart + 1;

        if ( ( $exonStart <= $absEnd ) && ( $exonEnd >= $absStart ) ) {
            ## overlapped
            if ( $strand eq "+" ) {
                $startPosiInExon = $absStart - $exonStart + 1;
                $startPosiInExon = 1 if ( $startPosiInExon < 1 );
                $endPosiInExon = $absEnd - $exonStart + 1;
                $endPosiInExon = $exonLength if ( $endPosiInExon > $exonLength );
            }
            elsif ( $strand eq "-" ) {
                $startPosiInExon = $exonEnd - $absEnd + 1;
                $startPosiInExon = 1 if ( $startPosiInExon < 1 );
                $endPosiInExon = $exonEnd - $absStart + 1;
                $endPosiInExon = $exonLength if ( $endPosiInExon > $exonLength );
            }

            my $relStart += $relExonStart + $startPosiInExon;
            my $relEnd += $relExonStart + $endPosiInExon;

            $absStartInExon = ( $absStart >= $exonStart ) ? $absStart : $exonStart; 
            $absEndInExon = ( $absEnd >= $exonEnd ) ? $exonEnd : $absEnd; 

            my $overlapString = join ( "\t", $absStartInExon, $absEndInExon, $transID, $relStart, $relEnd );
            push @{$ref_overlapRegion}, $overlapString;
        }

        $relExonStart += $exonLength;
    }

    1;
}

