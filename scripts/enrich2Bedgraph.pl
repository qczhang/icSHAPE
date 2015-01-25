#!/usr/bin/perl
#
use strict;
use warnings;

use lib "/home/qczhang/lib/perllib";
use Data::Dumper;
use Locale::Bio::IO qw( &readGTF &readFasta );

my %parameters = ( 
    gtfFile => "/home/qczhang/database/ensembl/current/mouse/gtf/Mus_musculus.GRCm38.74.gtf",
    gtfSource => "ensembl",
    fastaFile => "/home/qczhang/database/ensembl/current/mouse/gtf/transcriptome.fa"
);

my %readParameters = (
);

my $structureFile = shift;

my $ref_structure = &readStructureProfile ( $structureFile );
my $ref_annotation = readGTF ( $parameters{gtfFile}, 
    attribute => "transcript_id", 
    source => $parameters{gtfSource}, 
    verbose => 1, 
    skip => "Blat" ); smallfix ( $ref_annotation );
my ( $ref_seq ) = readFasta ( $parameters{fastaFile} );

foreach my $ensemblID ( keys %{$ref_structure} ) {
    next if ( $ensemblID =~ /rRNA/ );
    my @seq = split ( //, $ref_seq->{$ensemblID} );
    my $annoLenExon = &getExonLen ( $ref_annotation->{$ensemblID} );
    my $structLen = scalar ( @{$ref_structure->{$ensemblID}} );
    my $seqLen = scalar ( @seq );
    if ( ( $annoLenExon != $structLen ) or ( $structLen != $seqLen ) )  {
        print STDERR "Error! inconsistency within sequence, structure and/or annotation!\n";
        next;
    }
    
    my $ref_start = $ref_annotation->{$ensemblID}{exon}{start};
    my $ref_end = $ref_annotation->{$ensemblID}{exon}{end};
    my $exonNum = scalar( @{$ref_annotation->{$ensemblID}{exon}{start}} );

    ## maybe you also want to check whether exon strand is the same as overall strand, or you can trust annotations
    my $currentExon = 0;
    my $accumExonLen = 0;
    my $nextAccumExonLen = $ref_end->[$currentExon] - $ref_start->[$currentExon] + 1;
    for ( my $idxPos = 0; $idxPos <  scalar( @{$ref_structure->{$ensemblID}} ); $idxPos++ ) {
        if ( $idxPos >= $nextAccumExonLen ) {
            $currentExon++;
            $accumExonLen = $nextAccumExonLen;
            $nextAccumExonLen += $ref_end->[$currentExon] - $ref_start->[$currentExon] + 1;
        }

        if ( $ref_annotation->{$ensemblID}{strand} eq "+" ) {
            my $absPos = $idxPos - $accumExonLen + $ref_start->[$currentExon] - 1;
            print "chr", $ref_annotation->{$ensemblID}{exon}{seqName}[$currentExon], "\t", $absPos, "\t", $absPos+1; 
        }
        else {
            my $absPos = - $idxPos + $accumExonLen + $ref_end->[$currentExon] - 1;
            print "chr", $ref_annotation->{$ensemblID}{exon}{seqName}[$currentExon], "\t", $absPos, "\t", $absPos+1; 
        }
        print "\t", $ref_structure->{$ensemblID}[$idxPos], "\t", $ensemblID, "\t", $seq[$idxPos], "\n";
    }
}

## ---------------------
sub readStructureProfile 
{
    my $structureFile = shift;
    my %parameters = @_;

    my $ensemblIDcol = ( defined $parameters{ensemblIDcol} ) ? $parameters{ensemblIDcol} : 0;
    my $structureStartCol = ( defined $parameters{structureStartCol} ) ? $parameters{structureStartCol} : "3";

    open ( STR, $structureFile ) or die ( "Cannot open structure profile file $structureFile for reading!\n" );
    print STDERR "read in structure profile file $structureFile\n\t", `date`;
    my %trans_structure = ();  my $enrichScore = "NULL";  my $lineCount = 0;
    while ( my $line = <STR> ) {
        next if ( $line =~ /^#/ );
        $lineCount++;  print STDERR "  line $lineCount\n" if ( $lineCount % 1000 == 0 ); 

        chomp $line;
        my @data = split ( /\t/, $line );
        my $ensemblID = $data[$ensemblIDcol];
        $trans_structure{$ensemblID} = [ @data[$structureStartCol..$#data] ];
    }

    return \%trans_structure;
}

sub getExonLen 
{
    my $ref_transAnnotation = shift;

    ## iterate from the first exon on the 5'
    my $exonLen = 0;
    my @exonStart = ( sort {$a <=> $b} ( @{$ref_transAnnotation->{exon}{start}} ) );
    my @exonEnd = ( sort {$a <=> $b} ( @{$ref_transAnnotation->{exon}{end}} ) );

    for ( my $idxExon = 0; $idxExon < scalar ( @exonStart ); $idxExon++ ) {
        $exonLen += $exonEnd[$idxExon] - $exonStart[$idxExon] + 1;
    }

    return $exonLen;
}

sub get5primeLen 
{
    my $ref_transAnnotation = shift;

    ## iterate from the first exon on the 5'
    my $fivePrimeLen = 0;
    if ( defined $ref_transAnnotation->{start_codon} ) {
        if ( $ref_transAnnotation->{start_codon}{strand}[0] eq "+" ) {
            my $TSS = $ref_transAnnotation->{start_codon}{start}[0];
            my @exonStart = ( sort {$a <=> $b} ( @{$ref_transAnnotation->{exon}{start}} ) );
            my @exonEnd = ( sort {$a <=> $b} ( @{$ref_transAnnotation->{exon}{end}} ) );

            for ( my $idxExon = 0; $idxExon < scalar ( @exonStart ); $idxExon++ ) {
                if ( $exonEnd[$idxExon] >= $TSS ) {
                    $fivePrimeLen += $TSS - $exonStart[$idxExon];
                    last;
                }
                else {
                    $fivePrimeLen += $exonEnd[$idxExon] - $exonStart[$idxExon] + 1;
                }
            }
        }
        elsif ( $ref_transAnnotation->{start_codon}{strand}[0] eq "-" ) {
            my $TSS = $ref_transAnnotation->{start_codon}{end}[0];
            my @exonStart = ( sort {$b <=> $a} ( @{$ref_transAnnotation->{exon}{start}} ) );
            my @exonEnd = ( sort {$b <=> $a} ( @{$ref_transAnnotation->{exon}{end}} ) );

            for ( my $idxExon = 0; $idxExon < scalar ( @exonStart ); $idxExon++ ) {
                if ( $exonStart[$idxExon] <= $TSS ) {
                    $fivePrimeLen += $exonEnd[$idxExon] - $TSS;
                    last;
                }
                else {
                    $fivePrimeLen += $exonEnd[$idxExon] - $exonStart[$idxExon] + 1;
                }
            }
        }
    }

    return $fivePrimeLen;
}

sub get3primeLen 
{
    my $ref_transAnnotation = shift;

    ## iterate from the first exon on the 5'
    my $threePrimeLen = 0;
    if ( defined $ref_transAnnotation->{start_codon} ) {
        if ( $ref_transAnnotation->{stop_codon}{strand}[0] eq "-" ) {
            my $TES = $ref_transAnnotation->{stop_codon}{start}[0];
            my @exonStart = ( sort {$a <=> $b} ( @{$ref_transAnnotation->{exon}{start}} ) );
            my @exonEnd = ( sort {$a <=> $b} ( @{$ref_transAnnotation->{exon}{end}} ) );

            for ( my $idxExon = 0; $idxExon < scalar ( @exonStart ); $idxExon++ ) {
                if ( $exonEnd[$idxExon] >= $TES ) {
                    $threePrimeLen += $TES - $exonStart[$idxExon];
                    last;
                }
                else {
                    $threePrimeLen += $exonEnd[$idxExon] - $exonStart[$idxExon] + 1;
                }
            }
        }
        elsif ( $ref_transAnnotation->{stop_codon}{strand}[0] eq "+" ) {
            my $TES = $ref_transAnnotation->{stop_codon}{end}[0];
            my @exonStart = ( sort {$b <=> $a} ( @{$ref_transAnnotation->{exon}{start}} ) );
            my @exonEnd = ( sort {$b <=> $a} ( @{$ref_transAnnotation->{exon}{end}} ) );

            for ( my $idxExon = 0; $idxExon < scalar ( @exonStart ); $idxExon++ ) {
                if ( $exonStart[$idxExon] <= $TES ) {
                    $threePrimeLen += $exonEnd[$idxExon] - $TES;
                    last;
                }
                else {
                    $threePrimeLen += $exonEnd[$idxExon] - $exonStart[$idxExon] + 1;
                }
            }
        }
    }

    return $threePrimeLen;
}

sub smallfix
{
    my $ref_annotation = shift;
    $ref_annotation->{MM5_8SrRNA}{biotype} = [ "rRNA gene" ]; $ref_annotation->{MM18SrRNA}{biotype} = [ "rRNA gene" ]; $ref_annotation->{MM28SrRNA}{biotype} = [ "rRNA gene" ];
    $ref_annotation->{"lcl|gene=robert:MM18SrRNA"}{biotype} = [ "rRNA gene" ]; $ref_annotation->{"lcl|gene=robert:MM28SrRNA"}{biotype} = [ "rRNA gene" ];
}

