#! /usr/bin/perl
#
use strict;
use warnings;

use Getopt::Std;
use lib "$ENV{ICSHAPE}/module";
use IO qw( &readGTF &readFasta );

my $scoreFile = shift;
my $statFile = shift;
my $gtfFile = shift;
my $fastaFile = shift;

open ( SC, $scoreFile );
open ( STAT, ">$statFile" );

my %parameters = ( 
    gtfFile => $gtfFile,
    gtfSource => "ensembl"
);
my $ref_annotation = readGTF ( $parameters{gtfFile},
    attribute => "transcript_id",
    source => $parameters{gtfSource},
    verbose => 1,
    skip => "Blat" ); smallfix ( $ref_annotation );
my ( $ref_ensembl_seq ) = readFasta ( $fastaFile );

my $lineCount = 0;
my %allScore = ();
while ( my $line = <SC> ) {
    next if ( $line =~ /^#/ );
    chomp $line;

    my $total = 0; my $count = 0;
    my ( $ensemblID, $len, $rpkm, @score ) = split ( /\t/, $line );
    my $transSeq = uc ( $ref_ensembl_seq->{$ensemblID} );

    my $bioType = &getBioType ( $ref_annotation->{$ensemblID} );
    my $lenExon = 0; my $len5prime = 0; my $len3prime = 0;
    if ( $bioType eq "protein coding" ) {
        $len5prime = &get5primeLen ( $ref_annotation->{$ensemblID} );
        $len3prime = &get3primeLen ( $ref_annotation->{$ensemblID} );
        $lenExon = &getExonLen ( $ref_annotation->{$ensemblID} );
    }

    my $gini = &calcGini ( \@score, format => 1 );
    print STAT $ensemblID, "\t", $len, "\t", $rpkm, "\t", $gini, "\t", $bioType;
    if ( $bioType eq "protein coding" ) {
        my $gini5prime = &calcGini ( \@score, start => 0, end => $len5prime-1, format => 1 );
        my $ginicds = &calcGini ( \@score, start => $len5prime, end => $lenExon-$len3prime-1, format => 1 );
        my $gini3prime = &calcGini ( \@score, start => $lenExon-$len3prime, end => $lenExon-1, format => 1 );
    
        print STAT "\t", $gini5prime, "\t", $ginicds, "\t", $gini3prime;
    }
    print STAT "\n";

    my $sampleCount = int ( scalar ( @score ) / 100 );
    for ( my $idx = 0; $idx < $sampleCount; $idx++ ) {
        my $rand = int ( rand ( scalar ( @score ) ) );
        if ( $score[$rand] ne "NULL" ) {
            push @{$allScore{all}}, $score[$rand];

            my $base = substr ( $transSeq, $rand, 1 );
            push @{$allScore{$base}}, $score[$rand];
        }
    }

    $lineCount++;
    #last if ( $lineCount > 100 );
}

foreach my $class ( keys %allScore ) {
   my $gini = &calcGini ( $allScore{$class}, format => 1 );
   print STAT "overall\t", $class, "\t", $gini, "\n";
}

close SC;
close STAT;


## ------------------
sub calcGini 
{
    my $ref_array = shift;
    my %parameters = @_;

    my $len = scalar ( @{$ref_array} );
    my $start = defined ( $parameters{start} ) ? $parameters{start} : 0;
    my $end = defined ( $parameters{end} ) ? $parameters{end} : $len - 1;

    my $gini = "NULL";
    my @num = (); 
    for ( my $idx = $start; $idx < $end; $idx++ ) { 
        if ( $ref_array->[$idx] ne "NULL" ) {
            push @num, $ref_array->[$idx];
        }
    }
    @num = sort { $a <=> $b } @num;

    my $count = scalar ( @num );
    if ( $count > 10 ) {
        my $total = eval ( join ( "+", @num ) );
        if ( $total > 0 )  {
            my $giniB = 0;  my $accum = 0;
            for ( my $idx = 0; $idx < scalar(@num); $idx++ ) { 
                $accum += $num[$idx];
                #$giniB += ( ( $accum - $num[$idx]/2 ) / $total / $count ) ; 
                $giniB += ( $accum / $total / $count ) ; 
            }
            if ( $parameters{format} ) {
                $gini = sprintf ( "%.3f", 1 - 2*$giniB );
            }
            else {
                $gini = 1 - 2*$giniB;
            }
        }
    }

    return $gini;
}

sub getBioType 
{
    my $ref_annotation = shift;

    my $biotypeString = "not found";
    if ( defined $ref_annotation ) {
        if ( defined $ref_annotation->{gene_biotype} ) {
            $biotypeString = $ref_annotation->{gene_biotype};
        }

        $biotypeString =~ s/^_//; $biotypeString =~ s/_/ /g;
        $biotypeString =~ s/PSEUDO/pseudogene/;
        $biotypeString =~ s/^gene$/unclassified_gene/i; $biotypeString =~ s/unclassified gene$/unclassified_gene/i; $biotypeString =~ s/ gene$//;
    }

    return $biotypeString;
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
    if ( defined $ref_transAnnotation->{stop_codon} ) {
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
    $ref_annotation->{MM5_8SrRNA}{gene_biotype} = "rRNA gene"; $ref_annotation->{MM18SrRNA}{gene_biotype} = "rRNA gene"; $ref_annotation->{MM28SrRNA}{gene_biotype} = "rRNA gene";
    $ref_annotation->{"lcl|gene=robert:MM18SrRNA"}{gene_biotype} = "rRNA gene"; $ref_annotation->{"lcl|gene=robert:MM28SrRNA"}{gene_biotype} = "rRNA gene";
}
