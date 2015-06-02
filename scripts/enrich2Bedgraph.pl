#!/usr/bin/perl
#
use strict;
use warnings;

use Getopt::Std;
use lib "$ENV{ICSHAPE}/module";
use IO qw( &readGTF &readFasta );

use vars qw ($opt_h $opt_V $opt_D $opt_i $opt_o $opt_a $opt_g $opt_f $opt_c $opt_s );
&getopts('hVDi:o:a:g:f:c:s:');

my $usage = <<_EOH_;
## --------------------------------------
Convert transcript-wise file of enrichment structural scores to genome-wise bedgraph file

Command:
$0 -i transcript_profile -o output_bedgraph_file

# what it is:
 -i     profile for every transcript, in the format of:
        transcript_ID  [more info] position1_data  position2_data  ...
 -o     output bedgraph file

 -g     gtf annotation file (annotations from ensembl are suggested)
 -f     gtf file format (default: ensembl)
 -a     fasta file

 -s     genome size file - only these chromosomes will be kept

_EOH_
;

&main ();

sub main
{
    my %parameters = &init();

    my $ref_structure = &readStructureProfile ( $parameters{transcriptFile} );
    my $ref_annotation = readGTF ( $parameters{gtfFile}, 
        attribute => "transcript_id", 
        source => $parameters{gtfSource}, 
        verbose => 1, 
        skip => "Blat" ); smallfix ( $ref_annotation );
    my ( $ref_seq ) = readFasta ( $parameters{fastaFile} );
    my $ref_chr_size = readChrSize ( $parameters{genomeSize} ) if ( defined $parameters{genomeSize} );

    if ( defined $parameters{bedgraphFile} ) { open ( OUT, ">$parameters{bedgraphFile}" ); }
    foreach my $ensemblID ( keys %{$ref_structure} ) {
        next if ( $ensemblID =~ /rRNA/ );
        my $ensemblChr = $ref_annotation->{$ensemblID}{exon}{seqName}[0];
        next if ( not defined $ref_chr_size->{$ensemblChr} );

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
                if ( defined $parameters{bedgraphFile} ) { print OUT "chr", $ref_annotation->{$ensemblID}{exon}{seqName}[$currentExon], "\t", $absPos, "\t", $absPos+1; }
                else { print "chr", $ref_annotation->{$ensemblID}{exon}{seqName}[$currentExon], "\t", $absPos, "\t", $absPos+1; }
            }
            else {
                my $absPos = - $idxPos + $accumExonLen + $ref_end->[$currentExon] - 1;
                if ( defined $parameters{bedgraphFile} ) { print OUT "chr", $ref_annotation->{$ensemblID}{exon}{seqName}[$currentExon], "\t", $absPos, "\t", $absPos+1; }
                else { print "chr", $ref_annotation->{$ensemblID}{exon}{seqName}[$currentExon], "\t", $absPos, "\t", $absPos+1; }
            }
                if ( defined $parameters{bedgraphFile} ) { print OUT "\t", $ref_structure->{$ensemblID}[$idxPos], "\t", $ensemblID, "\t", $seq[$idxPos], "\n"; }
                else { print "\t", $ref_structure->{$ensemblID}[$idxPos], "\t", $ensemblID, "\t", $seq[$idxPos], "\n"; }
        }
    }
    if ( defined $parameters{bedgraphFile} ) { close OUT; }
}

## ---------------------
sub init
{
    my %parameters = ();
    die $usage if ( $opt_h || ( not $opt_i ) );

    $opt_V = 0 if ( not defined $opt_V );
    $opt_D = 0 if ( not defined $opt_D );
    my $pwd = `pwd`;  chomp $pwd;

    $parameters{transcriptFile} = $opt_i;
    if ( defined $opt_o ) { $parameters{bedgraphFile} = $opt_o; }

    if ( defined $opt_g ) { $parameters{gtfFile} = $opt_g; }
    else {
        print STDERR "Please specify genome annotation GTF file!\n";
        die $usage;
    }
    if ( defined $opt_f ) { $parameters{gtfSource} = $opt_f; }
    else { $parameters{gtfSource} = "ensembl"; }
    if ( defined $opt_a ) { $parameters{fastaFile} = $opt_a; }
    else {
        print STDERR "Please specify transcript fasta file!\n";
        die $usage;
    }
    if ( defined $opt_s ) { $parameters{genomeSize} = $opt_s; }

    return ( %parameters );
}


    my $ref_chr_size = readChrSize ( $parameters{genomeSize} ) if ( defined $parameters{genomeSize} );

sub readChrSize
{
    my $genomeSizeFile = shift;
    my %chr_size = ();

    open ( GS, $genomeSizeFile );
    while ( my $line = <GS> ) {
        chomp $line;
        my ( $chr, $size ) = split ( /\t/, $line );
        $chr =~ s/^[chr|Chr|CHR]//;

        $chr_size{$chr} = $size;
    }
    close GS;

    return \%chr_size;
}

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

