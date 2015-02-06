## --------------- package SHAPE-Seq ---------------
package icSHAPEutil;

use strict;
use warnings;
use Carp;
use File::Basename;

use base 'Exporter';
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK );

our $VERSION     = '0.01';
our @EXPORT      = ();
our @EXPORT_OK   = qw( readFasta extractFasta readBED readGTF readGenomeSize getExonLen get5primeLen get3primeLen getBioType );

sub readBED {
    my %parameters = @_;
}

sub readFasta {
    my $fastaFile = shift;

    my %fasta_len = ();
    my %fasta_sequence = ();
    my %fasta_alias = ();

    open ( FASTA, $fastaFile ) or die "Error in opening fasta file $fastaFile\n"; 
    print STDERR "read fasta sequence from $fastaFile\n";

    my $lineCount = 0;
    my $fastaID = "";
    my $fastaAlias = "";
    my $sequence = "";
    while ( my $line = <FASTA> ) {
        next if ( $line =~ /^#/ );

        chomp $line;
        if ( $line =~ /^>/ ) {
            if ( $fastaID ) {
                $fasta_len{$fastaID} = length ( $sequence );
                $fasta_sequence{$fastaID} = $sequence;
                $sequence = "";
                $fastaID = "";
            }

            $line =~ s/^>//;
            my @data = split ( /\s+/, $line );
            $fastaID = shift @data;
            foreach my $field ( @data ) {
                $fasta_alias{$field} = $fastaID;
            }
        }
        else {
            $sequence .= $line;
        }
    }
    close FASTA;

    if ( $fastaID ) {
        $fasta_len{$fastaID} = length ( $sequence );
        $fasta_sequence{$fastaID} = $sequence;
        $sequence = "";
        $fastaID = "";
    }

    return ( \%fasta_sequence, \%fasta_len, \%fasta_alias )
}

sub extractFasta 
{
    ## extract part of a sequence, closed interval, 1-indexed
    my $sequence = shift;  
    my ( $start, $end ) = @_;

    $start--;
    my $frag = "";
    if ( defined $sequence ) {
        my @wholeSequence = split ( //, $sequence );
        my $length = length ( $sequence );
        for ( my $idxPos = $start; $idxPos < $end; $idxPos++ ) {
            if ( ( $idxPos < 0 ) or ( $idxPos >= $length ) ) {
                $frag .= "-";
            }
            else {
                $frag .= $wholeSequence[$idxPos];
            }
        }
    }

    return $frag;
}

sub readGTF_ensembl {
    my $gtfFile = shift;
    my %parameters = @_;

    print STDERR "read genomic annotations from file: $gtfFile\n\t", `date` if ( $parameters{verbose} );
    my %annotation = ();
    my $lineCount = 0;
    open ( GTF, $gtfFile );
    while ( my $line = <GTF> ) {
        next if ( $line =~ /^#/ );
        chomp $line;
        $lineCount++;
        print STDERR "  line: $lineCount\n" if ( ( $parameters{verbose} ) and ( $lineCount % 100000 == 0 ) );
        my ( $seqName, $source, $feature, $start, $end, $score, $strand, $frame, $attribute, @fields ) = split ( /\t/, $line );
        #1       nonsense_mediated_decay stop_codon      4782680 4782682 .       -       0       gene_id "ENSMUSG00000033845"; transcript_id "ENSMUST00000045689"; exon_number "2"; gene_name "Mrpl15"; gene_biotype "protein_coding"; transcript_name "Mrpl15-003";
        my $location = $seqName . ":" . $start . "-" . $end;
        $attribute =~ s/;$//;

        if ( $parameters{attribute} ) {
            ## read GTF file in attribute centric way
            my @attr_values = split ( /; /, $attribute );
            my $foundAttribute = 0;
            foreach my $attr_value ( @attr_values ) {
                my ( $attr, $value ) = ( $attr_value =~ /(.+) "(.+)"$/ );
                if ( $attr eq $parameters{attribute} ) {
                    $foundAttribute = $value;
                    last;
                }
            }

            if ( $foundAttribute ) {
                push ( @{$annotation{$foundAttribute}{$feature}{strand}}, $strand );
                push ( @{$annotation{$foundAttribute}{$feature}{location}}, $location );
                push ( @{$annotation{$foundAttribute}{$feature}{seqName}}, $seqName );
                push ( @{$annotation{$foundAttribute}{$feature}{start}}, $start );
                push ( @{$annotation{$foundAttribute}{$feature}{end}}, $end );

                foreach my $attr_value ( @attr_values ) {
                    my ( $attr, $value ) = ( $attr_value =~ /(.+) "(.+)"$/ );
                    if ( $attr and ( $attr ne $parameters{attribute} ) ) {
                        if ( not defined $annotation{$foundAttribute}{$attr} ) {
                            $annotation{$foundAttribute}{$attr} = $value;
                        }
                        elsif ( $annotation{$foundAttribute}{$attr} !~ /$value/ ) {
                            $annotation{$foundAttribute}{$attr} .= ";" . $value;
                        }
                    }
                }
                if ( not defined $annotation{$foundAttribute}{strand} ) {
                    $annotation{$foundAttribute}{strand} = $strand;
                }
                elsif ( $annotation{$foundAttribute}{strand} !~ /\Q$strand\E/ ) {
                    $annotation{$foundAttribute}{strand} .= ";" . $strand;
                }
            }
        }
        elsif ( $parameters{location} ) {
            ## read GTF file in location centric way
            push @{$annotation{$feature}{$seqName}{source}}, $source;
            push ( @{$annotation{$feature}{$seqName}{start}}, $start );
            push ( @{$annotation{$feature}{$seqName}{end}}, $end );
            push ( @{$annotation{$feature}{$seqName}{strand}}, $strand );

            if ( $parameters{moreInfo} ) {
                push ( @{$annotation{$feature}{$seqName}{score}}, $score );

                my @attr_values = split ( /;/, $attribute );
                foreach my $attr_value ( @attr_values ) {
                    my ( $attr, $value ) = split ( /\s+/, $attr_value );
                    push @{$annotation{$feature}{$seqName}{$attr}}, $value;
                }
            }
        }
    }
    close GTF;

    return \%annotation;
}

sub readGTF {
    my $gtfFile = shift;
    my %parameters = @_;

    if ( $parameters{source} =~ /ensembl/i ) {
        return readGTF_ensembl ( $gtfFile, %parameters );
    }

    my %skipSource = ();
    if ( defined $parameters{skip} ) {
        %skipSource = map { $_ => 1 } ( split ( /;/, $parameters{skip} ) );
    }

    print STDERR "read genomic annotations from file: $gtfFile\n\t", `date` if ( $parameters{verbose} );
    my %annotation = ();
    my $lineCount = 0;
    open ( GTF, $gtfFile );
    while ( my $line = <GTF> ) {
        next if ( $line =~ /^#/ );
        chomp $line;
        $lineCount++;
        print STDERR "  line: $lineCount\n" if ( ( $parameters{verbose} ) and ( $lineCount % 100000 == 0 ) );
        my ( $seqName, $source, $feature, $start, $end, $score, $strand, $frame, $attribute, @fields ) = split ( /\t/, $line );
        #  3       ENSEMBL start_codon     24333068        24333070        .       +       0       Parent=ENSEMBL:ENSMUST00000077389;Dbxref=MGI:MGI:3645137,ENSEMBL:ENSMUSG00000057036
        next if ( defined $skipSource{$source} );
        my $location = $seqName . ":" . $start . "-" . $end;

        if ( $parameters{attribute} ) {
            ## read GTF file in attribute centric way
            if ( ( $parameters{source} =~ /ensembl/i ) or ( $parameters{source} =~ /mgi/i ) ) {
                ## ensembl database example:
                #  gene_id "ENSMUSG00000090025"; transcript_id "ENSMUST00000160944"; exon_number "1"; gene_name "Gm16088"; gene_biotype "pseudogene"; transcript_name "Gm16088-001"; exon_id "ENSMUSE00000848981";
                $attribute =~ s/;$//;
                my @attr_values = ();
                if ( $parameters{source} =~ /ensembl/i ) { 
                    @attr_values = split ( /; /, $attribute );
                }
                elsif ( $parameters{source} =~ /mgi/i )  {
                    @attr_values = split ( /;/, $attribute );
                }
                my $foundAttribute = 0;
                my $foundParent = 0;
                foreach my $attr_value ( @attr_values ) {
                    my ( $attr, $value ) = ( "", "" );
                    if ( $parameters{source} =~ /ensembl/i ) { 
                        ( $attr, $value ) = ( $attr_value =~ /(.+) "(.+)"$/ );
                    }
                    elsif ( $parameters{source} =~ /mgi/i )  {
                        ( $attr, $value ) = ( $attr_value =~ /(.+)=(.+)$/ );
                        $attr = lc ( $attr );
                    }

                    if ( $attr eq $parameters{attribute} ) {
                        $foundAttribute = $value;
                    }
                    if ( $attr eq "parent" ) {
                        $foundParent = $value;
                    }
                }
                if ( $foundAttribute ) {
                    push ( @{$annotation{$foundAttribute}{source}}, $source );
                    push ( @{$annotation{$foundAttribute}{feature}}, $feature );
                    push ( @{$annotation{$foundAttribute}{strand}}, $strand );
                    push ( @{$annotation{$foundAttribute}{location}}, $location );
                    push ( @{$annotation{$foundAttribute}{seqName}}, $seqName );
                    push ( @{$annotation{$foundAttribute}{start}}, $start );
                    push ( @{$annotation{$foundAttribute}{end}}, $end );

                    foreach my $attr_value ( @attr_values ) {
                        my ( $attr, $value ) = ( "", "" );
                        if ( $parameters{source} =~ /ensembl/i ) { 
                            ( $attr, $value ) = ( $attr_value =~ /(.+) "(.+)"$/ );
                        }
                        elsif ( $parameters{source} =~ /mgi/i )  {
                            ( $attr, $value ) = ( $attr_value =~ /(.+)=(.+)$/ );
                            $attr = lc ( $attr );
                        }
                        if ( $attr and ( $attr ne $parameters{attribute} ) ) {
                            push @{$annotation{$foundAttribute}{$attr}}, $value;
                        }
                    }
                }
                elsif ( $foundParent ) {
                    push ( @{$annotation{$foundParent}{$feature}{strand}}, $strand );
                    push ( @{$annotation{$foundParent}{$feature}{location}}, $location );
                    push ( @{$annotation{$foundParent}{$feature}{seqName}}, $seqName );
                    push ( @{$annotation{$foundParent}{$feature}{start}}, $start );
                    push ( @{$annotation{$foundParent}{$feature}{end}}, $end );

                    foreach my $attr_value ( @attr_values ) {
                        my ( $attr, $value ) = ( "", "" );
                        if ( $parameters{source} =~ /ensembl/i ) { 
                            ( $attr, $value ) = ( $attr_value =~ /(.+) "(.+)"$/ );
                        }
                        elsif ( $parameters{source} =~ /mgi/i )  {
                            ( $attr, $value ) = ( $attr_value =~ /(.+)=(.+)$/ );
                            $attr = lc ( $attr );
                        }
                        if ( $attr and ( $attr ne $parameters{attribute} ) ) {
                            next if ( $attr eq "parent" );
                            push @{$annotation{$foundParent}{$feature}{$attr}}, $value;
                        }
                    }
                }
            }
        }
        elsif ( $parameters{feature} ) {

        }
        else {
            ## read GTF file in location centric way
            push @{$annotation{$feature}{$seqName}{source}}, $source;
            push ( @{$annotation{$feature}{$seqName}{start}}, $start );
            push ( @{$annotation{$feature}{$seqName}{end}}, $end );
            push ( @{$annotation{$feature}{$seqName}{strand}}, $strand );

            if ( $parameters{moreInfo} ) {
                push ( @{$annotation{$feature}{$seqName}{score}}, $score );

                if ( $parameters{source} =~ /ensembl/i ) {
                    my @attr_values = split ( /;/, $attribute );
                    foreach my $attr_value ( @attr_values ) {
                        my ( $attr, $value ) = split ( /\s+/, $attr_value );
                        push @{$annotation{$feature}{$seqName}{$attr}}, $value;
                    }
                }
                elsif ( ( $parameters{source} =~ /flybase/i ) or ( $parameters{source} =~ /mgi/i ) ) {
                    my $id = "NULL";
                    if ( ( $attribute =~ /ID=(\S+?);/ ) or ( $attribute =~ /ID=(\S+?)$/ ) ) {
                        $id = $1;
                    }
                    push ( @{$annotation{$feature}{$seqName}{ID}}, $id );
                    my $name = "NULL";
                    if ( ( $attribute =~ /Name=(\S+?);/ ) or ( $attribute =~ /Name=(\S+?)$/ ) ) {
                        $name = $1;
                    }
                    push ( @{$annotation{$feature}{$seqName}{Name}}, $name );
                }
            }
        }
    }
    close GTF;

    return \%annotation;
}

sub readGenomeSize {
    my $genomeSizeFile = shift;

    open ( IN, $genomeSizeFile ) || die "Error in opening file $genomeSizeFile! Exit.\n";
    my %target_size = ();
    while ( my $line = <IN> ) {
        next if ( $line =~ /^#/ );

        chomp $line;
        my ( $target, $size ) = split ( /\t/, $line );
        $target_size{$target} = $size;
    }
    close IN;

    return \%target_size;
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

    return $fivePrimeLen;
}

sub get3primeLen 
{
    my $ref_transAnnotation = shift;

    ## iterate from the first exon on the 5'
    my $threePrimeLen = 0;
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

    return $threePrimeLen;
}

sub getBioType 
{
    my $transID = shift;
    my $ref_annotation = shift;

    my $biotypeString = "NULL";
    if ( defined $ref_annotation->{$transID} ) {
        if ( defined $ref_annotation->{$transID}{gene_biotype} ) {
            $biotypeString = $ref_annotation->{$transID}{gene_biotype};
        }

        $biotypeString =~ s/^_//; $biotypeString =~ s/_/ /g;
        $biotypeString =~ s/PSEUDO/pseudogene/;
        $biotypeString =~ s/^gene$/unclassified_gene/i; $biotypeString =~ s/unclassified gene$/unclassified_gene/i; $biotypeString =~ s/ gene$//;
    }

    return $biotypeString;
}



1;
