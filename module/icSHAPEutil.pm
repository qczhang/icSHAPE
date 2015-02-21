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
our @EXPORT_OK   = qw( readFasta extractFasta readBED readGTF readGTF_ensembl_new readGenomeSize getExonLen get5primeLen get3primeLen getBioType );

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

sub readGTF_ensembl_new {
    my $gtfFile = shift;
    my %parameters = @_;

    my %chr_size = (); my %gene_info = (); my %transcript_info = (); my %exon_info = ();

    my $lineCount = 0;
    open ( GTF, $gtfFile ) or die ( "Error in reading GTF file $gtfFile!\n" );
    print STDERR "read genomic annotations from file: $gtfFile\n\t", `date` if ( $parameters{verbose} );
    my $geneID = "";  my $transcriptID = "";  my $exonID = "";  my $exonNum = 0;
    while ( my $line = <GTF> ) {
        next if ( $line =~ /^#/ ); chomp $line;
        $lineCount++; print STDERR "  line: $lineCount\n" if ( ( $parameters{verbose} ) and ( $lineCount % 100000 == 0 ) );

        # 1     nonsense_mediated_decay     stop_codon      4782680     4782682     .       -       0       
        # gene_id "ENSMUSG00000033845"; transcript_id "ENSMUST00000045689"; exon_number "2"; 
        # gene_name "Mrpl15"; gene_biotype "protein_coding"; transcript_name "Mrpl15-003";
        # exon_id ...
        ## this is the statistics of feature field of mouse.74.gtf (in total 1242614 lines):
        #   432247 CDS
        #   628052 exon
        #    41170 start_codon
        #    41145 stop_codon
        ## this is the statistics of attribute field of mouse.74.gtf (in total 1242614 lines):
        #  1142614 gene_id
        #  1142614 gene_name
        #  1142614 gene_biotype
        #  1142614 transcript_id
        #  1142614 transcript_name
        #  1142614 exon_number
        #   628052 exon_id
        #   432247 protein_id
        my ( $chr, $class, $feature, $start, $end, $score, $strand, $frame, $attribute ) = split ( /\t/, $line );
        my @data = split ( /; /, $attribute );
        foreach my $field ( @data ) {
            my $index = index ( $field, " " );
            next if ( $index <= 0 );
            my $type = substr ( $field, 0, $index ); 
            if ( $type eq "gene_id" ) { $geneID = substr ( $field, $index+2, -1 ); }
            elsif ( $type eq "transcript_id" ) { $transcriptID = substr ( $field, $index+2, -1 ); }
            elsif ( $type eq "exon_id" ) { $exonID = substr ( $field, $index+2, -1 ); }
            elsif ( $type eq "exon_number" ) { $exonNum = substr ( $field, $index+2, -1 ); }
        }
=cut
        if ( ( $feature eq "gene" ) and ( $geneID ) ) {                   # not defined in ensembl
            if ( defined $gene_info{$geneID} ) { print STDERR "Warnning! skipping line $lineCount of repeated geneID: $geneID\n\t$line\n"; next; }
            $gene_info{$geneID}{chr} = $chr; $gene_info{$geneID}{strand} = $strand;
            $gene_info{$geneID}{start} = $start; $gene_info{$geneID}{end} = $end;
        }
        if ( ( $feature eq "transcript" ) and ( $transcriptID ) ) {       # not defined in ensembl
            if ( not $geneID ) { print STDERR "Warnning! skipping line $lineCount of no geneID:\n\t$line\n"; next; }
            if ( defined $transcript_info{$transcriptID} ) 
                { print STDERR "Warnning! skipping line $lineCount of repeated transcriptID: $transcriptID\n\t$line\n"; next; }
            push @{$gene_info{$geneID}{transcript}}, $transcriptID;
            $transcript_info{$transcriptID}{gene} = $geneID;
            $transcript_info{$transcriptID}{start} = $start; $transcript_info{$transcriptID}{end} = $end;
        }
=cut
        if ( $feature eq "exon" ) {                   # so far we only deal with exon
            if ( not $geneID ) { print STDERR "Warnning! skipping line $lineCount of no geneID:\n\t$line\n"; next; }
            if ( not $transcriptID ) { print STDERR "Warnning! skipping line $lineCount of no transcriptID:\n\t$line\n"; next; }
            if ( not $exonID ) { print STDERR "Warnning! skipping line $lineCount of no exonID:\n\t$line\n"; next; }

            if ( defined $exon_info{$exonID} ) {
                if ( ( $exon_info{$exonID}{start} != $start ) or ( $exon_info{$exonID}{end} != $end ) ) 
                    { print STDERR "Error! line $lineCount of inconsistent exonID $exonID:\n\t$line\n"; next; }
            }
            if ( defined $transcript_info{$transcriptID} ) {
                if ( not defined $gene_info{$geneID} ) 
                    { print STDERR "Warnning! in consistent transcript annotation of transcript $transcriptID in $line\n"; next; }
                if ( $transcript_info{$transcriptID}{gene} ne $geneID ) 
                    { print STDERR "Warnning! in consistent transcript annotation of transcript $transcriptID in $line\n"; next; }
                if ( defined $transcript_info{$transcriptID}{exon}{$exonNum} ) 
                    { print STDERR "Warnning! skipping line $lineCount of repeated exon_number $exonNum:\n\t$line\n"; next; }
            }
            if ( defined $gene_info{$geneID} ) {
                if ( ( $gene_info{$geneID}{chr} ne $chr ) or ( $gene_info{$geneID}{strand} ne $strand ) )
                    { print STDERR "Warnning! in consistent location annotation of gene $geneID in $line\n"; next; }
            }

            if ( not defined $chr_size{$chr} ) { $chr_size{$chr} = $end; }
            else { $chr_size{$chr} = $end if ( $end > $chr_size{$chr} ); }

            if ( not defined $gene_info{$geneID} ) {
                $gene_info{$geneID}{chr} = $chr; $gene_info{$geneID}{strand} = $strand;
                $gene_info{$geneID}{start} = $start; $gene_info{$geneID}{end} = $end;
            }
            else {
                $gene_info{$geneID}{start} = $start if ( $start < $gene_info{$geneID}{start} ); 
                $gene_info{$geneID}{end} = $end if ( $end > $gene_info{$geneID}{end} );
            }

            if ( not defined $transcript_info{$transcriptID} ) {
                $transcript_info{$transcriptID}{gene} = $geneID;
                $transcript_info{$transcriptID}{start} = $start; $transcript_info{$transcriptID}{end} = $end;
                push @{$gene_info{$geneID}{transcript}}, $transcriptID;
            }
            else {
                $transcript_info{$transcriptID}{start} = $start if ( $start < $transcript_info{$transcriptID}{start} ); 
                $transcript_info{$transcriptID}{end} = $end if ( $end > $transcript_info{$transcriptID}{end} );
            }

            $transcript_info{$transcriptID}{exon}{$exonNum} = $exonID;
            $exon_info{$exonID}{start} = $start; $exon_info{$exonID}{end} = $end;
        }
    }
    close GTF;

    return { 
        chr_size            => \%chr_size,
        gene_info           => \%gene_info, 
        transcript_info     => \%transcript_info, 
        exon_info           => \%exon_info 
    };
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
