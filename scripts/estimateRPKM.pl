#! /usr/bin/perl
#
#
use strict;
use warnings;
use File::Basename;
use Getopt::Std;
use Scalar::Util qw(looks_like_number);

use vars qw ($opt_h $opt_V $opt_D $opt_i $opt_o $opt_r $opt_m $opt_c );
&getopts('hVDi:o:r:m:c:');

my $usage = <<_EOH_;
## --------------------------------------
calculate expression values for all transcripts in FPKM from SAM files 

Command:
$0 -i input_sam_file -o output_expression_file

# what it is:
 -i     input sam files. If of paired ends, need to be sorted according to paired reads
 -o     output expression files in RPKM

# more options:
 -r     reference transcript list (only transcripts in this list are considered)
 -c     cutoff to output (defaut 0.1)
 -m     method ( not yet )

_EOH_
    ;

&main();

sub main {
    my %parameters = &init();

    my %transcript_len = ();
    my %transcript_rpkm = ();
    my %transcript_uniqRead = (); my %transcript_multiRead = ();

    my @samFiles = split ( /:/, $parameters{input} );
    my $firstSam = 1;   my $totalReadNum = 0;  my $totalUniqReadNum = 0;  my $totalMultiReadNum = 0;
    if ( $parameters{method} eq "even" ) {
        foreach my $samFile ( @samFiles ) {
            my ( $uniqReadNum, $multiReadNum, $readNum ) = &calcTranscriptAbundance_evenDistribute ( $samFile, \%transcript_len, \%transcript_uniqRead, \%transcript_multiRead, $firstSam );
            $totalUniqReadNum += $uniqReadNum;
            $totalMultiReadNum += $multiReadNum;
            $totalReadNum += $readNum;

            $firstSam = 0;
        }
        print STDERR "finally $totalReadNum hits of $totalUniqReadNum uniq hit reads and $totalMultiReadNum multi hit reads read in.\n";

        &calcRPKM ( \%transcript_rpkm, \%transcript_len, $totalReadNum, \%transcript_uniqRead, \%transcript_multiRead );
    }
    elsif ( $parameters{method} eq "mle" ) {
    }

    &output_rpkm ( $parameters{output}, $parameters{cutoff},  \%transcript_len, \%transcript_uniqRead, \%transcript_multiRead, \%transcript_rpkm );
    print STDERR "\nFinished calculating transcript abundance, check $parameters{output} for them in RPKM.\n"; 
}


sub init {
    my %parameters = ();

    die $usage if ( $opt_h || ( not defined $opt_i ) || ( not defined $opt_o ) );
    $opt_V = 0 if ( not defined $opt_V );
    $opt_D = 0 if ( not defined $opt_D );

    $parameters{input} = $opt_i;
    $parameters{output} = $opt_o;
    $parameters{method} = "even";
    $parameters{cutoff} = ( defined $opt_c ) ? $opt_c : 0.1 ;

    return ( %parameters );
}

sub calcTranscriptAbundance_evenDistribute {
    my ( $inputSam, $ref_transcript_len, $ref_transcript_uniqRead, $ref_transcript_multiRead, $firstSam )  = @_;
    print STDERR "Estimation of transcript expression level by equally distribute multiple hits...\n\t", `date`;

    open ( IN, $inputSam ) || die "Error in opening file $inputSam! Exit.\n";
    print STDERR "Read in SAM file $inputSam...\n\t", `date`;
    my $readOld = "";  my @hitIDs = ();  my $hitCount = 0;  
    my $totalCount = 0; my $multiHitCount = 0;  my $uniqHitCount = 0;
    while ( my $line = <IN> ) {
        if ( $line =~ /^@/ ) {
            if ( $firstSam ) {
                my ( $type, $seqID, $len ) = split ( /\t/, $line );
                if ( $type =~ /SQ/ ) {
                    $seqID =~ s/SN://;
                    $len =~ s/LN://; chomp $len;
                    $ref_transcript_len->{$seqID} = $len;
                    $ref_transcript_uniqRead->{$seqID} = 0;
                    $ref_transcript_multiRead->{$seqID} = 0;
                }
            }
        }
        else {
            my ( $read, $tag, $hit, $pos, $mapq, $cigar, $rnext, $pnext, $tlen, $seq, $qual ) = split ( /\t/, $line );

            next if ( ( not looks_like_number($tag) ) or ( not looks_like_number($pos) ) );  
            if ( ( $tag == 99 ) or ( $tag == 355 ) ) { next if ( $rnext ne "=" ); } ## PE
            elsif ( ( $tag == 0 ) or ( $tag == 256 ) ) { } ## SE
            else { next; }

            $totalCount++; print STDERR $totalCount, "\n\t", `date` if ( $totalCount % 1000000 == 0 );
            if ( $read ne $readOld ) {
                if ( $readOld ) {
                    if ( $hitCount > 1 ) {
                        foreach my $hit ( @hitIDs ) { 
                            next if ( not defined $ref_transcript_uniqRead->{$hit} );
                            $ref_transcript_multiRead->{$hit} += ( 1/$hitCount ); 
                        }
                        $multiHitCount++;
                    }
                    elsif ( $hitCount == 1 ) {
                        next if ( not defined $ref_transcript_uniqRead->{$hitIDs[0]} );
                        $ref_transcript_uniqRead->{$hitIDs[0]}++;
                        $uniqHitCount++;
                    }
                }
                $readOld = $read; @hitIDs = (); $hitCount = 0;
            }
            push @hitIDs, $hit; $hitCount++; 
        }
    }
    close IN;

    if ( $readOld ) {
        if ( $hitCount > 1 ) {
            foreach my $hit ( @hitIDs ) { 
                next if ( not defined $ref_transcript_uniqRead->{$hit} );
                $ref_transcript_multiRead->{$hit} += ( 1/$hitCount ); 
            }
            $multiHitCount++;
        }
        elsif ( $hitCount == 1 ) {
            next if ( not defined $ref_transcript_uniqRead->{$hitIDs[0]} );
            $ref_transcript_uniqRead->{$hitIDs[0]}++;
            $uniqHitCount++;
        }
    }
    
    print STDERR "\ttotal uniq reads $uniqHitCount\n\t", `date`;
    print STDERR "\ttotal multiple reads $multiHitCount\n\t", `date`;
    print STDERR "\ttotal reads $totalCount\n\t", `date`;

    return ( $uniqHitCount, $multiHitCount, $totalCount );
}

sub calcRPKM  {
    my ( $ref_transcript_rpkm, $ref_transcript_len, $totalRead, $ref_transcript_uniqRead, $ref_transcript_multiRead ) = @_;

    if ( defined $ref_transcript_multiRead ) {
        foreach my $tran ( keys %{$ref_transcript_len} ) {
            if ( defined $ref_transcript_uniqRead->{$tran} ) {
                if ( defined $ref_transcript_multiRead->{$tran} ) {
                    $ref_transcript_rpkm->{$tran} = ( ( $ref_transcript_uniqRead->{$tran} + $ref_transcript_multiRead->{$tran} ) * 1000000000 ) / ( $ref_transcript_len->{$tran} * $totalRead );
                }
                else {
                    $ref_transcript_rpkm->{$tran} = ( ( $ref_transcript_uniqRead->{$tran} ) * 1000000000 ) / ( $ref_transcript_len->{$tran} * $totalRead );
                }
            }
            else {
                if ( defined $ref_transcript_multiRead->{$tran} ) {
                    $ref_transcript_rpkm->{$tran} = ( ( $ref_transcript_multiRead->{$tran} ) * 1000000000 ) / ( $ref_transcript_len->{$tran} * $totalRead );
                }
            }
        }
    }
    else {
        foreach my $tran ( keys %{$ref_transcript_len} ) {
            if ( defined $ref_transcript_uniqRead->{$tran} ) {
                $ref_transcript_rpkm->{$tran} = ( ( $ref_transcript_uniqRead->{$tran} ) * 1000000000 ) / ( $ref_transcript_len->{$tran} * $totalRead );
            }
        }
    }

    1;
}

sub output_rpkm  {
    my ( $outputFile, $cutoff, $ref_transcript_len, $ref_transcript_uniqRead, $ref_transcript_multiRead, $ref_transcript_rpkm ) = @_;
    print STDERR "Output transcript abundance in RPKM to file $outputFile...\n\t", `date`;
    
    open ( OUT, ">>$outputFile" ) or die "Error opening $outputFile for writing transcript RPKM.\n";
    print OUT "#transcript\tlength\tunique_read_num\tmultiple_read_num\tRPKM\n";
    foreach my $tran ( sort { $ref_transcript_rpkm->{$b} <=> $ref_transcript_rpkm->{$a} } ( keys %{$ref_transcript_len} ) ) {
        last if ( $ref_transcript_rpkm->{$tran} < $cutoff );
        print OUT $tran, "\t", $ref_transcript_len->{$tran}, "\t", $ref_transcript_uniqRead->{$tran}, "\t", $ref_transcript_multiRead->{$tran}, "\t", $ref_transcript_rpkm->{$tran}, "\n";
    }
    close OUT;

    1;
}

