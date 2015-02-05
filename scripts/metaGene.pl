#! /usr/bin/perl
#
use strict;
use warnings;

use lib "/home/qczhang/lib/perllib";
use Data::Dumper;
use Getopt::Std;
use Locale::Bio::IO qw( &readGTF &readFasta );

use vars qw ($opt_h $opt_V $opt_D $opt_i $opt_o $opt_g $opt_f $opt_s $opt_t $opt_c $opt_C $opt_S $opt_p );
&getopts('hVDi:o:g:f:s:t:c:C:S:p:');

my $usage = <<_EOH_;
## --------------------------------------
meta gene analysis

Command:
$0 -i transcript_file -o metagene_file 

# what it is:
 -i     profile for every transcript
 -o     profile of meta gene

# more options:
 -g     gtf file 
        default: /home/qczhang/database/ensembl/current/mouse/gtf/Mus_musculus.GRCm38.74.gtf 
 -f     gtf file format
        default: ensembl

 -s     sequence file
        e.g. /home/qczhang/database/ensembl/current/mouse/gtf/Mus_musculus.GRCm38.74.dna_sm.primary_assembly.GRCm38.74.gtf.collapsed.shortHeader.rRNA.nr.shortHeader.fa
 -t     translation efficiency file
        e.g. /home/qczhang/shape-seq/new/hiseq1/riboprofiling_results/ingolia2011Cell.txt

 -c     transcript column
        default: 1
 -C     starting profile column
        default: 2 (if a column has multiple fields, only the first will be used)
 -p     rpkm column

 -S     statistics file (not used now)

_EOH_
;

my %parameters = &init();

my $transcriptFile = $parameters{transcriptFile};
my $transcriptCol = $parameters{transcriptCol};
my $profileCol = $parameters{profileCol};
my $rpkmCol = defined ( $parameters{rpkmCol} ) ? $parameters{rpkmCol} : 0;
my $metageneFile = $parameters{metageneFile};

my $gtfFile = $parameters{gtfFile};
my $gtfSource = $parameters{gtfSource};

my $ref_other_ucsc = undef; my $ref_trans_efficiency = undef;  my $ref_rb_mRNA = undef;
my @avgNum_highRibo = undef; my @avgProfile_highRibo = undef;
my @avgNum_lowRibo = undef; my @avgProfile_lowRibo = undef;
if ( defined $parameters{translationFile} ) { 
    $ref_other_ucsc = getIDmap2UCSC (); 
    ( $ref_trans_efficiency, $ref_rb_mRNA ) = readRiboEfficiency ( $parameters{translationFile} ); 
    @avgNum_highRibo = newArray ( 300 ); @avgProfile_highRibo = newArray ( 300 );
    @avgNum_lowRibo = newArray ( 300 ); @avgProfile_lowRibo = newArray ( 300 );
}
my $ref_trans_sequence = undef; my $ref_trans_len = undef; my $ref_trans_alias = undef;
my @avgNum_ac = undef; my @avgProfile_ac = undef;
if ( defined $parameters{sequenceFile} ) { 
    ( $ref_trans_sequence, $ref_trans_len, $ref_trans_alias ) = readFasta ( $parameters{sequenceFile} ); 
    @avgNum_ac = newArray ( 300 ); @avgProfile_ac = newArray ( 300 );
}

my $ref_annotation = readGTF ( $gtfFile, attribute => "transcript_id", source => $gtfSource, verbose => 1, skip => "Blat" ); smallfix ( $ref_annotation );
my @avgNum = newArray ( 300 ); my @avgProfile = newArray ( 300 );
my @avgNum_long = newArray ( 300 ); my @avgProfile_long = newArray ( 300 );
my @avgNum_short = newArray ( 300 ); my @avgProfile_short = newArray ( 300 );

my %uniqBioType = ();
my $lineCount = 0;
open ( RPKM, $transcriptFile );
while ( my $line = <RPKM> ) {
    next if ( $line =~ /^#/ );
    next if ( $line =~ /RTstop/ );
    $line =~ s/^\s+//; chomp $line;
    $lineCount++; #last if ( $lineCount > 10 );

    my @data = split ( /\s/, $line );
    my $transID = $data[$transcriptCol];
    my $rpkm = 0;  $rpkm = $data[$rpkmCol] if ( $rpkmCol );
    my @profile = @data[$profileCol..$#data];
    if ( defined $parameters{sequenceFile} )  {
        my $len = $data[1];
        if ( $len != $ref_trans_len->{$transID} ) {
            print STDERR "ERROR! length og $transID not consistent in sequence ($ref_trans_len->{$transID}) and RT stop ($len) files!\n";
            next;
        }
    }
    my $biotypeString = getBioType ( $transID, $ref_annotation, \%uniqBioType );
    print $transID, "\t", $biotypeString, "\n" if ( $opt_V );

    if ( $biotypeString eq "protein coding" ) {
        if ( ( defined $ref_annotation->{$transID}{start_codon} ) and ( defined $ref_annotation->{$transID}{stop_codon} ) ) {
            &printAndStat ( $parameters{statsFile}, $transID, $rpkm, \@profile, $ref_annotation, $ref_trans_sequence, \@avgNum, \@avgProfile, \@avgNum_ac, \@avgProfile_ac, \@avgNum_highRibo, \@avgProfile_highRibo, \@avgNum_lowRibo, \@avgProfile_lowRibo );
        }
        else {
            print STDERR "Warning! protein coding transcript $transID has no annotation!\n" if ( $opt_V );
        }
    }
}
close RPKM;

&printAverage ( $metageneFile, "Average", \@avgNum, \@avgProfile );
&printAverage ( $metageneFile, "Average_long", \@avgNum_long, \@avgProfile_long );
&printAverage ( $metageneFile, "Average_short", \@avgNum_short, \@avgProfile_short );
if ( defined $parameters{sequenceFile} ) { 
    &printAverage ( $metageneFile, "Average_AC", \@avgNum_ac, \@avgProfile_ac );
}
if ( defined $parameters{translationFile} ) { 
    &printAverage ( $metageneFile, "Average_highRibo", \@avgNum_highRibo, \@avgProfile_highRibo );
    &printAverage ( $metageneFile, "Average_lowRibo", \@avgNum_lowRibo, \@avgProfile_lowRibo );
}


## ------------------------------------
sub init
{
    my %parameters = ();
    die $usage if ( $opt_h || ( not $opt_i ) || ( not $opt_o ) );

    $opt_V = 0 if ( not defined $opt_V );
    $opt_D = 0 if ( not defined $opt_D );
    my $pwd = `pwd`;  chomp $pwd;

    $parameters{transcriptFile} = $opt_i;
    $parameters{metageneFile} = $opt_o;
    print STDERR `/bin/rm $parameters{metageneFile}` if ( -e $parameters{metageneFile} );
    if ( defined $opt_c ) { $parameters{transcriptCol} = $opt_c - 1; }
    else { $parameters{transcriptCol} = 0; }
    if ( defined $opt_C ) { $parameters{profileCol} = $opt_C - 1; }
    else { $parameters{profileCol} = 1; }
    if ( defined $opt_p ) { $parameters{rpkmCol} = $opt_p - 1; }

    if ( defined $opt_g ) { $parameters{gtfFile} = $opt_g; }
    else { $parameters{gtfFile} = "/home/qczhang/database/ensembl/current/mouse/gtf/Mus_musculus.GRCm38.74.gtf"; }
    if ( defined $opt_f ) { $parameters{gtfSource} = $opt_f; }
    else { $parameters{gtfSource} = "ensembl"; }

    if ( defined $opt_s ) { $parameters{sequenceFile} = $opt_s; }
    if ( defined $opt_S ) { $parameters{statsFile} = $opt_S; print STDERR `/bin/rm $parameters{statsFile}` if ( -e $parameters{statsFile} );  }
    if ( defined $opt_t ) { $parameters{translationFile} = $opt_t; }

    print Dumper \%parameters if ( $opt_D );
    return ( %parameters );
}


sub newArray 
{
    my $len = shift;
    my %parameters = @_;

    my @array = ();
    for ( my $idx = 0; $idx < 300; $idx++ ) {
        if ( defined $parameters{element} ) {
            $array[$idx] = $parameters{element};
        }
        else {
            $array[$idx] = 0;
        }
    }

    if ( ( defined $parameters{returnHash} ) and ( $parameters{returnHash} ) ) {
        return \@array;
    }
    else {
        return @array;
    }
}

sub printAndStat 
{
    my $statsFile = shift;
    my $transID = shift; my $rpkm = shift;  my $ref_profile = shift;
    my $ref_annotation = shift;  my $ref_trans_sequence = shift;
    my $ref_avgNum = shift;  my $ref_avgProfile = shift;
    my $ref_avgNum_ac = shift;  my $ref_avgProfile_ac = shift;
    my $ref_avgNum_highRibo = shift;  my $ref_avgProfile_highRibo = shift;
    my $ref_avgNum_lowRibo = shift;  my $ref_avgProfile_lowRibo = shift;

    my $exonLen = &getExonLen ( $ref_annotation->{$transID} );
    my $fivePrimeLen = &get5primeLen ( $ref_annotation->{$transID} );
    my $threePrimeLen = &get3primeLen ( $ref_annotation->{$transID} );
    my $CDSlen = $exonLen - $fivePrimeLen - $threePrimeLen;
    my $TSSpos = $fivePrimeLen;  my $TESpos = $fivePrimeLen + $CDSlen;
    ## check strand ## note that position are 0-indexed
    if ( defined $statsFile ) {
        open ( OUT, ">>$statsFile" );
        print OUT $transID, "\t", $rpkm, "\t", $exonLen, "\t", $fivePrimeLen, "\t", $CDSlen, "\t", $threePrimeLen;
    }

    my $riboEffi = "NULL";  my $mRNA = "NULL";
    if ( defined $parameters{translationFile} ) { ( $riboEffi, $mRNA ) = &getRiboEfficiency ( $transID, $ref_other_ucsc, $ref_trans_efficiency, $ref_rb_mRNA ); }
    if ( defined $statsFile ) { print OUT "\t", $riboEffi, "\t", $mRNA; }

    my @bases = undef;
    if ( defined $parameters{sequenceFile} ) { @bases = split ( //, uc($ref_trans_sequence->{$transID}) ); }

    if ( $fivePrimeLen < 50 ) {
        print STDERR "$transID 5 prime length less than 50!\n" if ( $opt_V );
        if ( defined $statsFile ) { for ( my $posIdx = 0; $posIdx < ( 50 - $fivePrimeLen ); $posIdx++ ) { print OUT "\tNULL"; } }
        for ( my $posIdx = 0; $posIdx < $TSSpos; $posIdx++ ) {
            if ( defined $statsFile ) { print OUT "\t", $ref_profile->[$posIdx]; }
            if ( $ref_profile->[$posIdx] ne "NULL" ) {
                $avgNum[50-$fivePrimeLen+$posIdx]++;
                $avgProfile[50-$fivePrimeLen+$posIdx] += $ref_profile->[$posIdx];

                if ( ( $riboEffi ne "NULL" ) and ( $CDSlen < 1500 ) and ( $CDSlen > 800 ) ) {
                    if ( $riboEffi < 0.4 ) {
                        $avgNum_lowRibo[50-$fivePrimeLen+$posIdx]++;
                        $avgProfile_lowRibo[50-$fivePrimeLen+$posIdx] += $ref_profile->[$posIdx];
                    }
                    elsif ( $riboEffi > 1 ) {
                        $avgNum_highRibo[50-$fivePrimeLen+$posIdx]++;
                        $avgProfile_highRibo[50-$fivePrimeLen+$posIdx] += $ref_profile->[$posIdx];
                    }
                }

                if ( $CDSlen < 400 ) {
                    $avgNum_short[50-$fivePrimeLen+$posIdx]++;
                    $avgProfile_short[50-$fivePrimeLen+$posIdx] += $ref_profile->[$posIdx];
                }
                elsif ( $CDSlen > 2000 ) {
                    $avgNum_long[50-$fivePrimeLen+$posIdx]++;
                    $avgProfile_long[50-$fivePrimeLen+$posIdx] += $ref_profile->[$posIdx];
                }

                if ( defined $parameters{sequenceFile} ) {
                    if ( ( $bases[$posIdx] eq "A" ) or ( $bases[$posIdx] eq "C" ) ) {
                        $avgNum_ac[50-$fivePrimeLen+$posIdx]++;
                        $avgProfile_ac[50-$fivePrimeLen+$posIdx] += $ref_profile->[$posIdx];
                    }
                }
            }
        }
    }
    else {
        for ( my $posIdx = ( $TSSpos - 50 ); $posIdx < $TSSpos; $posIdx++ ) {
            if ( defined $statsFile ) { print OUT "\t", $ref_profile->[$posIdx]; }
            if ( $ref_profile->[$posIdx] ne "NULL" ) {
                $avgNum[50-$fivePrimeLen+$posIdx]++;
                $avgProfile[50-$fivePrimeLen+$posIdx] += $ref_profile->[$posIdx];

                if ( ( $riboEffi ne "NULL" ) and ( $CDSlen < 1500 ) and ( $CDSlen > 800 ) ) {
                    if ( $riboEffi < 0.4 ) {
                        $avgNum_lowRibo[50-$fivePrimeLen+$posIdx]++;
                        $avgProfile_lowRibo[50-$fivePrimeLen+$posIdx] += $ref_profile->[$posIdx];
                    }
                    elsif ( $riboEffi > 1 ) {
                        $avgNum_highRibo[50-$fivePrimeLen+$posIdx]++;
                        $avgProfile_highRibo[50-$fivePrimeLen+$posIdx] += $ref_profile->[$posIdx];
                    }
                }

                if ( $CDSlen < 400 ) {
                    $avgNum_short[50-$fivePrimeLen+$posIdx]++;
                    $avgProfile_short[50-$fivePrimeLen+$posIdx] += $ref_profile->[$posIdx];
                }
                elsif ( $CDSlen > 2000 ) {
                    $avgNum_long[50-$fivePrimeLen+$posIdx]++;
                    $avgProfile_long[50-$fivePrimeLen+$posIdx] += $ref_profile->[$posIdx];
                }

                if ( defined $parameters{sequenceFile} ) {
                    if ( ( $bases[$posIdx] eq "A" ) or ( $bases[$posIdx] eq "C" ) ) {
                        $avgNum_ac[50-$fivePrimeLen+$posIdx]++;
                        $avgProfile_ac[50-$fivePrimeLen+$posIdx] += $ref_profile->[$posIdx];
                    }
                }
            }
        }
    }

    if ( $CDSlen < 100 ) {
        print STDERR "$transID CDS length less than 50!\n" if ( $opt_V );
        for ( my $posIdx = $TSSpos; $posIdx < ( $TSSpos + $CDSlen ); $posIdx++ ) {
            if ( defined $statsFile ) { print OUT "\t", $ref_profile->[$posIdx]; }
            if ( $ref_profile->[$posIdx] ne "NULL" ) {
                $avgNum[50-$fivePrimeLen+$posIdx]++;
                $avgProfile[50-$fivePrimeLen+$posIdx] += $ref_profile->[$posIdx];

                if ( ( $riboEffi ne "NULL" ) and ( $CDSlen < 1500 ) and ( $CDSlen > 800 ) ) {
                    if ( $riboEffi < 0.4 ) {
                        $avgNum_lowRibo[50-$fivePrimeLen+$posIdx]++;
                        $avgProfile_lowRibo[50-$fivePrimeLen+$posIdx] += $ref_profile->[$posIdx];
                    }
                    elsif ( $riboEffi > 1 ) {
                        $avgNum_highRibo[50-$fivePrimeLen+$posIdx]++;
                        $avgProfile_highRibo[50-$fivePrimeLen+$posIdx] += $ref_profile->[$posIdx];
                    }
                }

                if ( $CDSlen < 400 ) {
                    $avgNum_short[50-$fivePrimeLen+$posIdx]++;
                    $avgProfile_short[50-$fivePrimeLen+$posIdx] += $ref_profile->[$posIdx];
                }
                elsif ( $CDSlen > 2000 ) {
                    $avgNum_long[50-$fivePrimeLen+$posIdx]++;
                    $avgProfile_long[50-$fivePrimeLen+$posIdx] += $ref_profile->[$posIdx];
                }

                if ( defined $parameters{sequenceFile} ) {
                    if ( ( $bases[$posIdx] eq "A" ) or ( $bases[$posIdx] eq "C" ) ) {
                        $avgNum_ac[50-$fivePrimeLen+$posIdx]++;
                        $avgProfile_ac[50-$fivePrimeLen+$posIdx] += $ref_profile->[$posIdx];
                    }
                }
            }
        }
        if ( defined $statsFile ) {
            for ( my $posIdx = 0; $posIdx < ( 100 - $CDSlen ); $posIdx++ ) { print OUT "\tNULL"; }
            for ( my $posIdx = 0; $posIdx < ( 100 - $CDSlen ); $posIdx++ ) { print OUT "\tNULL"; }
        }
        for ( my $posIdx = ( $TESpos - 100 ); $posIdx < $TESpos; $posIdx++ ) {
            if ( defined $statsFile ) { print OUT "\t", $ref_profile->[$posIdx]; }
            if ( $ref_profile->[$posIdx] ne "NULL" ) {
                $avgNum[250-$fivePrimeLen-$CDSlen+$posIdx]++;
                $avgProfile[250-$fivePrimeLen-$CDSlen+$posIdx] += $ref_profile->[$posIdx];

                if ( ( $riboEffi ne "NULL" ) and ( $CDSlen < 1500 ) and ( $CDSlen > 800 ) ) {
                    if ( $riboEffi < 0.4 ) {
                        $avgNum_lowRibo[250-$fivePrimeLen-$CDSlen+$posIdx]++;
                        $avgProfile_lowRibo[250-$fivePrimeLen-$CDSlen+$posIdx] += $ref_profile->[$posIdx];
                    }
                    elsif ( $riboEffi > 1 ) {
                        $avgNum_highRibo[250-$fivePrimeLen-$CDSlen+$posIdx]++;
                        $avgProfile_highRibo[250-$fivePrimeLen-$CDSlen+$posIdx] += $ref_profile->[$posIdx];
                    }
                }

                if ( $CDSlen < 400 ) {
                    $avgNum_short[250-$fivePrimeLen-$CDSlen+$posIdx]++;
                    $avgProfile_short[250-$fivePrimeLen-$CDSlen+$posIdx] += $ref_profile->[$posIdx];
                }
                elsif ( $CDSlen > 2000 ) {
                    $avgNum_long[250-$fivePrimeLen-$CDSlen+$posIdx]++;
                    $avgProfile_long[250-$fivePrimeLen-$CDSlen+$posIdx] += $ref_profile->[$posIdx];
                }

                if ( defined $parameters{sequenceFile} ) {
                    if ( ( $bases[$posIdx] eq "A" ) or ( $bases[$posIdx] eq "C" ) ) {
                        $avgNum_ac[250-$fivePrimeLen-$CDSlen+$posIdx]++;
                        $avgProfile_ac[250-$fivePrimeLen-$CDSlen+$posIdx] += $ref_profile->[$posIdx];
                    }
                }
            }
        }
    }
    else {
        for ( my $posIdx = $TSSpos; $posIdx < ( $TSSpos + 100 ); $posIdx++ ) {
            if ( defined $statsFile ) { print OUT "\t", $ref_profile->[$posIdx]; }
            if ( $ref_profile->[$posIdx] ne "NULL" ) {
                $avgNum[50-$fivePrimeLen+$posIdx]++;
                $avgProfile[50-$fivePrimeLen+$posIdx] += $ref_profile->[$posIdx];

                if ( ( $riboEffi ne "NULL" ) and ( $CDSlen < 1500 ) and ( $CDSlen > 800 ) ) {
                    if ( $riboEffi < 0.4 ) {
                        $avgNum_lowRibo[50-$fivePrimeLen+$posIdx]++;
                        $avgProfile_lowRibo[50-$fivePrimeLen+$posIdx] += $ref_profile->[$posIdx];
                    }
                    elsif ( $riboEffi > 1 ) {
                        $avgNum_highRibo[50-$fivePrimeLen+$posIdx]++;
                        $avgProfile_highRibo[50-$fivePrimeLen+$posIdx] += $ref_profile->[$posIdx];
                    }
                }

                if ( $CDSlen < 400 ) {
                    $avgNum_short[50-$fivePrimeLen+$posIdx]++;
                    $avgProfile_short[50-$fivePrimeLen+$posIdx] += $ref_profile->[$posIdx];
                }
                elsif ( $CDSlen > 2000 ) {
                    $avgNum_long[50-$fivePrimeLen+$posIdx]++;
                    $avgProfile_long[50-$fivePrimeLen+$posIdx] += $ref_profile->[$posIdx];
                }

                if ( defined $parameters{sequenceFile} ) {
                    if ( ( $bases[$posIdx] eq "A" ) or ( $bases[$posIdx] eq "C" ) ) {
                        $avgNum_ac[50-$fivePrimeLen+$posIdx]++;
                        $avgProfile_ac[50-$fivePrimeLen+$posIdx] += $ref_profile->[$posIdx];
                    }
                }
            }
        }
        for ( my $posIdx = ( $TESpos - 100 ); $posIdx < $TESpos; $posIdx++ ) {
            if ( defined $statsFile ) { print OUT "\t", $ref_profile->[$posIdx]; }
            if ( $ref_profile->[$posIdx] ne "NULL" ) {
                $avgNum[250-$fivePrimeLen-$CDSlen+$posIdx]++;
                $avgProfile[250-$fivePrimeLen-$CDSlen+$posIdx] += $ref_profile->[$posIdx];

                if ( ( $riboEffi ne "NULL" ) and ( $CDSlen < 1500 ) and ( $CDSlen > 800 ) ) {
                    if ( $riboEffi < 0.4 ) {
                        $avgNum_lowRibo[250-$fivePrimeLen-$CDSlen+$posIdx]++;
                        $avgProfile_lowRibo[250-$fivePrimeLen-$CDSlen+$posIdx] += $ref_profile->[$posIdx];
                    }
                    elsif ( $riboEffi > 1 ) {
                        $avgNum_highRibo[250-$fivePrimeLen-$CDSlen+$posIdx]++;
                        $avgProfile_highRibo[250-$fivePrimeLen-$CDSlen+$posIdx] += $ref_profile->[$posIdx];
                    }
                }

                if ( $CDSlen < 400 ) {
                    $avgNum_short[250-$fivePrimeLen-$CDSlen+$posIdx]++;
                    $avgProfile_short[250-$fivePrimeLen-$CDSlen+$posIdx] += $ref_profile->[$posIdx];
                }
                elsif ( $CDSlen > 2000 ) {
                    $avgNum_long[250-$fivePrimeLen-$CDSlen+$posIdx]++;
                    $avgProfile_long[250-$fivePrimeLen-$CDSlen+$posIdx] += $ref_profile->[$posIdx];
                }

                if ( defined $parameters{sequenceFile} ) {
                    if ( ( $bases[$posIdx] eq "A" ) or ( $bases[$posIdx] eq "C" ) ) {
                        $avgNum_ac[250-$fivePrimeLen-$CDSlen+$posIdx]++;
                        $avgProfile_ac[250-$fivePrimeLen-$CDSlen+$posIdx] += $ref_profile->[$posIdx];
                    }
                }
            }
        }
    }

    if ( $threePrimeLen < 50 ) {
        print STDERR "$transID 3 prime length less than 50!\n" if ( $opt_V );
        for ( my $posIdx = $TESpos; $posIdx < ( $TESpos + $threePrimeLen ); $posIdx++ ) {
            if ( defined $statsFile ) { print OUT "\t", $ref_profile->[$posIdx]; }
            if ( $ref_profile->[$posIdx] ne "NULL" ) {
                $avgNum[250-$fivePrimeLen-$CDSlen+$posIdx]++;
                $avgProfile[250-$fivePrimeLen-$CDSlen+$posIdx] += $ref_profile->[$posIdx];

                if ( ( $riboEffi ne "NULL" ) and ( $CDSlen < 1500 ) and ( $CDSlen > 800 ) ) {
                    if ( $riboEffi < 0.4 ) {
                        $avgNum_lowRibo[250-$fivePrimeLen-$CDSlen+$posIdx]++;
                        $avgProfile_lowRibo[250-$fivePrimeLen-$CDSlen+$posIdx] += $ref_profile->[$posIdx];
                    }
                    elsif ( $riboEffi > 1 ) {
                        $avgNum_highRibo[250-$fivePrimeLen-$CDSlen+$posIdx]++;
                        $avgProfile_highRibo[250-$fivePrimeLen-$CDSlen+$posIdx] += $ref_profile->[$posIdx];
                    }
                }

                if ( $CDSlen < 400 ) {
                    $avgNum_short[250-$fivePrimeLen-$CDSlen+$posIdx]++;
                    $avgProfile_short[250-$fivePrimeLen-$CDSlen+$posIdx] += $ref_profile->[$posIdx];
                }
                elsif ( $CDSlen > 2000 ) {
                    $avgNum_long[250-$fivePrimeLen-$CDSlen+$posIdx]++;
                    $avgProfile_long[250-$fivePrimeLen-$CDSlen+$posIdx] += $ref_profile->[$posIdx];
                }

                if ( defined $parameters{sequenceFile} ) {
                    if ( ( $bases[$posIdx] eq "A" ) or ( $bases[$posIdx] eq "C" ) ) {
                        $avgNum_ac[250-$fivePrimeLen-$CDSlen+$posIdx]++;
                        $avgProfile_ac[250-$fivePrimeLen-$CDSlen+$posIdx] += $ref_profile->[$posIdx];
                    }
                }
            }
        }
        if ( defined $statsFile ) {
            for ( my $posIdx = 0; $posIdx < ( 50 - $threePrimeLen ); $posIdx++ ) { print OUT "\tNULL"; }
        }
    }
    else {
        for ( my $posIdx = $TESpos; $posIdx < ( $TESpos + 50 ); $posIdx++ ) {
            if ( defined $statsFile ) { print OUT "\t", $ref_profile->[$posIdx]; }
            if ( $ref_profile->[$posIdx] ne "NULL" ) {
                $avgNum[250-$fivePrimeLen-$CDSlen+$posIdx]++;
                $avgProfile[250-$fivePrimeLen-$CDSlen+$posIdx] += $ref_profile->[$posIdx];

                if ( ( $riboEffi ne "NULL" ) and ( $CDSlen < 1500 ) and ( $CDSlen > 800 ) ) {
                    if ( $riboEffi < 0.4 ) {
                        $avgNum_lowRibo[250-$fivePrimeLen-$CDSlen+$posIdx]++;
                        $avgProfile_lowRibo[250-$fivePrimeLen-$CDSlen+$posIdx] += $ref_profile->[$posIdx];
                    }
                    elsif ( $riboEffi > 1 ) {
                        $avgNum_highRibo[250-$fivePrimeLen-$CDSlen+$posIdx]++;
                        $avgProfile_highRibo[250-$fivePrimeLen-$CDSlen+$posIdx] += $ref_profile->[$posIdx];
                    }
                }

                if ( $CDSlen < 400 ) {
                    $avgNum_short[250-$fivePrimeLen-$CDSlen+$posIdx]++;
                    $avgProfile_short[250-$fivePrimeLen-$CDSlen+$posIdx] += $ref_profile->[$posIdx];
                }
                elsif ( $CDSlen > 2000 ) {
                    $avgNum_long[250-$fivePrimeLen-$CDSlen+$posIdx]++;
                    $avgProfile_long[250-$fivePrimeLen-$CDSlen+$posIdx] += $ref_profile->[$posIdx];
                }

                if ( defined $parameters{sequenceFile} ) {
                    if ( ( $bases[$posIdx] eq "A" ) or ( $bases[$posIdx] eq "C" ) ) {
                        $avgNum_ac[250-$fivePrimeLen-$CDSlen+$posIdx]++;
                        $avgProfile_ac[250-$fivePrimeLen-$CDSlen+$posIdx] += $ref_profile->[$posIdx];
                    }
                }
            }
        }
    }

    if ( defined $statsFile ) { print OUT "\n"; close OUT;  }
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

sub getRiboEfficiency 
{
    my $transID = shift;
    my $ref_other_ucsc = shift;
    my $ref_trans_efficiency = shift;
    my $ref_trans_mRNA = shift;

    my $riboEffiCount = 0; my $riboEffi = 0; my $mRNA = 0;
    my ( $transIDsim, @transIDappendix ) = split ( /\./, $transID );
    if ( defined $ref_other_ucsc->{$transIDsim} ) {
        my @ucscIDs = split ( /\|/, $ref_other_ucsc->{$transIDsim} );
        foreach my $uid ( @ucscIDs ) {
            if ( defined $ref_trans_efficiency->{$uid} ) {
                $riboEffiCount++;
                $riboEffi += $ref_trans_efficiency->{$uid};
                $mRNA += $ref_trans_mRNA->{$uid};
            }
        }
    }
    if ( $riboEffiCount ) {
        $riboEffi /= $riboEffiCount;
        $mRNA /= $riboEffiCount;
    }
    else {
        $riboEffi = "NULL";
        $mRNA = "NULL";
    }

    return ( $riboEffi, $mRNA );
}

sub readRiboEfficiency 
{
    my $riboEfficiencyFile = shift;

    my %trans_riboEffc = ();
    my %trans_mRNA = ();
    open ( RIBO, $riboEfficiencyFile ) or die "Error openning ribosomal profiling file $riboEfficiencyFile!\n";
    while ( my $line = <RIBO> ) {
        next if ( $line =~ /^#/ );

        chomp $line;
        my ( $ucscID, $riboDensity, $mRNADensity, $TEratio, $logTEratio, $gene, $desc ) = split ( /\t/, $line );
        $trans_riboEffc{$ucscID} = $logTEratio;
        $trans_mRNA{$ucscID} = $mRNADensity;
    }
    close RIBO;

    return ( \%trans_riboEffc, \%trans_mRNA );
}

sub getIDmap2UCSC 
{
    my %other_ucsc = ();

    my $ensembl2ucscFile = "/home/qczhang/projects2/databases/ucsc/mm9/ucsc_ensembl.map";
    open ( EU, $ensembl2ucscFile );
    while ( my $line = <EU> ) {
        next if ( $line =~ /^#/ );

        chomp $line;
        my ( $ucscID, $ensemblID ) = split ( /\t/, $line );
        if ( not defined $other_ucsc{$ensemblID} ) {
            $other_ucsc{$ensemblID} = $ucscID;
        }
        else {
            $other_ucsc{$ensemblID} .= "|" . $ucscID;
        }
    }
    close EU;

    my $refseq2ucscFile = "/home/qczhang/projects2/databases/ucsc/mm9/ucsc_refseq.map";
    open ( RU, $refseq2ucscFile );
    while ( my $line = <RU> ) {
        next if ( $line =~ /^#/ );

        chomp $line;
        my ( $ucscID, $refSeqID ) = split ( /\t/, $line );
        if ( not defined $other_ucsc{$refSeqID} ) {
            $other_ucsc{$refSeqID} = $ucscID;
        }
        else {
            $other_ucsc{$refSeqID} .= "|" . $ucscID;
        }
    }
    close RU;

    my $vega2ucscFile = "/home/qczhang/database/mgi/gtf/vega_others.map";
    open ( VU, $vega2ucscFile );
    while ( my $line = <VU> ) {
        next if ( $line =~ /^#/ );

        chomp $line;
        my ( $vegaID, $otherIDs ) = split ( /\t/, $line );
        my @otherIDList = split ( /\|/, $otherIDs );
        foreach my $otherID ( @otherIDList ) {
            my ( $otherIDsim, @otherIDappendix ) = split ( /\./, $otherID );
            if ( defined $other_ucsc{$otherIDsim} ) {
                if ( not defined $other_ucsc{$vegaID} ) {
                    $other_ucsc{$vegaID} = $other_ucsc{$otherIDsim};
                }
                else {
                    $other_ucsc{$vegaID} .= "|" . $other_ucsc{$otherIDsim};
                }
            }
        }
    }
    close VU;

    my $ncbiGene2ucscFile = "/home/qczhang/database/mgi/gtf/ncbiGene_others.map";
    open ( NU, $ncbiGene2ucscFile );
    while ( my $line = <NU> ) {
        next if ( $line =~ /^#/ );

        chomp $line;
        my ( $ncbiGeneID, $otherIDs ) = split ( /\t/, $line );
        my @otherIDList = split ( /\|/, $otherIDs );
        foreach my $otherID ( @otherIDList ) {
            my ( $otherIDsim, @otherIDappendix ) = split ( /\./, $otherID );
            if ( not $otherIDsim ) {
                print STDERR $line, "\n" if ( $opt_V );
            }
            if ( defined $other_ucsc{$otherIDsim} ) {
                if ( not defined $other_ucsc{$ncbiGeneID} ) {
                    $other_ucsc{$ncbiGeneID} = $other_ucsc{$otherIDsim};
                }
                else {
                    $other_ucsc{$ncbiGeneID} .= "|" . $other_ucsc{$otherIDsim};
                }
            }
        }
    }
    close NU;

    return \%other_ucsc;
}

sub getBioType 
{
    my $transID = shift;
    my $ref_annotation = shift;
    my $ref_uniqBioType = shift;

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

sub printAverage 
{
    my $outFile = shift;
    my $tag = shift;
    my $ref_numArray = shift;
    my $ref_avgArray = shift;

    open ( OUT, ">>$outFile" ) or die ( "Could not open file $outFile for reading!\n" );
    print OUT "$tag\tnum";
    for ( my $idx = 0; $idx < scalar ( @{$ref_numArray} ); $idx++ ) {
        print OUT "\t", $ref_numArray->[$idx];
    }
    print OUT "\n";

    print OUT "$tag\tenrich";
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

sub smallfix
{
    my $ref_annotation = shift;
    $ref_annotation->{MM5_8SrRNA}{biotype} = [ "rRNA gene" ]; $ref_annotation->{MM18SrRNA}{biotype} = [ "rRNA gene" ]; $ref_annotation->{MM28SrRNA}{biotype} = [ "rRNA gene" ];
    $ref_annotation->{"lcl|gene=robert:MM18SrRNA"}{biotype} = [ "rRNA gene" ]; $ref_annotation->{"lcl|gene=robert:MM28SrRNA"}{biotype} = [ "rRNA gene" ];
}
