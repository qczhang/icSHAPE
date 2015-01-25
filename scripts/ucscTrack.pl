#! /usr/bin/perl
#
use strict;
use warnings;

my $inFile = shift; ## must be sorted w.r.t position
my $uniqFile = shift;

my $oldChr = "",  my $oldStart = 0; my $oldScore = 0;  my $oldTrans = "";  my $oldNT = "";  my $posCount = 0;
open ( IN, $inFile );
open ( OUT, ">$uniqFile" );
my $lineCount = 0;
while ( my $line = <IN> ) {
    next if ( $line =~ /^#/ );
    chomp $line;

    $lineCount++;
    print STDERR "line: $lineCount\n" if ( $lineCount % 1000000 == 0 );
    my ( $chr, $start, $end, $score, $trans, $nt ) = split ( /\t/, $line );
    if ( ( $start != $oldStart ) or ( $chr ne $oldChr ) )  {
        $oldScore /= $posCount if ( $posCount );
        print OUT $oldChr, "\t", $oldStart, "\t", $oldStart+1, "\t", $oldScore, "\t", $oldTrans, "\t", $oldNT, "\n" if ( $oldStart );

        $oldChr = $chr;
        $oldStart = $start;
        $oldTrans = $trans;
        $oldNT = $nt;

        $oldScore = $score;
        if ( $score ne "NULL" ) { $posCount = 1; }
        else { $posCount = 0; }
    }
    else {
        if ( $score ne "NULL" ) {
            if ( $oldScore ne "NULL"  ) { $oldScore += $score; }
            else { $oldScore = $score; }
            $posCount++;
        }
        $oldTrans .= "," . $trans;
        $oldNT .= "," . $nt;
    }
}
print OUT $oldChr, "\t", $oldStart, "\t", $oldStart+1, "\t", $oldScore, "\t", $oldTrans, "\t", $oldNT, "\n";
close IN;
close OUT;
