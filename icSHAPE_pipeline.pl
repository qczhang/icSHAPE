#! /usr/bin/perl
# wrapper of icSHAPE pipeline
# copy right qiangfeng.zhang@gmail.com
# history: 0.01 
#   date: 01/06/2015

use strict;
use warnings;
use File::Basename;
use Getopt::Std;

my $_debug = 0;
my %config = ();

use vars qw ($opt_h $opt_V $opt_D $opt_i $opt_t $opt_o $opt_c );
&getopts('hVDi:t:o:c:');

my $usage = <<_EOH_;
## --------------------------------------

Command:
$0 -i input_file -t target_file 

# what it is:
 -i     input sequencing reads used as background
 -t     target sequencing reads of modified samples

# more options:
 -o     output directory
 -c     configuration file

_EOH_
;

&main ();

sub main
{
    &init ();
    if ( $_debug ) { foreach my $key ( keys %config ) { print $key, "\t", $config{$key}, "\n"; }  }

    my $outDir = $config{outDir};
    my @inputSeqFiles = split ( /:/, $config{input} ); my @targetSeqFiles = split ( /:/, $config{target} );
    my @inputSignalFile = (); my @targetSignalFile = ();
    foreach my $seqFile ( @inputSeqFiles, @targetSeqFiles ) {
        my ($fileName, $fileDir, $fileSuffix) = fileparse ( $seqFile, qr/\.[^.]*/ ); 
        my $seqCollapsed = $outDir . "/" . $fileName . ".rmdup.fastq";
        my $seqDatFasta = $outDir . "/" . $fileName . ".fa";
        if ( not -e "$seqCollapsed.done" ) {
            print STDERR "$config{COLLAPSEBIN} -U $seqFile -o $seqCollapsed -f $seqDatFasta\n";
            print STDERR `$config{COLLAPSEBIN} -U $seqFile -o $seqCollapsed -f $seqDatFasta`;
            if ( not $? ) {  print STDERR `touch $seqCollapsed.done`;  }
        }
        if ( not -e "$seqCollapsed.done" ) { die "Abort! icSHAPE pipeline die of unsuccessful read collapse.\n"; }

        my $seqTrimmed = $outDir . "/" . $fileName . ".trimmed.fastq";
        if ( not -e "$seqTrimmed.done" ) {
            print STDERR "$config{TRIMMER} -U $seqCollapsed -o $seqTrimmed -l $config{LEADINGTRIM} -t $config{TAILINGTRIM} -c $config{FASTQCODING} -a $config{ADAPTER} -m $config{TRIMMINLEN}\n";
            print STDERR `$config{TRIMMER} -U $seqCollapsed -o $seqTrimmed -l $config{LEADINGTRIM} -t $config{TAILINGTRIM} -c $config{FASTQCODING} -a $config{ADAPTER} -m $config{TRIMMINLEN}`;
            if ( not $? ) {  print STDERR `touch $seqTrimmed.done`;  }
        }
        if ( not -e "$seqTrimmed.done" ) { die "Abort! icSHAPE pipeline die of unsuccessful trimming.\n"; }

        my $samFile = $outDir . "/" . $fileName . ".sam";
        if ( not -e "$samFile.done" ) {
            my $alignOptions = ""; 
            $alignOptions .= "--" . $config{FASTQCODING} if ( defined $config{FASTQCODING} );
            $alignOptions .= " " . $config{MAPPINGOPTIONS} if ( defined $config{MAPPINGOPTIONS} );
            print STDERR "$config{ALIGNER} -U $seqTrimmed -S $samFile -x $config{MAPPINGREF} $alignOptions\n";
            print STDERR `$config{ALIGNER} -U $seqTrimmed -S $samFile -x $config{MAPPINGREF} $alignOptions`;
            if ( not $? ) {  print STDERR `touch $samFile.done`;  }
        }
        if ( not -e "$samFile.done" ) { die "Abort! icSHAPE pipeline die of unsuccessful alignment.\n"; }

        my $rpkmFile = $outDir . "/" . $fileName . ".rpkm";
        if ( not -e "$rpkmFile.done" ) {
            print STDERR "$config{ESTIMATERPKM} -i $samFile -o $rpkmFile\n";
            print STDERR `$config{ESTIMATERPKM} -i $samFile -o $rpkmFile`;
            if ( not $? ) {  print STDERR `touch $rpkmFile.done`;  }
        }
        if ( not -e "$rpkmFile.done" ) { die "Abort! icSHAPE pipeline die of unsuccessful abundance estimation.\n"; }

        # my $seqSumFile = $outDir . "/" . $fileName . ".seq.info";
        # print STDERR "$config{seqSummary} -i $samFile -o $seqSumFile -r $rpkmFile\n";
        # print STDERR `$config{seqSummary} -i $samFile -o $seqSumFile -r $rpkmFile`;

        my $rtFile = $outDir . "/" . $fileName . ".rt";
        if ( not -e "$rtFile.done" ) {
            print STDERR "$config{CALCRT} -i $samFile -o $rtFile -r $rpkmFile -c $config{MINLOAD}\n";
            print STDERR `$config{CALCRT} -i $samFile -o $rtFile -r $rpkmFile -c $config{MINLOAD}`;
            if ( not $? ) {  print STDERR `touch $rtFile.done`;  }
        }
        if ( not -e "$rtFile.done" ) { die "Abort! icSHAPE pipeline die of unsuccessful RT stop calculation.\n"; }

        push ( @inputSignalFile, $rtFile ) if grep { $_ eq $seqFile } @inputSeqFiles;
        push ( @targetSignalFile, $rtFile ) if grep { $_ eq $seqFile } @targetSeqFiles;
    }

    my $combinedInputSignalFile = $outDir . "/background.rt";
    if ( not -e "$combinedInputSignalFile.done" ) {
        if ( scalar @inputSignalFile > 1 ) {
            my $allInputSignalFile = "";
            foreach my $signalFile ( @inputSignalFile ) { $allInputSignalFile .= $signalFile . ":"; }
            $allInputSignalFile =~ s/:$//;
            print STDERR "$config{COMBINEBIN} -i $allInputSignalFile -o $combinedInputSignalFile";
            print STDERR `$config{COMBINEBIN} -i $allInputSignalFile -o $combinedInputSignalFile`;
        }
        else { print STDERR `ln -s $inputSignalFile[0] $combinedInputSignalFile` }
        if ( not $? ) {  print STDERR `touch $combinedInputSignalFile.done`; }
    }
    if ( not -e "$combinedInputSignalFile.done" ) { die "Abort! icSHAPE pipeline die of unsuccessful replicate combining.\n"; }

    my $normalizedInputSignalFile = $outDir . "/background.normalized.rt";
    if ( not -e "$normalizedInputSignalFile.done" ) {
        print STDERR "$config{NORMALIZEBIN} -i $combinedInputSignalFile -o $normalizedInputSignalFile -m $config{METHOD} -d $config{HEADTOSKIP} -l $config{TAILTOSKIP}\n";
        print STDERR `$config{NORMALIZEBIN} -i $combinedInputSignalFile -o $normalizedInputSignalFile -m $config{METHOD} -d $config{HEADTOSKIP} -l $config{TAILTOSKIP}`;
        if ( not $? ) {  print STDERR `touch $normalizedInputSignalFile.done`; }
    }
    if ( not -e "$normalizedInputSignalFile.done" ) { die "Abort! icSHAPE pipeline die of unsuccessful normalization.\n"; }

    my $combinedTargetSignalFile = $outDir . "/target.rt";
    if ( not -e "$combinedTargetSignalFile.done" ) {
        if ( scalar @targetSignalFile > 1 ) {
            my $allTargetSignalFile = "";
            foreach my $signalFile ( @targetSignalFile ) { $allTargetSignalFile .= $signalFile . ":"; }
            $allTargetSignalFile =~ s/:$//;
            print STDERR "$config{COMBINEBIN} -i $allTargetSignalFile -o $combinedTargetSignalFile";
            print STDERR `$config{COMBINEBIN} -i $allTargetSignalFile -o $combinedTargetSignalFile`;
        }
        else { print STDERR `ln -s $targetSignalFile[0] $combinedTargetSignalFile` }
        if ( not $? ) {  print STDERR `touch $combinedTargetSignalFile.done`; }
    }
    if ( not -e "$combinedTargetSignalFile.done" ) { die "Abort! icSHAPE pipeline die of unsuccessful replicate combining.\n"; }

    my $normalizedTargetSignalFile = $outDir . "/target.normalized.rt";
    if ( not -e "$normalizedTargetSignalFile.done" ) {
        print STDERR "$config{NORMALIZEBIN} -i $combinedTargetSignalFile -o $normalizedTargetSignalFile -m $config{METHOD} -d $config{HEADTOSKIP} -l $config{TAILTOSKIP}\n";
        print STDERR `$config{NORMALIZEBIN} -i $combinedTargetSignalFile -o $normalizedTargetSignalFile -m $config{METHOD} -d $config{HEADTOSKIP} -l $config{TAILTOSKIP}`;
        if ( not $? ) {  print STDERR `touch $normalizedTargetSignalFile.done`; }
    }
    if ( not -e "$normalizedTargetSignalFile.done" ) { die "Abort! icSHAPE pipeline die of unsuccessful normalization.\n"; }

    ## calculate icSHAPE enrichment scores and filter for valid ones
    my $icShapeAllFile = $outDir . "/icshape.tmp.out";
    if ( not -e "$icShapeAllFile.done" ) {
        print STDERR "$config{CALCENRICHBIN} -f $normalizedTargetSignalFile -b $normalizedInputSignalFile -o $icShapeAllFile -w $config{WINSOR} -x $config{DIVFACTOR} -y $config{SUBFACTOR}\n";
        print STDERR `$config{CALCENRICHBIN} -f $normalizedTargetSignalFile -b $normalizedInputSignalFile -o $icShapeAllFile -w $config{WINSOR} -x $config{DIVFACTOR} -y $config{SUBFACTOR}`;
        if ( not $? ) {  print STDERR `touch $icShapeAllFile.done`; }
    }
    if ( not -e "$icShapeAllFile.done" ) { die "Abort! icSHAPE pipeline die of unsuccessful icSHAPE score calculation.\n"; }

    my $icShapeFile = $outDir . "/icshape.out";
    if ( not -e "$icShapeFile.done" ) {
        print STDERR "$config{FILTERENRICH} -i $icShapeAllFile -o $icShapeFile -t $config{INPUTCOVERAGE} -T $config{TARGETHIT} -s $config{HEADTOSKIP} -e $config{TAILTOSKIP} \n";
        print STDERR `$config{FILTERENRICH} -i $icShapeAllFile -o $icShapeFile -t $config{INPUTCOVERAGE} -T $config{TARGETHIT} -s $config{HEADTOSKIP} -e $config{TAILTOSKIP}`;
        if ( not $? ) {  print STDERR `touch $icShapeFile.done`; }
    }
    if ( not -e "$icShapeFile.done" ) { die "Abort! icSHAPE pipeline die of unsuccessful icSHAPE score filtering.\n"; }

    ## generateTrack
    #my $icShapeBedgraph = $outDir . "/icshape.bedgraph";
    #my $icShapeBw = $outDir . "/icshape.bw";
    #print STDERR "$config{SHAPE2BEDGRAPH} -i $icShapeAllFile -o $icShapeBedgraph\n";
    #print STDERR "$config{bedGraphToBigWig} $icShapeBedgraph $config{GENOMESIZE} $icShapeBw";
    #enrich2Bedgraph.pl LIB_NAI-LIB_DMSO.PolyA.invivo.valid.enrich > LIB_NAI-LIB_DMSO.PolyA.invivo.bedgraph
    #sort -k1,1 -k2,3n LIB_NAI-LIB_DMSO.PolyA.invivo.bedgraph -o LIB_NAI-LIB_DMSO.PolyA.invivo.sorted.bedgraph
    #uniqueTrack.pl LIB_NAI-LIB_DMSO.PolyA.invivo.sorted.bedgraph LIB_NAI-LIB_DMSO.PolyA.invivo.sorted.uniq.bedgraph
    #cut -f1-4 LIB_NAI-LIB_DMSO.PolyA.invivo.sorted.uniq.bedgraph | grep -v NULL > LIB_NAI-LIB_DMSO.PolyA.invivo.sim.bedgraph
    #bedGraphToBigWig LIB_NAI-LIB_DMSO.PolyA.invivo.sim.bedgraph ~/database/ensembl/current/mouse/dna/genome.sm.chr.size LIB_NAI-LIB_DMSO.PolyA.invivo.sim.bw

    print STDERR "icSHAPE pipeline finished successfully!\n...check $icShapeFile for output icSHAPE scores.\n\t", `date`;
    1;
}

sub init
{
    die $usage if ( $opt_h || ( not defined $opt_i ) || ( not defined $opt_t ) );
    $opt_V = 0 if ( not defined $opt_V );
    $opt_D = 0 if ( not defined $opt_D );

    ## looking for default config
    my $configFile = "";
    if ( defined $opt_c ) { $configFile = $opt_c; }
    else {
        $configFile = ".config";
        if ( not -e $configFile ) { if ( defined $ENV{"ICSHAPE"} ) { $configFile = $ENV{"ICSHAPE"} . "/.config"; } }
    }

    &config_pipeline ( $configFile );

    $config{input} = $opt_i;
    $config{target} = $opt_t;
    if ( defined $opt_o ) { $config{outDir} = $opt_o; }
    else { $config{outDir} = "out$$"; }

    if ( not -e $config{outDir} ) {
        print STDERR `mkdir $config{outDir}`;
        die "Cannot output to $config{outDir}!\n" if ( not -e $config{outDir} );
    }

    1;
}

sub config_pipeline
{
    my $configFile = shift;

    open ( CONFIG, $configFile ) or ( die "Cannot configure pipeline to run!\n" );
    while ( my $line = <CONFIG> ) {
        next if ( $line =~ /^#/ );
        chomp $line;

        my ( $key, $value ) = ( $line =~ /^(\S+)\s+(.+)$/ );
        $value =~ s/^\s+//; $value =~ s/\s+$//;
        $value =~ s/^"//; $value =~ s/"$//;
        $config{$key} = $value;
    }
    close CONFIG;

    1;
}
