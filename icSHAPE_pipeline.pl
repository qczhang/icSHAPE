#! /usr/bin/perl
# wrapper of icSHAPE pipeline
# copy right qiangfeng.zhang@gmail.com
# history: 0.01 
#   date: 01/06/2015

use strict;
use warnings;
use File::Basename;
use Getopt::Std;

my $_debug = 1;
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
    foreach my $key ( keys %config ) { print $key, "\t", $config{$key}, "\n"; }

    my $outDir = $config{outDir};
    my @inputSeqFile = split ( /:/, $config{input} ); my @targetSeqFiles = split ( /:/, $config{target} );
    my @inputSignalFile = (); my @targetSignalFile = ();
    foreach my $seqFile ( @inputSeqFile, @targetSeqFiles ) {
        my ($fileName, $fileDir, $fileSuffix) = fileparse ( $seqFile, qr/\.[^.]*/ ); 
        my $seqCollapsed = $outDir . "/" . $fileName . ".rmdup.fastq";
        my $seqDatFasta = $outDir . "/" . $fileName . ".fa";
        print STDERR "$config{COLLAPSEBIN} -U $seqFile -o $seqCollapsed -f $seqDatFasta\n";
#        print STDERR `$config{COLLAPSEBIN} -U $seqFile -o $seqCollapsed -f $seqDatFasta`;
        if ( not $? ) {  print STDERR `touch $seqCollapsed.done`;  }

        my $seqTrimmed = $outDir . "/" . $fileName . ".trimmed.fastq";
        print STDERR "$config{TRIMMER} -U $seqCollapsed -o $seqTrimmed -l $config{LEADINGTRIM} -t $config{TAILINGTRIM} -c $config{FASTQCODING} -a $config{ADAPTER}\n";
#        print STDERR `$config{TRIMMER} -U $seqCollapsed -o $seqTrimmed -l $config{LEADINGTRIM} -t $config{TAILINGTRIM} -c $config{FASTQCODING} -a $config{ADAPTER}`;
        if ( not $? ) {  print STDERR `touch $seqTrimmed.done; /bin/rm $seqCollapsed.done`;  }

        my $samFile = $outDir . "/" . $fileName . ".sam";
        my $alignOptions = ""; 
        $alignOptions .= "-" . $config{FASTQCODING} if ( defined $config{FASTQCODING} );
        $alignOptions .= " " . $config{MAPPINGOPTIONS} if ( defined $config{MAPPINGOPTIONS} );
        print STDERR "$config{ALIGNER} -U $seqTrimmed -o $samFile -x $config{MAPPINGREF} $alignOptions\n";
#        print STDERR `$config{ALIGNER} -U $seqTrimmed -o $samFile -x $config{MAPPINGREF} $alignOptions`;
        if ( not $? ) {  print STDERR `touch $samFile.done; /bin/rm $seqTrimmed.done`;  }

        my $rpkmFile = $outDir . "/" . $fileName . ".rpkm";
        print STDERR "$config{ESTIMATERPKM} -i $samFile -o $rpkmFile\n";
#        print STDERR `$config{ESTIMATERPKM} -i $samFile -o $rpkmFile`;
        if ( not $? ) {  print STDERR `touch $rpkmFile.done; /bin/rm $samFile.done`;  }

        # my $seqSumFile = $outDir . "/" . $fileName . ".seq.info";
        # print STDERR "$config{seqSummary} -i $samFile -o $seqSumFile -r $rpkmFile\n";
        # print STDERR `$config{seqSummary} -i $samFile -o $seqSumFile -r $rpkmFile`;

        my $rtFile = $outDir . "/" . $fileName . ".rt";
        print STDERR "$config{CALCRT} -i $samFile -o $rtFile -r $rpkmFile\n";
#        print STDERR `$config{CALCRT} -i $samFile -o $rtFile -r $rpkmFile`;
        if ( not $? ) {  print STDERR `touch $rtFile.done; /bin/rm $rpkmFile.done`;  }
        push ( @inputSignalFile, $rtFile );
        push ( @targetSignalFile, $rtFile );
    }

    my $combinedInputSignalFile = $outDir . "/background.rt";
    my $normalizedInputSignalFile = $outDir . "/background.normalized.rt";
    if ( scalar @inputSignalFile > 1 ) {
        my $allInputSignalFile = "";
        foreach my $signalFile ( @inputSignalFile ) { $allInputSignalFile .= $signalFile . ":"; }
        $allInputSignalFile =~ s/:$//;
        print STDERR "$config{COMBINEBIN} -a $allInputSignalFile -o $combinedInputSignalFile";
#        print STDERR `$config{COMBINEBIN} -a $allInputSignalFile -o $combinedInputSignalFile`;
    }
    else { print STDERR `ln -s $inputSignalFile[0] $combinedInputSignalFile` }
    print STDERR "$config{NORMALIZEBIN} -i $combinedInputSignalFile -o $normalizedInputSignalFile -m $config{METHOD} -d $config{HEADTOSKIP} -l $config{TAILTOSKIP}\n";
#    print STDERR `$config{NORMALIZEBIN} -i $combinedInputSignalFile -o $normalizedInputSignalFile -m $config{METHOD} -d $config{HEADTOSKIP} -l $config{TAILTOSKIP}`;

    my $combinedTargetSignalFile = $outDir . "/target.rt";
    my $normalizedTargetSignalFile = $outDir . "/target.normalized.rt";
    if ( scalar @targetSignalFile > 1 ) {
        my $allTargetSignalFile = "";
        foreach my $signalFile ( @targetSignalFile ) { $allTargetSignalFile .= $signalFile . ":"; }
        $allTargetSignalFile =~ s/:$//;
        print STDERR "$config{COMBINEBIN} -a $allTargetSignalFile -o $combinedTargetSignalFile";
        #       print STDERR `$config{COMBINEBIN} -a $allTargetSignalFile -o $combinedTargetSignalFile`;
    }
    else { print STDERR `ln -s $targetSignalFile[0] $combinedTargetSignalFile` }
    print STDERR "$config{NORMALIZEBIN} -i $combinedTargetSignalFile -o $normalizedTargetSignalFile -m $config{METHOD} -d $config{HEADTOSKIP} -l $config{TAILTOSKIP}\n";
#    print STDERR `$config{NORMALIZEBIN} -i $combinedTargetSignalFile -o $normalizedTargetSignalFile -m $config{METHOD} -d $config{HEADTOSKIP} -l $config{TAILTOSKIP}`;

    ## calculate enrichment scores and filter valid ones
    my $enrichFile = $outDir . "/icshape.out"; my $enrichAllFile = $outDir . "/icshape.tmp.out";
    print STDERR "$config{CALCENRICHBIN} -f $normalizedTargetSignalFile -b $normalizedInputSignalFile -o $enrichAllFile -w $config{WINSOR} -x $config{DIVFACTOR} -y $config{SUBFACTOR}\n";
#    print STDERR `$config{CALCENRICHBIN} -f $normalizedTargetSignalFile -b $normalizedInputSignalFile -o $enrichAllFile -w $config{WINSOR} -x $config{DIVFACTOR} -y $config{SUBFACTOR}`;
    print STDERR "$config{FILTERENRICH} -a $enrichAllFile -o $enrichFile\n";
#    print STDERR `$config{FILTERENRICH} -a $enrichAllFile -o $enrichFile`;

    ## generateTrack

    1;
}

sub init
{
    die $usage if ( $opt_h || ( not defined $opt_i ) || ( not defined $opt_t ) );
    $opt_V = 0 if ( not defined $opt_V );
    $opt_D = 0 if ( not defined $opt_D );

    ## looking for default config
    my $configFile = ".config";
    if ( defined $opt_c ) { $configFile = $opt_c; }
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

    if ( not defined $configFile ) {
        ## search for configuration file
        if ( -e ".config" ) { $configFile = ".config"; }
        elsif ( ( defined $ENV{"ICSHAPECONFIG"} ) and ( -e $ENV{"ICSHAPECONFIG"} ) ) { $configFile = $ENV{"ICSHAPECONFIG"}; }
    }

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
