README
=======

NAME
##------
icSHAPE pipeline
    pipeline that calculates RNA base reactivity scores from icSHAPE experiments 

SYNOPSIS
##------
perl icSHAPE_pipeline.pl -i input_fastq_file1:input_fastq_file2:... -t target_fastq_file1:target_fastq_file2:... -o output_directory -c configuration_file

INSTALLATION
##------
1, check out the pipeline from https://github.com/qczhang/icSHAPE and put into any directory (referred as installation_directory in the following) you want.
2, define enviroment variables by add the following to your bash profiles:
    export ICSHAPE=installation_directory
    export BOWTIELIB=bowtie_index_directory
   Alternatively, you can edit the .config file so the pipeline can find scripts and programs used in the pipeline without defining the above enviroment variables
3, edit .config file for additional configurations that specify how the pipeline can find the correct scripts and programs, e.g., java, bowtie, etc.

PREREQUISITES 
##------
1, java
2, bowtie
3, ucsc_tools if generating ucsc tracks

DESCRIPTION
##------
    This pipeline calculates RNA base reactivity scores from icSHAPE experiments. It performs the following steps:
1, PCR duplicates removal
    scripts and programs used: scripts/readCollapse.pl bin/readCollapse
    related options in configuration file: 
    COLLAPSEBIN     $ICSHAPE/scripts/readCollapse.pl 
    COLLAPSEFASTA   collapse.fa
2, adapter trimming and possible barcode removal
    scripts and programs used: scripts/trimming.pl bin/simpleTrim bin/trimmomatic-0.30.jar
    related options in configuration file: 
    FASTQCODING     phred33
    JAVABIN         /usr/java/latest/bin/java
    TRIMMER         $ICSHAPE/scripts/trimming.pl
    TRIMLOG         trimming.log
    ADAPTER         $ICSHAPE/data/TruSeq2-PE.fa
    LEADINGTRIM     13
    TAILINGTRIM     0
    TRIMMINLEN      25
3, reads mapping
    scripts and programs used: custom installation of bowtie 
    related options in configuration file: 
    FASTQCODING     phred33
    ALIGNER         /srv/gs1/software/bowtie/2.0.5/bowtie2
    MAPPINGREF      $BOWTIELIB/mouse/ensembl.transcriptome
    MAPPINGOPTIONS  "--non-deterministic --time"
4, transcript abundance estimation
    scripts and programs used: scripts/estimateRPKM.pl
    related options in configuration file: 
    ESTIMATERPKM    $ICSHAPE/scripts/estimateRPKM.pl
5, RT stop calculation
    scripts and programs used: scripts/calcRT.pl
    related options in configuration file: 
    CALCRT          $ICSHAPE/scripts/calcRT.pl
    MINLOAD         5
6, replicate combining
    scripts and programs used: scripts/combineRTreplicates.pl
    related options in configuration file: 
    COMBINEBIN      $ICSHAPE/scripts/combineRTreplicates.pl
7, normalization
    scripts and programs used: scripts/normalizeRTfile.pl
    NORMALIZEBIN    $ICSHAPE/scripts/normalizeRTfile.pl
    HEADTOSKIP      32
    TAILTOSKIP      32
    METHOD          mean:vigintile2
8, calculation of reactivity score
    scripts and programs used: scripts/calcEnrich.pl
    CALCENRICHBIN   $ICSHAPE/scripts/calcEnrich.pl
    WINSOR          factor5:scaling1
    DIVFACTOR       10
    SUBFACTOR       0.25
9, filter and seletion valid reactivity score
    scripts and programs used: scripts/filterEnrich.pl
    FILTERENRICH    $ICSHAPE/scripts/filterEnrich.pl
    INPUTCOVERAGE   200
    TARGETHIT       2
    HEADTOSKIP      5
    TAILTOSKIP      30
10, generate UCSC track (not implemented)

Possible additional steps:
0, split library by barcodes
    scripts and programs to use: scripts/splitFastq.pl
    SYNOPSIS
    paired-ends reads:
    perl scripts/splitFastq.pl -1 fastq_PE_reads_1 -2 fastq_PE_reads_2 -l barcode1:lib_name1::barcode2:lib_name2 -b barcode_position:barcode_length
    single-ends reads:
    perl scripts/splitFastq.pl -U fastq_SE_reads -l barcode1:lib_name1::barcode2:lib_name2 -b barcode_position:barcode_length





