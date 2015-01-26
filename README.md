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
    This pipeline calculates RNA base reactivity scores from icSHAPE experiments. 


