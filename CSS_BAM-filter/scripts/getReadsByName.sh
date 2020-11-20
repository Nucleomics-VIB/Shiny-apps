#!/bin/bash

# getReadsByName.sh: 
# parse a PB CCS BAM file
# extract reads from a list of read names
# save results to a second BAM file
#
# Stéphane Plaisance - VIB-NC-BITS Nov-19-2020 v1.0
#
# visit our Git: https://github.com/Nucleomics-VIB

if [ $# -ne 2 ]
then
	echo "# provide a bam file and list of read names to be kept!"
	exit 1
else
	inbam=$1; # Uploads/data.bam
	rlist=$2; # Uploads/reads_list.txt
	outbam=${inbam%.bam}_filtered.bam; # Uploads/data_filtered.bam
fi

# edit her teh path to your copy of PICARD

PICARD=/opt/biotools/picard

java -jar $PICARD/picard.jar FilterSamReads \
  --INPUT ${inbam} \
  --FILTER includeReadList \
  --READ_LIST_FILE ${rlist} \
  --OUTPUT ${outbam} \
  --VALIDATION_STRINGENCY SILENT
