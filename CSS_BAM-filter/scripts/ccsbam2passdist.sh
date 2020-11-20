#!/bin/bash

# ccsbam2passdist.sh: 
# parse a PB CCS BAM file
# extract molecule ID, read length, ccs pass np:i:# ccs accuracy rq:i:#
# save results to a text table (TSV) for stats in R
#
# Stéphane Plaisance - VIB-NC-BITS Nov-19-2020 v1.0
#
# visit our Git: https://github.com/Nucleomics-VIB

# required:
# Samtools for parsing bam
# GNU awk for parsing sam data

if [ -z "${1}" ]
then
	echo "# provide a bam file to be parsed!"
	exit 1
else
	inbam=$1
	out=$2
fi

# check if requirements are present
PATH=$PATH:/opt/biotools/bin
$( hash samtools 2>/dev/null ) || ( echo "# samtools not found in PATH"; exit 1 )
$( hash awk 2>/dev/null ) || ( echo "# awk not found in PATH"; exit 1 )

samtools view ${inbam} | \
awk 'BEGIN{FS="\t"; OFS=","; print "readID","MovieName","ZMW","npass","Accuracy","len"}
	{split($1,hd,"/");
	if( $14 ~ /np:i/ ){
		split($14,np,":");
		split($15,rq,":");
		print $1,hd[1],hd[2],np[3],rq[3],length($10) 
		}
	}' > ${out}
