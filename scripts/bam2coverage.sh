#!/bin/bash

BAM=$1
CHROM_SIZES=$2

#module load mugqic/igvtools
module load mugqic/bedtools

#igvtools count $BAM stdout $CHROM_SIZES | wig2bed | sortBed
bam2bed -d < $BAM | awk -v OFS='\t' '{ print $1,$2,$3 }' | sort-bed -
