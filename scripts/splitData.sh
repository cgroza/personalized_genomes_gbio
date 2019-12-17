#!/bin/bash
HAPLOID=$1
HAPLOID_FASTA=../$HAPLOID/$HAPLOID.fa

module load mugqic/bedtools
module load mugqic/vcftools

# for each reference peak, write overlapping $HAPLOID peak. Uniq the lines to record $HAPLOID peaks that overlap multiple peaks only once.
intersectBed -wa -a ${HAPLOID}_peaks.broadPeak -b reference_peaks.broadPeak | sort-bed - > intersected.broadPeak

# get a list of peaks that are in the indiviudualized referece, hg19 reference but not in both

subtractBed  -A -b reference_peaks.broadPeak -a ${HAPLOID}_peaks.broadPeak | sort-bed - > I-Ref.broadPeak
subtractBed  -A -a reference_peaks.broadPeak -b ${HAPLOID}_peaks.broadPeak | sort-bed - > Ref-I.broadPeak

intersectBed -c -a intersected.broadPeak -b SNPs.bed > intersected_SNPs_counts.bed
intersectBed -c -a intersected.broadPeak -b Indels.bed > intersected_indels_counts.bed

intersectBed -c -a I-Ref.broadPeak -b SNPs.bed > I-Ref_SNPs_counts.bed
intersectBed -c -a I-Ref.broadPeak -b Indels.bed > I-Ref_indels_counts.bed

intersectBed -c -a Ref-I.broadPeak -b SNPs.bed > Ref-I_SNPs_counts.bed
intersectBed -c -a Ref-I.broadPeak -b Indels.bed > Ref-I_indels_counts.bed

