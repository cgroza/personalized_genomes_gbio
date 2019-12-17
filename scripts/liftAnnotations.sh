#!/bin/bash

module load mugqic/bedtools mugqic/vcftools mugqic/ucsc

scripts=$GS_HOME/scripts
haploid=$1
g=$2
varcalls=$3
chipseq_path=$4
sample_name=$5
pwd=$6

echo liftAnnotation Current Directory for $haploid:  $(pwd)
# sort variations by SNPs, indels, and liftOver to haploid
for varcall in $varcalls/*.vcf.gz
do
V=$(basename $varcall)
echo Doing $varcall
vcftools --indv $g --gzvcf $varcall --keep-only-indels --stdout --recode | \
	perl $scripts/haploid_variations.pl $haploid | vcf2bed | liftOver -tab -bedPlus=5 stdin chains/$haploid.chain  bed_comparison_$haploid/$V.Indels.bed  bed_comparison_$haploid/$V.unlifted_Indels.txt

vcftools --indv $g --gzvcf $varcall --remove-indels --stdout --recode | \
	perl $scripts/haploid_variations.pl $haploid | vcf2bed | liftOver -tab -bedPlus=5 stdin chains/$haploid.chain bed_comparison_$haploid/$V.SNPs.bed bed_comparison_$haploid/$V.unlifted_SNPs.txt

done

echo Merging
# merge variations from all chromosomes and sort the final bed
cat bed_comparison_$haploid/*.Indels.bed | sortBed  > bed_comparison_$haploid/Indels.bed
cat bed_comparison_$haploid/*.SNPs.bed  | sortBed >  bed_comparison_$haploid/SNPs.bed

echo Cleaning
# remove intermediate temp files
rm bed_comparison_$haploid/*.Indels.bed bed_comparison_$haploid/*.SNPs.bed #bed_comparison_$haploid/*unlifted*

echo Removing sex chromosomes
# remove sex chromosmes from peak calls and place them in the comparison directory for the haploid
cat peak_call/$haploid/${haploid}_peaks.broadPeak | $scripts/removeSexChr.sh > bed_comparison_$haploid/${haploid}_peaks.broadPeak

echo Lifting reference peaks
# liftOver the reference peaks and place them in the comparison directory for the haploid
liftOver -tab -bedPlus=3 peak_call/reference/reference_peaks.broadPeak chains/$haploid.chain stdout bed_comparison_$haploid/unlifted_reference_peaks.txt | $scripts/removeSexChr.sh > bed_comparison_$haploid/reference_peaks.broadPeak

#liftOver -tab -bedPlus=4 reference_read_counts.bed chains/$haploid.chain stdout /dev/null | $scripts/removeSexChr.sh > bed_comparison_$haploid/reference_read_counts.bed

echo Lifting repeats
# lift repeat maksker annotations
liftOver -gff $pwd/annotations/Alus.gff chains/$haploid.chain bed_comparison_$haploid/Alus.gff bed_comparison_$haploid/unlifted_Alus.txt
liftOver -gff $pwd/annotations/LINEs.gff chains/$haploid.chain bed_comparison_$haploid/LINEs.gff bed_comparison_$haploid/unlifted_LINEs.txt
liftOver -gff $pwd/annotations/hg19.repeats.gff chains/$haploid.chain bed_comparison_$haploid/Repeats.gff bed_comparison_$haploid/unlifted_Repeats.txt

if [ -f ALU.final_comp.vcf ] && [ -f LINE1.final_comp.vcf ]
then
# lift ALU, L1 insertions
echo Lifting ALU, LINE1 insertions
$scripts/renameChrom.pl ALU.final_comp.vcf | vcf2bed | liftOver -tab -bedPlus=5 stdin chains/$haploid.chain stdout bed_comparison_$haploid/unlifted_AluInsertions.txt | sortBed > bed_comparison_$haploid/AluInsertions.bed
$scripts/renameChrom.pl LINE1.final_comp.vcf | vcf2bed | liftOver -tab -bedPlus=5 stdin chains/$haploid.chain stdout bed_comparison_$haploid/unlifted_LINEInsertions.txt | sortBed > bed_comparison_$haploid/LINEInsertions.bed

else
liftOver -tab -bedPlus=4 $g.alu.bed chains/$haploid.chain stdout bed_comparison_$haploid/unlifted_AluInsertions.txt | \
 shiftBed  -s -282 -i /dev/stdin -g bed_comparison_$haploid/Haploid.chrom.sizes | sortBed > bed_comparison_$haploid/AluInsertions.bed

liftOver -tab -bedPlus=4 $g.line1.bed chains/$haploid.chain stdout bed_comparison_$haploid/unlifted_LINEInsertions.txt | \
 shiftBed  -s -6019 -i /dev/stdin -g bed_comparison_$haploid/Haploid.chrom.sizes | sortBed > bed_comparison_$haploid/LINEInsertions.bed

fi
echo Lifting reference coverage
liftOver reference_coverage.bed chains/$haploid.chain stdout bed_comparison_$haploid/unlifted_Hg19Coverage.bed | sort-bed --max-mem 20G - > bed_comparison_$haploid/Hg19Coverage.bed
