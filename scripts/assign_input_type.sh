#!/bin/bash

SAMPLE=$1

TYPE=$(grep $SAMPLE $GS_HOME/scripts/Blueprint_Epivar_public_results_release.index | cut -f27)

MONOCYTE=$GS_HOME/BluePrint/inputs/monocytes/monocytes.fastq.gz
GRANULOCYTE=$GS_HOME/BluePrint/inputs/granulocytes/granulocytes.fastq.gz
NAIVETCELL=$GS_HOME/BluePrint/inputs/naiveTcells/naiveTcells.fastq.gz
case "$TYPE" in

	*monocyte*)
	echo $MONOCYTE
;;
	*"T cell"*)
	echo $NAIVETCELL
;;
	*neutrophil*)
	echo $GRANULOCYTE
;;
	*)
	# for NA12878
	echo $GS_HOME/sorts/fastq/ENCFF002ECP.fastq.gz
;;
esac


