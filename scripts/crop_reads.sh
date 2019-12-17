#!/bin/bash
READS=$1
CROPPED=$2
CROP_SIZE=$3
module load mugqic/trimmomatic

java -Xmx20G -jar $TRIMMOMATIC_JAR SE $READS $CROPPED CROP:$CROP_SIZE

