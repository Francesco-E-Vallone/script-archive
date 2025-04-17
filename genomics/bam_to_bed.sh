#!/bin/bash

#dir
BASE_DIR="/Users/FV/Scienze Mediche Dropbox/Pers.immuno/fe.vallone/other/metabolic adaptation in chronic lymphocytic leukemia (AF, SS)/cut&tag/analysis"

#recursively find all .bam files
find "$BASE_DIR" -type f -name "*.bam" | while read -r bam_file; do
    #output .bed file in the same directory as the .bam
    bed_file="${bam_file%.bam}.bed"

    #conversion
    echo "Converting: $bam_file → $bed_file"
    bedtools bamtobed -i "$bam_file" > "$bed_file"
done