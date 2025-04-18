#!/bin/bash

#path to your bigWig files
bw_dir="/Users/FV/Scienze Mediche Dropbox/Pers.immuno/fe.vallone/other/metabolic adaptation in chronic lymphocytic leukemia (AF, SS)/cut&tag/analysis/bw"
cd "$bw_dir"

#path to chr sizes file
chrom_sizes="/Users/FV/Scienze Mediche Dropbox/Pers.immuno/fe.vallone/other/metabolic adaptation in chronic lymphocytic leukemia (AF, SS)/cut&tag/analysis/bw/fixed_chrom.sizes.hg38"

#create output directory for merged files
mkdir -p merged

#loop through unique condition+histone pairs
for combo in $(ls *.bw | sed -E 's/_rep[0-9]+\.bw//' | sort | uniq); do
    base="$combo"
    replicates=$(ls ${base}_rep*.bw 2>/dev/null)

    #count number of replicates
    rep_count=$(echo "$replicates" | wc -l)

    if [[ $rep_count -eq 3 ]]; then
        echo "Merging $rep_count replicates for: $base"

        #merge replicates into bedGraph
        bigWigMerge $replicates "${base}_merged.bedGraph"

        #sort bedGraph
        sort -k1,1 -k2,2n "${base}_merged.bedGraph" > "${base}_merged.sorted.bedGraph"

        #filter out contigs not listed in chrom.sizes
        grep -v -E "^(GL|KI|chrUn|HSCHR|JH|MT|EBV|NC|NT)" "${base}_merged.sorted.bedGraph" > "${base}_merged.sorted.filtered.bedGraph"

        #check if filtered file is not empty
        if [[ -s "${base}_merged.sorted.filtered.bedGraph" ]]; then
            #convert to bigWig and move to merged/ directory
            bedGraphToBigWig "${base}_merged.sorted.filtered.bedGraph" "$chrom_sizes" "merged/${base}_merged.bw"

            #only clean up if .bw was successfully created
            if [[ -f "merged/${base}_merged.bw" ]]; then
                rm "${base}_merged.bedGraph" "${base}_merged.sorted.bedGraph" "${base}_merged.sorted.filtered.bedGraph"
                echo "Merged bigWig saved: merged/${base}_merged.bw"
            else
                echo "Conversion to bigWig failed for $base : .bw file not created"
            fi
        else
            echo "Filtered bedGraph is empty for $base : skipping conversion"
        fi
    else
        echo "Skipping $base : expected 3 replicates, found $rep_count"
    fi
done