#!/bin/bash
source "$(dirname "$0")/../config.sh"
# 기준 개체 36(중국 18, 일본 18) BAM 목록과 메타데이터를 각 단계가 읽는 위치에 만든다.
set -eu
M=$(dirname "$0")/../metadata/reference_individuals.tsv
CAL=$KH_RESULTS/china_ref/calib; JP=$KH_RESULTS/japan_ref; mkdir -p $CAL $JP
awk -F'\t' 'NR>1 && $2=="China" {print "'$REFBAM_DIR'/"$1"_dedup.bam"}' $M > $CAL/bams_cal.txt
awk -F'\t' 'BEGIN{print "sample\tgroup"} NR>1 && $2=="China" {print $1"\tChina:"$4}' $M > $CAL/calib_meta.tsv
awk -F'\t' 'BEGIN{print "Unified_ID\tRegion\tgrp\tgb\tbam"} NR>1 && $2=="Japan" {print $1"\t"$4"\t"$3"\tNA\t'$REFBAM_DIR'/"$1"_dedup.bam"}' $M > $JP/japan_anchor_selected.tsv
cp $CAL/bams_cal.txt $JP/bams_anchor.frozen.txt; cut -f5 $JP/japan_anchor_selected.tsv | tail -n +2 >> $JP/bams_anchor.frozen.txt
echo "China $(wc -l < $CAL/bams_cal.txt)  Japan $(($(wc -l < $JP/japan_anchor_selected.tsv)-1))"
