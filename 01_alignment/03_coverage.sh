#!/bin/bash
source "$(dirname "$0")/../config.sh"
# 군락별 16 염색체 평균 깊이·참조 피복률 (samtools coverage) → $KH_RESULTS/qc/all_coverage.tsv (Table 1)
set -eu
mkdir -p $KH_RESULTS/qc
do_one(){ S=$1
  samtools coverage $BAM_DIR/$S.bam | awk -v s="$S" 'NR>1 && $1 ~ /^NC_0838/ {d+=$7*$3; c+=$6*$3; len+=$3} END{printf "%s\t%.1f\t%.2f\n", s, d/len, c/len}' > $KH_RESULTS/qc/$S.cov; }
export -f do_one; export BAM_DIR KH_RESULTS
xargs -P 6 -I{} bash -c 'do_one {}' < "$SAMPLES"
cat $KH_RESULTS/qc/*.cov | sort > $KH_RESULTS/qc/all_coverage.tsv
