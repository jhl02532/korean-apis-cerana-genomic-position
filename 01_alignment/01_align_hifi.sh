#!/bin/bash
source "$(dirname "$0")/../config.sh"
# 군락 HiFi 읽기 → AcerK_1.0 (minimap2 map-hifi, samtools sort/index)
# RAW_DIR 에는 NABIC 파일명(<movie>.<sample>.<barcode>.fastq.gz)이 그대로 있다고 가정한다.
set -eu
mkdir -p $BAM_DIR $BAM_DIR/logs
one(){ S=$1; FQ=$(ls $RAW_DIR/*.$S.*.fastq.gz | head -1); bam=$BAM_DIR/$S.bam
  [ -s $bam.bai ] && return 0
  minimap2 -t 12 -ax map-hifi -R "@RG\tID:$S\tSM:$S" $REF $FQ 2>>$BAM_DIR/logs/$S.log | \
    samtools sort -@ 4 -m 2G -o $bam - 2>>$BAM_DIR/logs/$S.log && samtools index -@ 4 $bam; }
export -f one; export REF RAW_DIR BAM_DIR
xargs -P 4 -I{} bash -c 'one {}' < "$SAMPLES"
