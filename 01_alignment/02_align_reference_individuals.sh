#!/bin/bash
source "$(dirname "$0")/../config.sh"
# 공개 기준 개체(metadata/reference_individuals.tsv)의 짧은 읽기 → AcerK_1.0
# 설정은 원래 분석과 같다: bwa mem -M → samtools sort → GATK MarkDuplicates.
# 읽기는 <accession>_1.fastq.gz / _2.fastq.gz 로 $REFBAM_DIR/fastq 에 받아 둔다 (예: fasterq-dump).
set -eu
FQ=$REFBAM_DIR/fastq; mkdir -p $REFBAM_DIR/tmp
one(){ A=$1
  [ -s $REFBAM_DIR/${A}_dedup.bam ] && return 0
  bwa mem -t 14 -M -R "@RG\tID:$A\tSM:$A\tPL:ILLUMINA" $REF $FQ/${A}_1.fastq.gz $FQ/${A}_2.fastq.gz | \
    samtools sort -@ 14 -m 4G -T $REFBAM_DIR/tmp/$A -o $REFBAM_DIR/$A.sorted.bam -
  gatk --java-options "-Xmx20g" MarkDuplicates -I $REFBAM_DIR/$A.sorted.bam -O $REFBAM_DIR/${A}_dedup.bam \
    -M $REFBAM_DIR/${A}_dup_metrics.txt --CREATE_INDEX true --TMP_DIR $REFBAM_DIR/tmp --VALIDATION_STRINGENCY SILENT
  [ -s $REFBAM_DIR/${A}_dedup.bai ] && ln -sf ${A}_dedup.bai $REFBAM_DIR/${A}_dedup.bam.bai; }
export -f one; export REF FQ REFBAM_DIR
tail -n +2 $(dirname "$0")/../metadata/reference_individuals.tsv | cut -f1 | xargs -P 6 -I{} bash -c 'one {}'
