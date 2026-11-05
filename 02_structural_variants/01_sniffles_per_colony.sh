#!/bin/bash
source "$(dirname "$0")/../config.sh"
# 군락별 Sniffles2 .snf (다중 시료 병합은 02_joint_sv_merge.sh)
OUT=$KH_RESULTS/snf
SNIF=sniffles
mkdir -p $OUT $OUT/logs
one() {
  B=$1; S=$(basename "$B" .bam)
  [ -s "$OUT/$S.snf" ] && return 0
  $SNIF --input "$B" --snf "$OUT/$S.snf" --reference "$REF" --sample-id "$S" \
        --threads 4 --minsvlen 50 --mapq 20 --output-rnames \
        > "$OUT/logs/$S.sniffles.log" 2>&1
}
export -f one; export OUT REF SNIF
sed "s#^#$BAM_DIR/#; s#\$#.bam#" "$SAMPLES" | xargs -P 4 -I{} bash -c 'one "$@"' _ {}
echo SNF_DONE
