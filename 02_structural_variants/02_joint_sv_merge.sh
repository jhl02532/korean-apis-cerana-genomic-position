#!/bin/bash
source "$(dirname "$0")/../config.sh"
# 19군락 .snf → Sniffles2 다중 시료 병합 → 유전형 행렬 (원래 분석 드라이버 스크립트의 SV 부분)
set -u
R=$KH_RESULTS; mkdir -p $R/recall
ls $R/snf/*.snf | grep -v tmp_ | grep -F -f <(sed 's#$#.snf#; s#^#/#' "$SAMPLES") > $R/recall/snf_list.tsv
echo "SNF $(wc -l < $R/recall/snf_list.tsv)"
sniffles --input $R/recall/snf_list.tsv --vcf $R/recall/merged_sv19.vcf --reference $REF --threads 32 > $R/recall/merge.log 2>&1
bcftools query -l $R/recall/merged_sv19.vcf | tr '\n' '\t' > $R/gt_header.txt; echo >> $R/gt_header.txt
bcftools query -f '%CHROM\t%POS\t%SVTYPE\t%SVLEN\t%SUPPORT\t%FILTER[\t%GT\t%DR\t%DV]\n' $R/recall/merged_sv19.vcf > $R/gt_matrix2.tsv
bcftools query -f '%CHROM\t%POS\t%SVTYPE\t%SVLEN\t%SUPPORT\t%FILTER[\t%GT]\n' $R/recall/merged_sv19.vcf > $R/gt_matrix.tsv
echo "SV_DONE records $(wc -l < $R/gt_matrix2.tsv) samples $(awk -F'\t' '{print NF}' $R/gt_header.txt)"
