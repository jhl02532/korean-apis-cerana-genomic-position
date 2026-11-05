#!/bin/bash
source "$(dirname "$0")/../config.sh"
# 19군락 합동 SNP 호출(16 염색체, 5 Mb 조각 병렬) → 필터 → MAF 0.05 → 2 kb 거리 솎기
# (원래 분석 드라이버 스크립트의 SNP 부분)
set -u
R=$KH_RESULTS; mkdir -p $R/recall/vcf $R/snpstruct
PG=$(dirname "$(command -v bcftools)")
sed "s#^#$BAM_DIR/#; s#\$#.bam#" "$SAMPLES" > $R/recall/bamlist.txt; echo "BAM $(wc -l < $R/recall/bamlist.txt)"
awk '$1 ~ /^NC_0838(5[2-9]|6[0-7])\.1$/ {for(s=1;s<=$2;s+=5000000){e=s+4999999; if(e>$2)e=$2; print $1":"s"-"e}}' $REF.fai > $R/recall/regions.txt
echo "REGIONS $(wc -l < $R/recall/regions.txt)"
doreg(){ G=$1; T=$(echo $G | tr ':-' '__')
  $PG/bcftools mpileup -f $REF -b $R/recall/bamlist.txt -r $G -q 20 -Q 20 -d 100 -I -Ou 2>/dev/null | $PG/bcftools call -m -v -Oz -o $R/recall/vcf/$T.vcf.gz 2>/dev/null; }
export -f doreg; export REF R PG
cat $R/recall/regions.txt | xargs -P ${THREADS:-32} -I{} bash -c 'doreg {}'
ls $R/recall/vcf/*.vcf.gz | xargs -P 16 -I{} $PG/bcftools index -f -t {}
ls $R/recall/vcf/*.vcf.gz | sort -V > $R/recall/vcf_list.txt
$PG/bcftools concat -a -f $R/recall/vcf_list.txt --threads 16 -Oz -o $R/snpstruct/all.vcf.gz && $PG/bcftools index -f -t $R/snpstruct/all.vcf.gz
$PG/bcftools view -m2 -M2 -v snps -i 'QUAL>=30' -Ou $R/snpstruct/all.vcf.gz | $PG/bcftools view -e 'F_MISSING>0.05' -Oz -o $R/snpstruct/filt.vcf.gz && $PG/bcftools index -f -t $R/snpstruct/filt.vcf.gz
echo "SNP_FILT $($PG/bcftools view -H $R/snpstruct/filt.vcf.gz | wc -l)"
plink2 --vcf $R/snpstruct/filt.vcf.gz --allow-extra-chr --set-missing-var-ids @:# --maf 0.05 --make-bed --out $R/snpstruct/snp19 > $R/snpstruct/plink1.log 2>&1
echo "SNP_MAF $(wc -l < $R/snpstruct/snp19.bim)"
plink2 --bfile $R/snpstruct/snp19 --allow-extra-chr --bp-space 2000 --make-bed --out $R/snpstruct/thin2000 > $R/snpstruct/thin2000.log 2>&1
echo "SNP_THIN $(wc -l < $R/snpstruct/thin2000.bim)"
