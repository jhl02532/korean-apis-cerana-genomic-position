#!/bin/bash
source "$(dirname "$0")/../config.sh"
# 전장 SNP(97,609 거리 솎음 사이트)에서 군락 19 + 기준 개체 36 결합 유전형 판정 → 결합 PCA 입력
set -u
P=$KH_ROOT; N=$KH_ROOT; O=$KH_RESULTS/joint_pca_snp; mkdir -p $O/chunks $O/parts; cd $O
PG=$(dirname "$(command -v bcftools)")
$PG/bcftools --version | head -1 > bcftools_version.txt
# 사이트·대립유전자 (bim 위치 + filt.vcf.gz 의 REF/ALT)
awk '{print $1"\t"$4}' $KH_RESULTS/snpstruct/thin2000.bim > sites.pos
$PG/bcftools query -R sites.pos -f '%CHROM\t%POS\t%REF,%ALT\n' $KH_RESULTS/snpstruct/filt.vcf.gz | awk -F'\t' '$3!~/,.*,/' | sort -k1,1 -k2,2n -u > alleles.tsv
bgzip -f alleles.tsv; tabix -f -s1 -b2 -e2 alleles.tsv.gz
# BAM 목록: 군락 19, 중국 18 (calib_meta China*), 일본 18 (japan_anchor_selected)
sed "s#^#$BAM_DIR/#; s#\$#.bam#" "$SAMPLES" > bams.txt
awk -F'\t' 'NR>1 && $2 ~ /^China/ {print $1}' $KH_RESULTS/china_ref/calib/calib_meta.tsv > cn.ids
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="Unified_ID") c=i; next} {print $c}' $KH_RESULTS/japan_ref/japan_anchor_selected.tsv > jp.ids
grep -F -f cn.ids $KH_RESULTS/china_ref/calib/bams_cal.txt >> bams.txt
grep -F -f jp.ids $KH_RESULTS/japan_ref/bams_anchor.frozen.txt >> bams.txt
echo "BAMS $(wc -l < bams.txt) CN $(wc -l < cn.ids) JP $(wc -l < jp.ids) SITES $(zcat alleles.tsv.gz | wc -l)"
# 염색체 × 6 조각
zcat alleles.tsv.gz | cut -f1,2 > all.pos
awk '{print > "chunks/"$1".all"}' all.pos
for f in chunks/*.all; do c=$(basename $f .all); n=$(wc -l < $f); k=$(( (n+5)/6 )); split -l $k -d -a 1 $f chunks/${c}__; done
ls chunks | grep "__" > chunks.list
geno(){ K=$1; C=${K%%__*}
  $PG/bcftools mpileup -r $C -R chunks/$K -f $REF -b bams.txt -a AD,DP -q 20 -Q 20 -d 250 -I 2>/dev/null \
   | $PG/bcftools call -m -C alleles -T alleles.tsv.gz -Oz -o parts/$K.vcf.gz 2>/dev/null; echo "  $K $(date +%T)"; }
export -f geno; export REF PG
cat chunks.list | xargs -P 92 -I{} bash -c 'geno {}'
$PG/bcftools concat -Oz -o u.vcf.gz $(sed "s|^|parts/|;s|$|.vcf.gz|" chunks.list) 2>/dev/null
$PG/bcftools sort -Oz -o joint.vcf.gz u.vcf.gz 2>/dev/null
$PG/bcftools query -l joint.vcf.gz | tr '\n' '\t' > samples.txt
$PG/bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' joint.vcf.gz > ad.tsv
echo "GENO_DONE $(wc -l < ad.tsv) $(date +%T)"
