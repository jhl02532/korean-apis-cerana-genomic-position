#!/bin/bash
source "$(dirname "$0")/../config.sh"
# 일본 기준 개체 18개를 Korea–China 사이트(china_ref sites_sub)에서 판정 — 06_genotype_chinese_reference.sh 와 같은 설정
set -u
P=$KH_RESULTS
D=$P/china_ref; O=$P/japan_ref/jp_on_kc; mkdir -p $O/parts; cd $O
PG=$(dirname "$(command -v bcftools)")
cut -f5 $P/japan_ref/japan_anchor_selected.tsv | tail -n +2 > bams_jp.txt
awk '{print > "parts/"$1".all"}' $D/sites_sub.pos
for f in parts/*.all; do c=$(basename $f .all); n=$(wc -l < $f); k=$(( (n+5)/6 )); split -l $k -d -a 1 $f parts/${c}__; done
ls parts | grep "__" > chunks.list
g(){ K=$1; C=${K%%__*}; $PG/bcftools mpileup -r $C -R parts/$K -f $REF -b bams_jp.txt -a AD -q 20 -Q 20 -d 200 -I 2>/dev/null \
   | $PG/bcftools call -m -C alleles -T $D/alleles_sub.tsv.gz -Oz -o parts/$K.vcf.gz 2>/dev/null; }
export -f g; export REF PG D
cat chunks.list | xargs -P 96 -I{} bash -c 'g {}'
$PG/bcftools concat -Oz -o u.vcf.gz $(sed "s|^|parts/|;s|$|.vcf.gz|" chunks.list) 2>/dev/null; $PG/bcftools sort -Oz -o jp_kc.vcf.gz u.vcf.gz 2>/dev/null
$PG/bcftools query -l jp_kc.vcf.gz | tr '\n' '\t' > samples.txt
$PG/bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' jp_kc.vcf.gz > ad.tsv
echo "JP_KC_DONE $(wc -l < ad.tsv)"
