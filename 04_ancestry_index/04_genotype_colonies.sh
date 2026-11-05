#!/bin/bash
source "$(dirname "$0")/../config.sh"
set -u
# 19군락을 진단 부위에서 판정. 03_build_diagnostic_sites.py 의 307,737 부위 중 위치 순으로 8번째마다 (38,468)
D=$KH_RESULTS/china_ref
cd $D
awk 'NR%8==1' china_diag_alleles.tsv > sub_alleles.raw
cut -f1,2 sub_alleles.raw > sites_sub.pos
awk -v OFS='\t' '{print $1,$2,$3","$4}' sub_alleles.raw | bgzip -c > alleles_sub.tsv.gz; tabix -f -s1 -b2 -e2 alleles_sub.tsv.gz
mkdir -p parts
sed "s#^#$BAM_DIR/#; s#\$#.bam#" "$SAMPLES" > bams.txt
one() {
  C=$1
  awk -v c=$C '$1==c' sites_sub.pos > parts/$C.pos
  [ -s parts/$C.pos ] || return 0
  bcftools mpileup -r $C -R parts/$C.pos -f $REF -b bams.txt -a AD -q 20 -Q 20 -d 200 -I 2>/dev/null \
    | bcftools call -m -C alleles -T alleles_sub.tsv.gz -Oz -o parts/$C.vcf.gz 2>/dev/null
  echo "  $C done $(date +%H:%M:%S)"
}
export -f one; export REF D
cut -f1 sites_sub.pos | uniq | xargs -P 8 -I{} bash -c 'one "$@"' _ {}
bcftools concat -Oz -o china_diag_geno.vcf.gz $(cut -f1 sites_sub.pos | uniq | sed 's|^|parts/|;s|$|.vcf.gz|') 2>/dev/null
bcftools index -f china_diag_geno.vcf.gz
echo "GENO_DONE $(bcftools view -H china_diag_geno.vcf.gz | wc -l) sites"
bcftools query -l china_diag_geno.vcf.gz | tr '\n' '\t' > geno_samples.txt
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' china_diag_geno.vcf.gz > ad_matrix.tsv
