#!/bin/bash
source "$(dirname "$0")/../config.sh"
set -u
D=$KH_RESULTS/china_ref
CAL=$D/calib

mkdir -p $CAL/parts; cd $CAL
one() {
  C=$1
  awk -v c=$C '$1==c' $D/sites_sub.pos > parts/$C.pos
  [ -s parts/$C.pos ] || return 0
  bcftools mpileup -r $C -R parts/$C.pos -f $REF -b bams_cal.txt -a AD -q 20 -Q 20 -d 200 -I 2>/dev/null \
    | bcftools call -m -C alleles -T $D/alleles_sub.tsv.gz -Oz -o parts/$C.vcf.gz 2>/dev/null
  echo "  $C done $(date +%H:%M:%S)"
}
export -f one; export REF D CAL
cut -f1 $D/sites_sub.pos | uniq | xargs -P 8 -I{} bash -c 'one "$@"' _ {}
bcftools concat -Oz -o calib_geno.vcf.gz $(cut -f1 $D/sites_sub.pos | uniq | sed 's|^|parts/|;s|$|.vcf.gz|') 2>/dev/null
bcftools index -f calib_geno.vcf.gz
bcftools query -l calib_geno.vcf.gz | tr '\n' '\t' > calib_samples.txt
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' calib_geno.vcf.gz > calib_ad.tsv
echo "CALIB_DONE $(wc -l < calib_ad.tsv) sites"
