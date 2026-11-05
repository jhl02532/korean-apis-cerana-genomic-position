#!/bin/bash
source "$(dirname "$0")/../config.sh"
D=$KH_RESULTS/china_ref

cd $D
one() {
  L=$1
  [ -s $L.vcf ] && return 0
  minimap2 -cx asm10 --cs -t 12 $REF $L.fna.gz 2> $L.mm2.log \
    | sort -k6,6 -k8,8n \
    | k8 "$(command -v paftools.js)" call -L 10000 -l 1000 -f $REF - \
    > $L.vcf 2> $L.call.log
  echo "$L $(grep -vc '^#' $L.vcf) variants"
}
export -f one; export REF D
printf '%s\n' CN_Hainan CN_Sichuan CN_Hubei CN_Qinghai | xargs -P 2 -I{} bash -c 'one "$@"' _ {}
echo CALL_DONE
