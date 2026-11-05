#!/bin/bash
source "$(dirname "$0")/../config.sh"
D=$KH_RESULTS/china_ref
mkdir -p $D && cd $D
get() {  # $1=accession $2=asmname $3=label
  A=$1; N=$2; L=$3
  P=$(echo $A | sed 's/_//;s/\(...\)\(...\)\(...\)\(...\).*/\1\/\2\/\3\/\4/' | sed 's/^GCA//')
  U="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/${P#/}/${A}_${N}/${A}_${N}_genomic.fna.gz"
  echo "[$L] $U"
  curl -sS -m 900 -o ${L}.fna.gz "$U" && echo "  ok $(du -h ${L}.fna.gz | cut -f1)" || echo "  FAIL $L"
}
# 표지는 BioSample 채집지 (Methods §6)
get GCA_011100585.1 ASM1110058v1 CN_Hainan
get GCA_020141525.2 ASM2014152v2 CN_Sichuan
get GCA_029531695.1 ASM2953169v1 CN_Hubei
get GCA_039655025.1 ASM3965502v1 CN_Qinghai
echo DL_DONE
