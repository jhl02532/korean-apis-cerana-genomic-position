#!/bin/bash
source "$(dirname "$0")/../config.sh"
# 읽기 기반 유전자 보유: 19군락 HiFi BAM 에서 RefSeq 염색체 코딩 유전자의 CDS 평균 깊이 (mosdepth)
set -eu
P=$KH_ROOT
OUT=$KH_RESULTS/read_gene_presence; mkdir -p $OUT/md
GFF=$KH_REFSEQ_DIR/RS2026.coding.gff3   # 01_refseq_to_coding_gff.py 산출물
MD=mosdepth
# CDS → 유전자 ID (mRNA Parent → gene)
awk -F'\t' '$1~/^NC_0838/ && $3=="mRNA"{match($9,/ID=[^;]+/);m=substr($9,RSTART+3,RLENGTH-3);match($9,/Parent=[^;]+/);g=substr($9,RSTART+7,RLENGTH-7);print m"\t"g}' $GFF > $OUT/mrna2gene.tsv
# 유전자별로 CDS 구간 병합 (유전자 경계를 섞지 않도록)
awk -F'\t' 'NR==FNR{g[$1]=$2;next} $1~/^NC_0838/ && $3=="CDS"{match($9,/Parent=[^;]+/);m=substr($9,RSTART+7,RLENGTH-7);split(m,mm,",");if(mm[1] in g)print $1"\t"$4-1"\t"$5"\t"g[mm[1]]}' $OUT/mrna2gene.tsv $GFF \
 | sort -k4,4 -k1,1 -k2,2n | awk 'BEGIN{OFS="\t"} {if($4==pg && $1==pc && $2<=pe){if($3>pe)pe=$3} else {if(pg!="")print pc,ps,pe,pg; pc=$1;ps=$2;pe=$3;pg=$4}} END{print pc,ps,pe,pg}' \
 | sort -k1,1 -k2,2n > $OUT/cds_by_gene.bed
echo "genes: $(cut -f4 $OUT/cds_by_gene.bed | sort -u | wc -l)  intervals: $(wc -l < $OUT/cds_by_gene.bed)"
sed "s#^#$BAM_DIR/#; s#\$#.bam#" "$SAMPLES" > $OUT/bams.txt
run(){ b=$1; s=$(basename $b .bam); [ -s $OUT/md/$s.regions.bed.gz ] && return 0
  $MD -t 4 -n --fast-mode --by $OUT/cds_by_gene.bed $OUT/md/$s $b && echo "$s done $(date +%T)"; }
export -f run; export MD OUT
cat $OUT/bams.txt | xargs -P 8 -I{} bash -c 'run {}'
echo READ_GENE_DONE
