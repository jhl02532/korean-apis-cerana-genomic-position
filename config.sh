# 모든 셸 스크립트가 source 하는 경로 설정. 자기 환경에 맞게 고친다.
# Python 스크립트는 같은 변수를 환경변수로 읽으므로, 실행 전에 `source config.sh` 를 한다.
export KH_ROOT=${KH_ROOT:-$HOME/kh19}                    # 작업 폴더
export KH_RESULTS=${KH_RESULTS:-$KH_ROOT/03_results}      # 모든 산출물
export KH_REFSEQ_DIR=${KH_REFSEQ_DIR:-$KH_ROOT/refseq}    # RefSeq GFF 와 변환본
export REF=${REF:-$KH_ROOT/ref/GCF_029169275.1_AcerK_1.0_genomic.fna}   # AcerK_1.0 (samtools faidx 필요)
export RAW_DIR=${RAW_DIR:-$KH_ROOT/raw_hifi}              # NABIC 에서 받은 HiFi fastq.gz
export BAM_DIR=${BAM_DIR:-$KH_ROOT/bam}                   # 군락 HiFi 정렬 BAM (<sample>.bam)
export REFBAM_DIR=${REFBAM_DIR:-$KH_ROOT/ref_bam}         # 기준 개체 정렬 BAM (<accession>_dedup.bam)
export SAMPLES=${SAMPLES:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/metadata/samples_19.txt}
export THREADS=${THREADS:-32}
