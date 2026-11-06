# Genomic position of managed Korean *Apis cerana*

Analysis code for:

> Lee J.-H., Kim H., Frunze O., Nah G. and Kwon H. W. Genomic position of managed Korean *Apis cerana* relative to Chinese and Japanese lineages. *Journal of Apiculture* (in review).

PacBio HiFi reads from 19 managed Korean *A. cerana* colonies (one worker each) are aligned to the Korean reference genome AcerK_1.0. The code then:

- calls structural variants (SVs) and SNPs;
- computes an ancestry index at diagnostic sites defined from four Chinese assemblies;
- places the colonies relative to 36 public Chinese and Japanese individuals;
- tests whether declared strain names and apiaries locate colonies on the resulting axis;
- assesses gene content.

## Data

| Data | Source |
|---|---|
| HiFi reads of the 19 colonies (FASTQ) | NABIC (https://nabic.rda.go.kr) (to be deposited upon publication) |
| Reference genome AcerK_1.0 | NCBI GCF_029169275.1 |
| RefSeq annotation | NCBI GCF_029169275.2, release RS_2026_08 |
| Chinese assemblies for diagnostic sites | GCA_011100585.1 (Hainan), GCA_020141525.2 (Sichuan), GCA_029531695.1 (Hubei), GCA_039655025.1 (Qinghai). These are downloaded by `04_ancestry_index/01_*.sh` |
| 36 reference individuals (short reads) | SRA/DRA accessions in `metadata/reference_individuals.tsv` (Supplementary Table S1) |

## Layout and run order

Run the folders in order, and the scripts in each folder in numerical order. Each folder corresponds to a section of the Materials and Methods.

| Folder | Methods section | Content |
|---|---|---|
| `00_setup/` | – | BAM lists and metadata for the reference individuals |
| `01_alignment/` | 4, 6 | HiFi reads to AcerK_1.0 (minimap2); reference individuals to AcerK_1.0 (bwa mem -M, GATK MarkDuplicates); depth and coverage (Table 1) |
| `02_structural_variants/` | 4, 5 | Sniffles2 per colony; multi-sample merge; genotype matrix and filters; Jaccard PCoA |
| `03_snps/` | 5 | Joint SNP calling, filtering, MAF and distance thinning; PCA and distances |
| `04_ancestry_index/` | 6 | Diagnostic sites from the Chinese assemblies; genotyping of colonies and reference individuals; index with 1 Mb block bootstrap; depth-matched check |
| `05_joint_pca/` | 7 | Joint genotyping of colonies and reference individuals; pseudo-haploid PCA |
| `06_structure_and_labels/` | 5, 9 | See the list below |
| `07_gene_content/` | 8 | RefSeq GFF conversion; CDS depth (mosdepth); gene presence and threshold sensitivity |
| `figures/` | – | Figs 1–4 and Supplementary Figs S1–S3 |
| `metadata/` | 1 | Colony list and declared metadata (Table 1); reference individuals |
| `results/` | – | Small derived tables read by the figure scripts (see below) |

`06_structure_and_labels/` contains:

- clustering, silhouette and gap tests, with silhouette calibration;
- distance-based RDA;
- declared-strain tests and *Hallabeol* cohesion;
- apiary variance components;
- label tests and minimum detectable effect;
- SV callability sensitivity. Run it as `python 06_structure_and_labels/06_sv_callability_sensitivity.py 1.0 <existing_outdir>`.

### Before running

1. Install the tools listed in `environment.yml`.
2. Edit the paths in `config.sh`.
3. Prepare the reference genome: `samtools faidx $REF` and `bwa index $REF`.
4. Download these inputs:
   - the RefSeq GFF `GCF_029169275.2_AcerK_1.0_genomic.gff.gz` into `$KH_REFSEQ_DIR`;
   - the HiFi FASTQ files from NABIC into `$RAW_DIR`;
   - the reference-individual FASTQ files, e.g. with `fasterq-dump`, into `$REFBAM_DIR/fastq`.
5. Run `source config.sh` in bash before running Python scripts. The Python scripts read `KH_ROOT`, `KH_RESULTS` and `KH_REFSEQ_DIR` from the environment. Shell scripts source `config.sh` themselves.

### Figures

The figures need only `results/`. From the repository root, run `python figures/make_fig3.py` and so on; output goes to `figures/out/`. The figure scripts reproduce the published figures pixel for pixel; see `environment.yml` for the plotting versions used for that check.

### Origin of `results/`

| File | Produced by (file name in `$KH_RESULTS`) |
|---|---|
| `two_axes_colonies.tsv`, `two_axes_anchors.tsv` | `04_ancestry_index/08_index_bootstrap_and_summary.py` |
| `joint_pca/pca_pseudohap_oriented.tsv`, `joint_pca/summary.json` | `05_joint_pca/02_joint_pca.py` (`joint_pca_snp/`) |
| `snp_vs_sv_structure.tsv`, `silhouette_fst_calibration.tsv` | `06_structure_and_labels/01_structure_summary.py` (`summary_tables/`) |
| `strain_effect_tests.tsv` | `06_structure_and_labels/03_declared_strain_tests.py` (`label_tests/`) |
| `rs_gene_presence_read.tsv` | `07_gene_content/03_gene_presence.py` (`read_gene_presence/gene_presence_read.tsv`) |
| `sites_19.tsv`, `colony_labels_index.tsv` | Columns of the tables above combined with `metadata/colonies.tsv` (index = `china_index_poly`) |
| `ref_meta.txt` | `metadata/reference_individuals.tsv` (group and locality) |

## What was changed from the analysis as run

The scripts are those that produced the published results. Changes were limited to the following:

- **Paths.** Server paths were replaced by the variables in `config.sh`.
- **Inputs.** Input lists were restricted to the 19 colonies (`metadata/samples_19.txt`).
- **Removed code.** Code for analyses not reported in the paper was removed.
- **Split drivers.** Two driver scripts were split into `02_structural_variants/02_joint_sv_merge.sh` and `03_snps/01_call_snps.sh`.
- **Gene content.** `07_gene_content/` reads the RefSeq annotation directly, and `03_gene_presence.py` also computes the threshold sensitivity reported in Results section 5.
- **Small steps.** Steps that had been run by hand were added as code:
  - subsampling of diagnostic sites in `04/04`;
  - allele-depth export in `04/04`;
  - PCA sign orientation in `05/02`.

These scripts are new in this release:

- `00_setup/make_reference_lists.sh`;
- `01_alignment/01_align_hifi.sh`, `02_align_reference_individuals.sh` and `03_coverage.sh`. They use the same commands and settings as the original alignment and coverage runs.
- `06_structure_and_labels/05_label_tests.py`. It reproduces the label statistics in Results section 4, including the minimum detectable difference of 0.169.

Comments in the original scripts are in Korean.

## Glossary of internal names

| Name | Meaning |
|---|---|
| `china_index_poly`, `KC` | The ancestry index |
| `two_axes_*` | Index tables for colonies and reference individuals |
| `china_ref/` | Diagnostic-site files |
| `Halla` | *Hallabeol* |

## Licence

MIT (see `LICENSE`).
