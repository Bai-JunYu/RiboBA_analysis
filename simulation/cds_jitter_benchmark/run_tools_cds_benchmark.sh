#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
DEFAULT_PUBLIC_INPUTS_DIR="${PROJECT_DIR}/public_inputs"
PUBLIC_INPUTS_DIR="${1:-${PUBLIC_INPUTS_DIR:-${DEFAULT_PUBLIC_INPUTS_DIR}}}"
SIM_DIR="${PUBLIC_INPUTS_DIR}/data/simulate1/simu_fastq"

require_dir() {
  local path="$1"
  [[ -d "${path}" ]] || { echo "Error: required directory not found: ${path}" >&2; exit 1; }
}

activate_conda_env() {
  local env_name="$1"
  if ! command -v conda >/dev/null 2>&1; then
    echo "Error: conda command not found (required for env ${env_name})." >&2
    exit 1
  fi
  # shellcheck disable=SC1091
  source "$(conda info --base)/etc/profile.d/conda.sh"
  conda activate "${env_name}"
}

if [[ ! -d "${PUBLIC_INPUTS_DIR}" ]]; then
  echo "Error: PUBLIC_INPUTS_DIR not found: ${PUBLIC_INPUTS_DIR}" >&2
  echo "Usage: bash run_tools_cds_benchmark.sh [PUBLIC_INPUTS_ROOT]" >&2
  exit 1
fi

require_dir "${SIM_DIR}"
require_dir "${PUBLIC_INPUTS_DIR}/result/prepare/for_price/index/hsa"
require_dir "${PUBLIC_INPUTS_DIR}/result/prepare/for_price/indexgenome"
require_dir "${PUBLIC_INPUTS_DIR}/result/prepare/ribocode"
require_dir "${PUBLIC_INPUTS_DIR}/result/prepare/riborf/transcriptome_data"
require_dir "${PUBLIC_INPUTS_DIR}/result/prepare/riborf/bowtie2_index"

cd "${SIM_DIR}" || exit 1
mkdir -p orfrater price riboba ribocode riborf ribotish

# price
for fastq in ./*.fastq; do
  sample_file="$(basename "${fastq}")"
  sample="${sample_file%_SRR*}"
  STAR --runThreadN 15 --genomeDir "${PUBLIC_INPUTS_DIR}/result/prepare/for_price/index/hsa" \
    --readFilesIn "${fastq}" --outSAMtype BAM SortedByCoordinate --outSAMattributes MD NH --alignEndsType EndToEnd \
    --outFileNamePrefix "price/${sample}" \
    --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 \
    --outFilterMismatchNoverLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --chimSegmentMin 20
done

cd price || exit 1
for bam in ./*Aligned.sortedByCoord.out.bam; do
  sambamba sort -o "${bam%Aligned*}_sorted.bam" "${bam}"
done

for bam in ./*sorted.bam; do
  "${PUBLIC_INPUTS_DIR}/software/Gedi_1.0.6/gedi" -e Price -reads "${bam}" \
    -genomic "${PUBLIC_INPUTS_DIR}/result/prepare/for_price/indexgenome/Homo_sapiens.GRCh38.109.oml" \
    -prefix "${SIM_DIR}/price/${bam%_sorted*}"
done

# orfrater
activate_conda_env py2.7
mkdir -p ../orfrater
for bam in ./*sorted.bam; do
  bam_base="$(basename "${bam}")"
  sample="${bam_base%_sorted.bam}"
  sample_dir="../orfrater/${sample}"
  mkdir -p "${sample_dir}"
  mv "${bam}" "${sample_dir}/${sample}_sorted.bam"
  for bai in "${bam}.bai" "${bam%.bam}.bai"; do
    if [[ -f "${bai}" ]]; then
      mv "${bai}" "${sample_dir}/${sample}_sorted.bam.bai"
      break
    fi
  done
done

cd ../orfrater || exit 1
for sample_dir in ./*/; do
  sample="$(basename "${sample_dir%/}")"
  cd "${sample_dir}" || exit 1
  bam="${sample}_sorted.bam"

  "${PUBLIC_INPUTS_DIR}/software/ORF-RATER/prune_transcripts.py" \
    --inbed "${PUBLIC_INPUTS_DIR}/result/prepare/orf_rater/Homo_sapiens.GRCh38.109.bed" \
    -p 15 \
    -v "${PUBLIC_INPUTS_DIR}/data/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa" \
    "${bam}" > prune.log

  "${PUBLIC_INPUTS_DIR}/software/ORF-RATER/make_tfams.py" -v > make_tfams.log

  "${PUBLIC_INPUTS_DIR}/software/ORF-RATER/find_orfs_and_types.py" \
    "${PUBLIC_INPUTS_DIR}/data/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa" \
    --codons ATG -p 15 -v > find_ORF.log

  "${PUBLIC_INPUTS_DIR}/software/ORF-RATER/psite_trimmed.py" "${bam}" --minrdlen 28 --maxrdlen 41 \
    --tallyfile tallies.txt \
    --cdsbed "${PUBLIC_INPUTS_DIR}/result/prepare/orf_rater/Homo_sapiens.GRCh38.109.bed" \
    -p 15 -v > psite.log

  "${PUBLIC_INPUTS_DIR}/software/ORF-RATER/regress_orfs.py" "${bam}" -p 1 -v > regress.log

  "${PUBLIC_INPUTS_DIR}/software/ORF-RATER/rate_regression_output.py" \
    ./ -p 15 --CSV rate_regression.csv -v > rate.log

  "${PUBLIC_INPUTS_DIR}/software/ORF-RATER/make_orf_bed.py" --minrating 0
  cd .. || exit 1
done

# ribotish
activate_conda_env py3.7
mkdir -p ../ribotish
rsync -aR --include='*/' --include='*sorted.bam' --include='*sorted.bam.bai' --exclude='*' ./ ../ribotish
cd ../ribotish || exit 1
for sample_dir in ./*/; do
  sample="$(basename "${sample_dir%/}")"
  cd "${sample_dir}" || exit 1
  bam="${sample}_sorted.bam"
  ribotish quality -b "${bam}" \
    -g "${PUBLIC_INPUTS_DIR}/data/gtf/Homo_sapiens.GRCh38.109.gtf" \
    -p 15 -v -l 19,42 --th 0.34

  ribotish predict -b "${bam}" \
    -g "${PUBLIC_INPUTS_DIR}/data/gtf/Homo_sapiens.GRCh38.109.gtf" \
    -f "${PUBLIC_INPUTS_DIR}/data/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa" \
    -p 10 -v --framebest --fpth 1 --fspth 1 --fsqth 1 \
    -o pred.txt
  cd .. || exit 1
done

# ribocode
cd "${SIM_DIR}" || exit 1
for fastq in ./*.fastq; do
  sample_file="$(basename "${fastq}")"
  sample="${sample_file%_SRR*}"
  STAR --runThreadN 10 \
    --genomeDir "${PUBLIC_INPUTS_DIR}/result/prepare/for_price/index/hsa" \
    --readFilesIn "${fastq}" \
    --outFilterMismatchNmax 2 \
    --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts \
    --outFilterMultimapNmax 1 --alignEndsType EndToEnd \
    --outFileNamePrefix "ribocode/${sample}"
done

cd ribocode || exit 1
for bam in ./*Aligned.toTranscriptome.out.bam; do
  sambamba sort -o "${bam%Aligned*}_sorted.bam" "${bam}"
done

for bam in ./*sorted.bam; do
  sample="${bam%_sorted*}"
  metaplots -a "${PUBLIC_INPUTS_DIR}/result/prepare/ribocode" \
    -r "${bam}" -s no -f0_percent 0.34 -o "${sample}" -m 19 -M 42 -pv1 0.01 -pv2 0.01

  RiboCode -a "${PUBLIC_INPUTS_DIR}/result/prepare/ribocode" \
    -c "${sample}_pre_config.txt" -l no -b -o "${sample}" -p 1
done

# ribORF
perl "${PUBLIC_INPUTS_DIR}/software/RibORF.2.0/ORFannotate.pl" \
  -g "${PUBLIC_INPUTS_DIR}/data/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa" \
  -t "${PUBLIC_INPUTS_DIR}/data/gtf/Homo_sapiens.GRCh38.109.genePred.txt" \
  -o "${PUBLIC_INPUTS_DIR}/result/prepare/riborf_atg" \
  -s ATG -l 10

activate_conda_env py2.7
cd "${SIM_DIR}" || exit 1
for fastq in ./*.fastq; do
  sample_file="$(basename "${fastq}")"
  sample="${sample_file%_SRR*}"
  mkdir -p "${SIM_DIR}/riborf/${sample}"
  tophat2 --no-novel-juncs --library-type fr-unstranded --read-realign-edit-dist 0 --segment-length 20 --min-anchor-length 5 --min-intron-length 50 --num-threads 5 \
    --GTF "${PUBLIC_INPUTS_DIR}/data/gtf/Homo_sapiens.GRCh38.109.gtf" \
    --transcriptome-index="${PUBLIC_INPUTS_DIR}/result/prepare/riborf/transcriptome_data/known" \
    --output-dir "${SIM_DIR}/riborf/${sample}" \
    --no-convert-bam -p 10 --no-coverage-search \
    "${PUBLIC_INPUTS_DIR}/result/prepare/riborf/bowtie2_index/hsa" \
    "${fastq}"
done

cd riborf || exit 1
for sample_dir in ./*/; do
  cd "${sample_dir}" || exit 1
  mkdir -p ribORF_result
  perl "${PUBLIC_INPUTS_DIR}/software/RibORF.2.0/offsetCorrect.pl" \
    -r accepted_hits.sam \
    -p offset.corretion.parameters.txt \
    -o corrected.mapping.sam

  perl "${PUBLIC_INPUTS_DIR}/software/RibORF.2.0/ribORF.pl" -p 0 \
    -f corrected.mapping.sam \
    -c "${PUBLIC_INPUTS_DIR}/result/prepare/riborf_atg/candidateORF.genepred.txt" \
    -o ribORF_result
  cd .. || exit 1
done

# riboba
cd "${SIM_DIR}" || exit 1
: > cds.txt
for fastq in ./*.fastq; do
  sample_file="$(basename "${fastq}")"
  sample="${sample_file%_SRR*}"
  echo "${sample_file}" >> cds.txt
  bowtie -p 12 -v 1 -a -m 10 --best --strata -S --no-unal \
    -x "${PUBLIC_INPUTS_DIR}/result/prepare/index/tx_cds/cds" \
    "${fastq}" \
    > "riboba/${sample}.sam"
done

cd riboba || exit 1
for sam in ./*.sam; do
  sambamba view -S -f bam "${sam}" | \
    sambamba sort -o "${sam%.*}_txsorted.bam" /dev/stdin
done
