#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
DEFAULT_PUBLIC_INPUTS_DIR="${PROJECT_DIR}/public_inputs"
DEFAULT_OUTPUT_ROOT="${PROJECT_DIR}/../RiboBA_analysis_data/processed_data/tools_orf_prediction"
PUBLIC_INPUTS_DIR="${1:-${PUBLIC_INPUTS_DIR:-${DEFAULT_PUBLIC_INPUTS_DIR}}}"
OUTPUT_ROOT="${2:-${OUTPUT_ROOT:-${DEFAULT_OUTPUT_ROOT}}}"

if [[ ! -d "${PUBLIC_INPUTS_DIR}" ]]; then
  echo "Error: PUBLIC_INPUTS_DIR not found: ${PUBLIC_INPUTS_DIR}" >&2
  echo "Usage: bash run_orf_rater.sh [PUBLIC_INPUTS_ROOT] [OUTPUT_ROOT]" >&2
  exit 1
fi

ORFRATER_ROOT="${OUTPUT_ROOT}/orf_rater"
PREP_DIR="${ORFRATER_ROOT}/prepare"
RESULT_DIR="${ORFRATER_ROOT}/results"
SOURCE_DIR="${OUTPUT_ROOT}/price/intermediate/map_tx"

mkdir -p "${PREP_DIR}" "${RESULT_DIR}"

if [[ ! -d "${SOURCE_DIR}" ]]; then
  echo "Error: source BAM dir not found: ${SOURCE_DIR}" >&2
  echo "Run run_price.sh first (same OUTPUT_ROOT)." >&2
  exit 1
fi

BED_FILE="${PREP_DIR}/Homo_sapiens.GRCh38.109.bed"
"${PUBLIC_INPUTS_DIR}/software/gtfToGenePred" -ignoreGroupsWithoutExons \
  "${PUBLIC_INPUTS_DIR}/data/gtf/Homo_sapiens.GRCh38.109.gtf" stdout | \
  "${PUBLIC_INPUTS_DIR}/software/genePredToBed" stdin "${BED_FILE}"

for bam in "${SOURCE_DIR}"/*sorted.bam; do
  bam_base="$(basename "${bam}")"
  sample="${bam_base%_sorted.bam}"
  sample_dir="${RESULT_DIR}/${sample}"
  mkdir -p "${sample_dir}"
  cp -f "${bam}" "${sample_dir}/${sample}_sorted.bam"

  for bai in "${bam}.bai" "${bam%.bam}.bai"; do
    if [[ -f "${bai}" ]]; then
      cp -f "${bai}" "${sample_dir}/${sample}_sorted.bam.bai"
      break
    fi
  done
done

declare -A MIN_RDLEN=(
  [SRR8449566]=28 [SRR23242344]=28 [SRR23242345]=28 [SRR23242346]=30 [SRR23242347]=31
  [SRR2433794]=23 [SRR7073124]=36 [SRR8449567]=24 [SRR8449568]=19
)
declare -A MAX_RDLEN=(
  [SRR8449566]=35 [SRR23242344]=32 [SRR23242345]=33 [SRR23242346]=41 [SRR23242347]=39
  [SRR2433794]=33 [SRR7073124]=38 [SRR8449567]=34 [SRR8449568]=30
)

for sample_dir in "${RESULT_DIR}"/*/; do
  [[ -d "${sample_dir}" ]] || continue
  sample="$(basename "${sample_dir%/}")"
  cd "${sample_dir}" || exit 1
  bam="${sample}_sorted.bam"
  if [[ ! -f "${bam}" ]]; then
    echo "Warning: skip ${sample}, missing ${bam}" >&2
    cd "${RESULT_DIR}" || exit 1
    continue
  fi

  "${PUBLIC_INPUTS_DIR}/software/ORF-RATER/prune_transcripts.py" \
    --inbed "${BED_FILE}" -p 15 \
    -v "${PUBLIC_INPUTS_DIR}/data/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa" \
    "${bam}" > prune.log

  "${PUBLIC_INPUTS_DIR}/software/ORF-RATER/make_tfams.py" -v > make_tfams.log

  "${PUBLIC_INPUTS_DIR}/software/ORF-RATER/find_orfs_and_types.py" \
    "${PUBLIC_INPUTS_DIR}/data/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa" \
    --codons ATG CTG GTG TTG ACG AAG AGG ATA ATT ATC -p 15 -v > find_ORF.log

  if [[ -n "${MIN_RDLEN[$sample]:-}" && -n "${MAX_RDLEN[$sample]:-}" ]]; then
    "${PUBLIC_INPUTS_DIR}/software/ORF-RATER/psite_trimmed.py" "${bam}" \
      --minrdlen "${MIN_RDLEN[$sample]}" --maxrdlen "${MAX_RDLEN[$sample]}" \
      --tallyfile tallies.txt --cdsbed "${BED_FILE}" -p 15 -v > psite.log
  else
    echo "Warning: no min/max read length configured for ${sample}, skip psite_trimmed.py" >&2
  fi

  "${PUBLIC_INPUTS_DIR}/software/ORF-RATER/regress_orfs.py" "${bam}" -p 1 -v > regress.log
  "${PUBLIC_INPUTS_DIR}/software/ORF-RATER/rate_regression_output.py" ./ -p 15 --CSV rate_regression.csv -v > rate.log
  "${PUBLIC_INPUTS_DIR}/software/ORF-RATER/make_orf_bed.py" --minrating 0
  "${PUBLIC_INPUTS_DIR}/software/ORF-RATER/quantify_orfs.py" "${bam}" -p 15 -v --CSV quant.csv --minrating 0

  cd "${RESULT_DIR}" || exit 1
done

echo "ORF-RATER done. Outputs: ${ORFRATER_ROOT}" >&2
