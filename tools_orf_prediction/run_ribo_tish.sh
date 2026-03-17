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
  echo "Usage: bash run_ribo_tish.sh [PUBLIC_INPUTS_ROOT] [OUTPUT_ROOT]" >&2
  exit 1
fi

RIBOTISH_ROOT="${OUTPUT_ROOT}/ribotish"
INPUT_DIR="${RIBOTISH_ROOT}/results"
mkdir -p "${INPUT_DIR}"

cd "${INPUT_DIR}" || exit 1

for sample_dir in ./*/; do
  [[ -d "${sample_dir}" ]] || continue
  sample="$(basename "${sample_dir%/}")"
  cd "${sample_dir}" || exit 1
  bam="${sample}_sorted.bam"
  if [[ ! -f "${bam}" ]]; then
    echo "Warning: skip ${sample}, missing ${bam}" >&2
    cd .. || exit 1
    continue
  fi

  ribotish quality -b "${bam}" \
    -g "${PUBLIC_INPUTS_DIR}/data/gtf/Homo_sapiens.GRCh38.109.gtf" \
    -p 15 -v -l 19,42 --th 0.34

  ribotish predict --alt -b "${bam}" \
    -g "${PUBLIC_INPUTS_DIR}/data/gtf/Homo_sapiens.GRCh38.109.gtf" \
    -f "${PUBLIC_INPUTS_DIR}/data/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa" \
    -p 15 -v --framebest \
    -o pred.txt
  cd .. || exit 1
done

echo "RiboTISH done. Outputs: ${RIBOTISH_ROOT}" >&2
