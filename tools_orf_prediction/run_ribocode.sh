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
  echo "Usage: bash run_ribocode.sh [PUBLIC_INPUTS_ROOT] [OUTPUT_ROOT]" >&2
  exit 1
fi

RIBOCODE_ROOT="${OUTPUT_ROOT}/ribocode"
PREP_DIR="${RIBOCODE_ROOT}/prepare"
MAP_TX_DIR="${RIBOCODE_ROOT}/intermediate/map_tx"
RESULT_DIR="${RIBOCODE_ROOT}/results"
STAR_INDEX_DIR="${OUTPUT_ROOT}/price/prepare/index/hsa"

if [[ ! -d "${STAR_INDEX_DIR}" ]]; then
  echo "Error: STAR index not found: ${STAR_INDEX_DIR}" >&2
  echo "Run run_price.sh first (or provide matching OUTPUT_ROOT)." >&2
  exit 1
fi

mkdir -p "${PREP_DIR}" "${MAP_TX_DIR}" "${RESULT_DIR}"

prepare_transcripts -g "${PUBLIC_INPUTS_DIR}/data/gtf/Homo_sapiens.GRCh38.109.gtf" \
  -f "${PUBLIC_INPUTS_DIR}/data/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa" \
  -o "${PREP_DIR}"

for fq in "${PUBLIC_INPUTS_DIR}"/data/p4/lnc_circRNA/*norrna.fq; do
  sample_file="$(basename "${fq}")"
  sample="${sample_file%_*}"
  STAR --runThreadN 15 \
    --genomeDir "${STAR_INDEX_DIR}" \
    --readFilesIn "${fq}" \
    --outFilterMismatchNmax 2 \
    --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts \
    --outFilterMultimapNmax 1 --alignEndsType EndToEnd \
    --outFileNamePrefix "${MAP_TX_DIR}/${sample}_"
done

cd "${MAP_TX_DIR}" || exit 1

for bam in ./*Aligned.toTranscriptome.out.bam; do
  sambamba sort -o "${bam%Ali*}_sorted.bam" "${bam}"
done

for bam in ./*sorted.bam; do
  sample="$(basename "${bam%_*}")"
  metaplots -a "${PREP_DIR}" \
    -r "${bam}" -s no -f0_percent 0.4 -o "${sample}" -m 19 -M 42

  cfg="${sample}_pre_config.txt"
  if [[ ! -f "${cfg}" ]]; then
    echo "Warning: skip ${sample}, missing ${cfg}" >&2
    continue
  fi

  RiboCode -a "${PREP_DIR}" \
    -c "${cfg}" -l no -b -A CTG,GTG,TTG,ACG,AAG,AGG,ATA,ATT,ATC -o "${RESULT_DIR}/${sample}"
done

echo "RiboCode done. Outputs: ${RIBOCODE_ROOT}" >&2
