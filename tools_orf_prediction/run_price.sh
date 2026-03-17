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
  echo "Usage: bash run_price.sh [PUBLIC_INPUTS_ROOT] [OUTPUT_ROOT]" >&2
  exit 1
fi

PRICE_ROOT="${OUTPUT_ROOT}/price"
PREP_DIR="${PRICE_ROOT}/prepare"
MAP_TX_DIR="${PRICE_ROOT}/intermediate/map_tx"
RESULT_DIR="${PRICE_ROOT}/results"

mkdir -p "${PREP_DIR}/index/hsa" "${PREP_DIR}/indexgenome" "${MAP_TX_DIR}" "${RESULT_DIR}"

STAR --runThreadN 15 --runMode genomeGenerate --genomeSAsparseD 2 \
  --genomeDir "${PREP_DIR}/index/hsa" \
  --genomeFastaFiles "${PUBLIC_INPUTS_DIR}/data/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa" \
  --sjdbGTFfile "${PUBLIC_INPUTS_DIR}/data/gtf/Homo_sapiens.GRCh38.109.gtf" \
  --sjdbOverhang 30

for fq in "${PUBLIC_INPUTS_DIR}"/data/p4/lnc_circRNA/*fq; do
  sample_file="$(basename "${fq}")"
  sample="${sample_file%_*}"
  STAR --runThreadN 15 --genomeDir "${PREP_DIR}/index/hsa" \
    --readFilesIn "${fq}" --outSAMtype BAM SortedByCoordinate --outSAMattributes MD NH --alignEndsType EndToEnd \
    --outFileNamePrefix "${MAP_TX_DIR}/${sample}_" \
    --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
    --outFilterMismatchNmax 2 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --chimSegmentMin 20
done

cd "${MAP_TX_DIR}" || exit 1

for bam in ./*.bam; do
  sambamba sort -o "${bam%_*}_sorted.bam" "${bam}"
done

OML="${PREP_DIR}/indexgenome/Homo_sapiens.GRCh38.109.oml"
"${PUBLIC_INPUTS_DIR}/software/Gedi_1.0.6/gedi" -e IndexGenome \
  -s "${PUBLIC_INPUTS_DIR}/data/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa" \
  -a "${PUBLIC_INPUTS_DIR}/data/gtf/Homo_sapiens.GRCh38.109.gtf" \
  -f "${PREP_DIR}/indexgenome" \
  -o "${OML}" \
  -nobowtie -nokallisto -nostar

for bam in ./*_sorted.bam; do
  name="$(basename "${bam%_*}")"
  "${PUBLIC_INPUTS_DIR}/software/Gedi_1.0.6/gedi" -e Price -reads "${bam}" \
    -genomic "${OML}" \
    -prefix "${RESULT_DIR}/${name}"
done

for cit in "${RESULT_DIR}"/*.orfs.cit; do
  "${PUBLIC_INPUTS_DIR}/software/Gedi_1.0.6/gedi" \
    -e ViewCIT -m bed -name 'd.getType()' "${cit}" > "${cit%.*}.orfs.bed"
done

echo "PRICE done. Outputs: ${PRICE_ROOT}" >&2
