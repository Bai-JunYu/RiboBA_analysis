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
  echo "Usage: bash run_ribORF.sh [PUBLIC_INPUTS_ROOT] [OUTPUT_ROOT]" >&2
  exit 1
fi

RIBORF_ROOT="${OUTPUT_ROOT}/riborf"
PREP_DIR="${RIBORF_ROOT}/prepare"
RESULT_DIR="${RIBORF_ROOT}/results"

mkdir -p "${PREP_DIR}/bowtie2_index" "${PREP_DIR}/transcriptome_data" "${RESULT_DIR}"

bowtie2-build -f "${PUBLIC_INPUTS_DIR}/data/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa" \
  "${PREP_DIR}/bowtie2_index/hsa"

tophat2 -G "${PUBLIC_INPUTS_DIR}/data/gtf/Homo_sapiens.GRCh38.109.gtf" \
  --transcriptome-index="${PREP_DIR}/transcriptome_data/known" \
  "${PREP_DIR}/bowtie2_index/hsa"

for fq in "${PUBLIC_INPUTS_DIR}"/data/p4/lnc_circRNA/*fq; do
  sample_file="$(basename "${fq}")"
  sample="${sample_file%_*}"
  mkdir -p "${RESULT_DIR}/${sample}"
  tophat2 --no-novel-juncs --library-type fr-unstranded --read-realign-edit-dist 0 --segment-length 20 --min-anchor-length 5 --min-intron-length 50 --num-threads 5 \
    --GTF "${PUBLIC_INPUTS_DIR}/data/gtf/Homo_sapiens.GRCh38.109.gtf" \
    --transcriptome-index="${PREP_DIR}/transcriptome_data/known" \
    --output-dir "${RESULT_DIR}/${sample}" \
    --no-convert-bam -p 15 --no-coverage-search \
    "${PREP_DIR}/bowtie2_index/hsa" \
    "${fq}"
done

"${PUBLIC_INPUTS_DIR}/software/gtfToGenePred" \
  "${PUBLIC_INPUTS_DIR}/data/gtf/Homo_sapiens.GRCh38.109.gtf" \
  "${PREP_DIR}/Homo_sapiens.GRCh38.109.genePred.txt"

perl "${PUBLIC_INPUTS_DIR}/software/RibORF.2.0/ORFannotate.pl" \
  -g "${PUBLIC_INPUTS_DIR}/data/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa" \
  -t "${PREP_DIR}/Homo_sapiens.GRCh38.109.genePred.txt" \
  -o "${PREP_DIR}" \
  -s ATG/CTG/GTG/TTG/AAG/ACG/AGG/ATA/ATC/ATT -l 15

CANDIDATE_ORF="${PREP_DIR}/candidateORF.genepred.txt"
[[ -f "${CANDIDATE_ORF}" ]] || { echo "Error: missing ${CANDIDATE_ORF}" >&2; exit 1; }

cd "${RESULT_DIR}" || exit 1
for sample_dir in ./*/; do
  [[ -d "${sample_dir}" ]] || continue
  cd "${sample_dir}" || exit 1
  : > offset.corretion.parameters.txt
  perl "${PUBLIC_INPUTS_DIR}/software/RibORF.2.0/readDist.pl" -f accepted_hits.sam -g \
    "${CANDIDATE_ORF}" -o ./plot
  mkdir -p ribORF_result
  perl "${PUBLIC_INPUTS_DIR}/software/RibORF.2.0/offsetCorrect.pl" \
    -r accepted_hits.sam \
    -p offset.corretion.parameters.txt \
    -o corrected.mapping.sam

  perl "${PUBLIC_INPUTS_DIR}/software/RibORF.2.0/ribORF.pl" \
    -f corrected.mapping.sam \
    -c "${CANDIDATE_ORF}" \
    -o ribORF_result
  cd .. || exit 1
done

echo "RibORF done. Outputs: ${RIBORF_ROOT}" >&2
