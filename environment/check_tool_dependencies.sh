#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
DEFAULT_PUBLIC_INPUTS_DIR="${PROJECT_DIR}/public_inputs"
PUBLIC_INPUTS_DIR="${1:-${PUBLIC_INPUTS_DIR:-${DEFAULT_PUBLIC_INPUTS_DIR}}}"

need_cmds=(
  conda rsync bowtie bowtie2-build STAR samtools sambamba perl tophat2
  prepare_transcripts metaplots RiboCode ribotish
)

need_files=(
  "${PUBLIC_INPUTS_DIR}/software/Gedi_1.0.6/gedi"
  "${PUBLIC_INPUTS_DIR}/software/gtfToGenePred"
  "${PUBLIC_INPUTS_DIR}/software/genePredToBed"
  "${PUBLIC_INPUTS_DIR}/software/RibORF.2.0/ORFannotate.pl"
  "${PUBLIC_INPUTS_DIR}/software/ORF-RATER/prune_transcripts.py"
)

fail=0

echo "[Check] commands on PATH"
for c in "${need_cmds[@]}"; do
  if command -v "${c}" >/dev/null 2>&1; then
    echo "  OK   ${c} -> $(command -v "${c}")"
  else
    echo "  MISS ${c}"
    fail=1
  fi
done

echo "[Check] files under PUBLIC_INPUTS_DIR=${PUBLIC_INPUTS_DIR}"
for f in "${need_files[@]}"; do
  if [[ -e "${f}" ]]; then
    echo "  OK   ${f}"
  else
    echo "  MISS ${f}"
    fail=1
  fi
done

if [[ ${fail} -ne 0 ]]; then
  echo "Dependency check failed." >&2
  exit 1
fi

echo "Dependency check passed."
