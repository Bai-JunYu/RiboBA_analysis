# RiboBA-analysis

Reproducible code repository for RiboBA manuscript analyses, including simulation benchmarks, external tool comparisons, downstream validation, and figure generation.

## What is included
- Analysis scripts under `simulation/`, `tools_orf_prediction/`, and `downstream_validation/`.
- Main and supplementary plotting scripts under `figures/`.
- Configuration snapshots under `config/`.
- Environment and version captures under `environment/`.
- Reproduction notes under `docs/quickstart.md`.

## Data location
- This repository stores code only.
- Primary analysis data are expected under `../../RiboBA_analysis_data` (relative to repo root).
- Large raw/reference inputs (FASTQ/FASTA/GTF etc.) should be downloaded by users from public sources.
- Optionally set `RIBOBA_DATA_DIR` to override the external data root.

## Quick start
See `docs/quickstart.md`.

## Conda environments
- Tool environment (Python 2.7): `environment/conda_env_tools_py27.yml` (`py2.7`, used by ORF-RATER and RibORF)
- Tool environment (Python 3.7): `environment/conda_env_tools_py37.yml` (`py3.7`, used by RiboTISH, RiboCode, and PRICE dependencies)
- Tool/software version registry: `environment/tool_versions.md`
- Dependency checker: `environment/check_tool_dependencies.sh [PUBLIC_INPUTS_ROOT]`

## Notes
- Full raw-data reruns require external datasets/tools plus user-downloaded public inputs.
- Figure wrappers auto-resolve figure-ready `.rds` from `RIBOBA_DATA_DIR` or default to `../../RiboBA_analysis_data`.
