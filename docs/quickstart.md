# Quickstart

This repository is code-only.

## 1) Directory layout

From repo root (`RiboBA-analysis`):

- code: `./`
- external data root (default): `../RiboBA_analysis_data/`

Expected figure-ready inputs:

- `../RiboBA_analysis_data/figure_ready_data/*.rds`

## 2) Setup environments

```bash
cd /path/to/RiboBA-analysis

conda env create -f environment/conda_env_tools_py27.yml
conda env create -f environment/conda_env_tools_py37.yml
```

## 3) Check external dependencies

`RiboBase` is required by ORF-tool workflows.

```bash
bash environment/check_tool_dependencies.sh /path/to/RiboBase
```

## 4) Run figures (example)

```bash
cd /path/to/RiboBA-analysis
Rscript figures/fig2.R
Rscript figures/fig3.R
```

Outputs are written to `figures/*.pdf`.

## 5) Run ORF prediction tools (example)

```bash
cd /path/to/RiboBA-analysis

# [arg1] RiboBase root, [arg2] output root
bash tools_orf_prediction/run_price.sh /path/to/RiboBase ../RiboBA_analysis_data/processed_data/tools_orf_prediction
```

Each tool writes to:

- `../RiboBA_analysis_data/processed_data/tools_orf_prediction/<tool>/`

## 6) Reproducibility notes

- All scripts are expected to run with relative paths from repo root.
- Optional override for data root: `RIBOBA_DATA_DIR=/abs/path/to/RiboBA_analysis_data`.
- Tool parameter snapshots are under `config/tool_params/*.yml`.
