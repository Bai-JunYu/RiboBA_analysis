# Tool Versions and Requirements

Pinned/validated dependencies for `tools_orf_prediction/*` and `simulation/*/run_tools_*_benchmark.sh`.

## Conda-managed tools

### Python 2.7 tool environment (`environment/conda_env_tools_py27.yml`)
- python: 2.7.15
- tophat2: 2.1.1
- bowtie: 1.2.3
- bowtie2: 2.2.5

### Python 3.7 tool environment (`environment/conda_env_tools_py37.yml`)
- python: 3.7.12
- ribocode: 1.2.15
- ribotish: 0.2.7
- samtools: 1.6
- sambamba: 1.0.0
- bowtie: 1.3.1
- bowtie2: 2.5.1
- STAR: 2.7.10b
- cutadapt: 4.4
- pigz: conda package

## Additional required commands (system / external install)
- `rsync` (used in simulation benchmark tool runners)
- `conda` (scripts activate `py2.7` / `py3.7`)

## Tools expected under `${PUBLIC_INPUTS_DIR}/software`
- Gedi 1.0.6 binary:
  - `${PUBLIC_INPUTS_DIR}/software/Gedi_1.0.6/gedi`
- ORF-RATER scripts:
  - `${PUBLIC_INPUTS_DIR}/software/ORF-RATER/*`
- RibORF 2.0 scripts:
  - `${PUBLIC_INPUTS_DIR}/software/RibORF.2.0/*`
- UCSC converters:
  - `${PUBLIC_INPUTS_DIR}/software/gtfToGenePred`
  - `${PUBLIC_INPUTS_DIR}/software/genePredToBed`

## RiboCode binaries required on PATH
- `prepare_transcripts`
- `metaplots`
- `RiboCode`

## Environment files
- py2.7 tools: `environment/conda_env_tools_py27.yml`
- py3.7 tools: `environment/conda_env_tools_py37.yml`
