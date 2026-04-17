# amr2microreact

Convert [AMRFinderPlus](https://github.com/ncbi/amr) output into [Microreact](https://microreact.org)-compatible metadata for interactive visualization with phylogenetic trees.

## Web App

**No install required** — use the hosted web app:

**[amrfindertomr.streamlit.app](https://amrfindertomr.streamlit.app)**

Upload your AMRFinderPlus TSV files, preview the metadata table, and download the Microreact-ready CSV directly from your browser.

## What it does

Takes one or more AMRFinderPlus TSV output files and produces a single CSV with:

- **`id`** column — sample identifiers matching your tree tip labels
- **Drug class summary columns** — comma-separated gene lists per antimicrobial class (e.g. `BETA-LACTAM`, `AMINOGLYCOSIDE`)
- **Gene presence/absence columns** — `yes`/`no` for each detected gene across all samples

Upload this CSV alongside a Newick tree to [Microreact](https://microreact.org) to get interactive metadata blocks with your phylogeny.

## Command-line usage

For batch processing or integration into pipelines, use the CLI directly.

### Setup

Requires [pixi](https://pixi.sh):

```bash
pixi install
```

This installs Python, pandas, and AMRFinderPlus.

For tree generation, [mashtree](https://github.com/lskatz/mashtree) is used via Docker:

```bash
docker pull staphb/mashtree
```

### 1. Run AMRFinderPlus on your assemblies

```bash
# Single sample
pixi run amrfinder --nucleotide sample.fasta --name sample_id -o sample_amr.tsv

# Batch (parallel)
for f in assemblies/*.fasta; do
  name=$(basename "$f" .fasta)
  pixi run amrfinder --nucleotide "$f" --name "$name" -o "amr_results/${name}_amr.tsv" &
done
wait
```

### 2. Convert to Microreact format

```bash
# From a directory of AMRFinderPlus outputs
pixi run python amr2microreact.py -i amr_results/ -o microreact_metadata.csv

# From specific files
pixi run python amr2microreact.py -i sample1_amr.tsv sample2_amr.tsv -o metadata.csv

# AMR genes only (exclude stress/virulence)
pixi run python amr2microreact.py -i amr_results/ -o metadata.csv --amr-only
```

### 3. Generate a tree (optional)

```bash
docker run -v $(pwd)/assemblies:/data staphb/mashtree \
  mashtree /data/*.fasta > tree.nwk
```

### 4. Upload to Microreact

Go to [microreact.org](https://microreact.org), upload your `microreact_metadata.csv` and `tree.nwk`, and explore your AMR data interactively.

## Output format

| id | BETA.LACTAM | AMINOGLYCOSIDE | blaTEM.1 | aadA1 | ... |
|----|-------------|----------------|----------|-------|-----|
| sample1 | blaTEM-1 | aadA1 | yes | yes | ... |
| sample2 | NA | aadA1,aph(6)-Id | no | yes | ... |

- Drug class columns contain comma-separated gene names or `NA`
- Gene columns contain `yes` or `no`
- The `id` column must match tree tip labels for Microreact to link them

## License

MIT
