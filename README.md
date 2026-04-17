<p align="center">
  <img src="logo.svg" alt="amr2microreact logo" width="120"/>
</p>

# amr2microreact

[![CI](https://github.com/ghruproject/amrfindertomr/actions/workflows/ci.yml/badge.svg)](https://github.com/ghruproject/amrfindertomr/actions/workflows/ci.yml)

Convert AMR tool output into [Microreact](https://microreact.org)-compatible metadata for interactive visualization with phylogenetic trees.

**Supports:** [AMRFinderPlus](https://github.com/ncbi/amr) | [ABRicate](https://github.com/tseemann/abricate) | [ResFinder](https://cge.food.dtu.dk/services/ResFinder/) | [CARD RGI](https://github.com/arpcard/rgi) — format auto-detected.

## Web App

**No install required** — use the hosted web app:

**[amrfindertomr-web.streamlit.app](https://amrfindertomr-web.streamlit.app)**

Upload your AMR output files, adjust filters, preview the metadata table, download the CSV, or create a Microreact project directly from your browser.

## Features

- **Multi-format support** — auto-detects AMRFinderPlus, ABRicate, ResFinder, and CARD RGI from column headers. Mix formats in a single run.
- **Microreact colour columns** — auto-generates `__colour` columns (red = resistance genes present, green = absent) for instant visual mapping.
- **Filtering** — filter by scope (core/plus), drug class, minimum coverage %, and minimum identity %.
- **Microreact API integration** — create Microreact projects directly from the CLI or web app. Upload metadata + tree and get back a shareable URL.
- **Scales to 400+ samples** — handles large datasets efficiently.

## Output format

| id | BETA.LACTAM | BETA.LACTAM__colour | blaTEM.1 | blaTEM.1__colour | ... |
|----|-------------|---------------------|----------|------------------|-----|
| sample1 | blaTEM-1 | #E53935 | yes | #E53935 | ... |
| sample2 | NA | #43A047 | no | #EEEEEE | ... |

- Drug class columns: comma-separated gene names or `NA`
- Gene columns: `yes` / `no`
- `__colour` columns: hex colours for Microreact rendering
- `id` column must match tree tip labels

## Command-line usage

### Setup

```bash
# Install pixi (https://pixi.sh), then:
pixi install
```

### Quick start

```bash
# Convert AMRFinderPlus outputs (auto-detected)
pixi run python amr2microreact.py -i amr_results/ -o metadata.csv

# Mix formats - auto-detected
pixi run python amr2microreact.py -i amrfinder.tsv abricate.tsv rgi.txt -o metadata.csv

# Filter: AMR only, min 80% coverage, min 90% identity
pixi run python amr2microreact.py -i results/ -o metadata.csv \
  --amr-only --min-coverage 80 --min-identity 90

# Filter by scope and drug class
pixi run python amr2microreact.py -i results/ -o metadata.csv \
  --scope core --classes BETA-LACTAM AMINOGLYCOSIDE

# Disable colour columns
pixi run python amr2microreact.py -i results/ -o metadata.csv --no-colours
```

### Create a Microreact project via API

You can create a Microreact project directly from the command line. Get your API access token from [your Microreact account settings](https://microreact.org/my-account).

```bash
# Pass API key directly
pixi run python amr2microreact.py -i results/ -o metadata.csv \
  --tree tree.nwk \
  --microreact-api-key YOUR_TOKEN \
  --project-name "My AMR Study"

# Or use an environment variable
export MICROREACT_API_KEY=your_token
pixi run python amr2microreact.py -i results/ -o metadata.csv --tree tree.nwk
```

The CSV is always saved locally. If an API key is provided, the project is also created on Microreact and the URL is printed to stderr.

### Generate a tree (optional)

```bash
docker run -v $(pwd)/assemblies:/data staphb/mashtree \
  bash -c "mashtree /data/*.fasta" > tree.nwk
```

## Examples

The `examples/` directory contains 3 AMRFinderPlus output files for testing:

```bash
pixi run python amr2microreact.py -i examples/ -o test_output.csv
```

## Development

```bash
pixi install
pixi run pytest tests/ -v
```

## License

MIT
