#!/usr/bin/env python3
"""
amr2microreact.py - Convert AMRFinderPlus output(s) to Microreact-compatible CSV.

Reads one or more AMRFinderPlus TSV output files and produces a CSV with:
  - id column (sample name, matching tree tip labels)
  - Drug class summary columns (comma-separated gene lists per class)
  - Individual gene presence/absence columns (yes/no)

Usage:
    python amr2microreact.py -i amr_output_dir/ -o microreact_metadata.csv
    python amr2microreact.py -i file1.tsv file2.tsv -o microreact_metadata.csv
"""

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path


def parse_amrfinder(filepath: Path) -> list[dict]:
    """Parse a single AMRFinderPlus TSV output file."""
    rows = []
    with open(filepath) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            rows.append(row)
    return rows


def extract_sample_name(filepath: Path, rows: list[dict]) -> str:
    """Get sample name from the Name column or fall back to filename."""
    if rows:
        name = rows[0].get("Name", "")
        if name:
            return name
    # Fall back to filename without extension
    return filepath.stem.replace(".result", "").replace("_amr", "")


def collect_inputs(input_paths: list[str]) -> list[Path]:
    """Resolve input paths - can be files or directories."""
    files = []
    for p in input_paths:
        path = Path(p)
        if path.is_dir():
            files.extend(sorted(path.glob("*.tsv")))
            files.extend(sorted(path.glob("*.txt")))
        elif path.is_file():
            files.append(path)
        else:
            print(f"Warning: {p} not found, skipping", file=sys.stderr)
    return files


def build_metadata(input_files: list[Path], amr_only: bool = False):
    """Build the metadata table from AMRFinderPlus outputs.

    Returns (rows, drug_classes, all_genes) where:
      - rows: list of dicts with sample data
      - drug_classes: sorted list of drug class names
      - all_genes: sorted list of all gene symbols
    """
    samples = {}  # sample_name -> {class -> set(genes), genes -> set}
    all_classes = set()
    all_genes = set()

    for filepath in input_files:
        rows = parse_amrfinder(filepath)
        if not rows:
            continue

        sample_name = extract_sample_name(filepath, rows)
        if sample_name not in samples:
            samples[sample_name] = {
                "classes": defaultdict(set),
                "genes": set(),
                "filepath": filepath,
            }

        for row in rows:
            scope = row.get("Scope", "")
            element_type = row.get("Type", "")
            gene = row.get("Element symbol", "")
            drug_class = row.get("Class", "")
            subclass = row.get("Subclass", "")

            if not gene or gene == "NA":
                continue

            # Optionally filter to AMR only (skip STRESS, VIRULENCE)
            if amr_only and element_type not in ("AMR",):
                continue

            samples[sample_name]["genes"].add(gene)
            all_genes.add(gene)

            # Use subclass if available and different from class
            if drug_class and drug_class != "NA":
                samples[sample_name]["classes"][drug_class].add(gene)
                all_classes.add(drug_class)

    drug_classes = sorted(all_classes)
    gene_list = sorted(all_genes)

    return samples, drug_classes, gene_list


def sanitize_column_name(name: str) -> str:
    """Make column names Microreact-friendly (replace problematic chars)."""
    return name.replace("(", ".").replace(")", ".").replace("'", ".").replace(
        '"', "."
    ).replace("-", ".").replace("/", ".")


def write_csv(
    samples: dict,
    drug_classes: list[str],
    gene_list: list[str],
    output_path: Path,
):
    """Write the Microreact-compatible CSV."""
    # Build header
    header = ["id"]
    # Drug class summary columns
    for dc in drug_classes:
        header.append(sanitize_column_name(dc))
    # Individual gene columns
    for gene in gene_list:
        header.append(sanitize_column_name(gene))

    rows_out = []
    for sample_name in sorted(samples.keys()):
        sdata = samples[sample_name]
        row = [sample_name]

        # Drug class summaries (comma-separated gene lists)
        for dc in drug_classes:
            genes_in_class = sdata["classes"].get(dc, set())
            row.append(",".join(sorted(genes_in_class)) if genes_in_class else "NA")

        # Individual gene presence
        for gene in gene_list:
            row.append("yes" if gene in sdata["genes"] else "no")

        rows_out.append(row)

    with open(output_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(rows_out)

    print(f"Wrote {len(rows_out)} samples to {output_path}", file=sys.stderr)
    print(
        f"  {len(drug_classes)} drug class columns, {len(gene_list)} gene columns",
        file=sys.stderr,
    )


def main():
    parser = argparse.ArgumentParser(
        description="Convert AMRFinderPlus output to Microreact metadata CSV"
    )
    parser.add_argument(
        "-i",
        "--input",
        nargs="+",
        required=True,
        help="Input AMRFinderPlus TSV file(s) or directory containing them",
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output CSV file path",
    )
    parser.add_argument(
        "--amr-only",
        action="store_true",
        help="Only include AMR genes (exclude STRESS, VIRULENCE)",
    )
    args = parser.parse_args()

    input_files = collect_inputs(args.input)
    if not input_files:
        print("Error: No input files found", file=sys.stderr)
        sys.exit(1)

    print(f"Processing {len(input_files)} AMRFinderPlus output files...", file=sys.stderr)

    samples, drug_classes, gene_list = build_metadata(input_files, args.amr_only)
    write_csv(samples, drug_classes, gene_list, Path(args.output))


if __name__ == "__main__":
    main()
