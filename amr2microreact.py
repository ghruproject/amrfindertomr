#!/usr/bin/env python3
from __future__ import annotations

"""
amr2microreact.py - Convert AMR tool output(s) to Microreact-compatible CSV.

Supports: AMRFinderPlus, ABRicate, ResFinder, CARD RGI.
Auto-detects input format from column headers.

Reads one or more output files and produces a CSV with:
  - id column (sample name, matching tree tip labels)
  - Drug class summary columns (comma-separated gene lists per class)
  - Drug class __colour columns for Microreact colouring
  - Individual gene presence/absence columns (yes/no)

Usage:
    python amr2microreact.py -i amr_output_dir/ -o microreact_metadata.csv
    python amr2microreact.py -i file1.tsv file2.tsv -o microreact_metadata.csv
"""

import argparse
import base64
import csv
import json
import os
import sys
from collections import defaultdict
from pathlib import Path
from urllib.request import Request, urlopen
from urllib.error import HTTPError


# ---------------------------------------------------------------------------
# Colour scheme for Microreact
# ---------------------------------------------------------------------------

# Colours for drug class summary columns
COLOUR_HAS_GENES = "#E53935"   # red - resistance genes present
COLOUR_NO_GENES = "#43A047"    # green - no resistance genes
# Colours for gene presence columns
COLOUR_PRESENT = "#E53935"     # red
COLOUR_ABSENT = "#EEEEEE"     # light grey


# ---------------------------------------------------------------------------
# Format detection and parsing
# ---------------------------------------------------------------------------

# Signature columns used to identify each format
FORMAT_SIGNATURES = {
    "amrfinderplus": {"Element symbol", "Element name", "Scope"},
    "abricate": {"#FILE", "GENE", "%COVERAGE", "%IDENTITY"},
    "resfinder": {"Resistance gene", "Phenotype", "Accession no."},
    "rgi": {"Best_Hit_ARO", "Drug Class", "Best_Hit_Identity"},
}


def detect_format(filepath: Path) -> str:
    """Auto-detect the AMR tool format from column headers."""
    with open(filepath) as f:
        header_line = f.readline().strip()

    # Handle tab-separated and comma-separated
    if "\t" in header_line:
        columns = set(header_line.split("\t"))
    else:
        columns = set(header_line.split(","))

    for fmt, signature in FORMAT_SIGNATURES.items():
        if signature.issubset(columns):
            return fmt

    return "unknown"


def parse_tsv(filepath: Path) -> tuple[list[dict], str]:
    """Parse a TSV/CSV file and return (rows, delimiter)."""
    with open(filepath) as f:
        first_line = f.readline()
    delimiter = "\t" if "\t" in first_line else ","
    with open(filepath) as f:
        reader = csv.DictReader(f, delimiter=delimiter)
        return list(reader), delimiter


def normalize_row(row: dict, fmt: str, filepath: Path) -> dict | None:
    """Normalize a row from any supported format into a common structure.

    Returns a dict with keys: sample, gene, drug_class, element_type,
    coverage, identity, scope. Returns None if the row should be skipped.
    """
    if fmt == "amrfinderplus":
        gene = row.get("Element symbol", "")
        if not gene or gene == "NA":
            return None
        return {
            "sample": row.get("Name", "") or filepath.stem.replace(".result", "").replace("_amr", ""),
            "gene": gene,
            "drug_class": row.get("Class", "NA"),
            "element_type": row.get("Type", ""),
            "coverage": _parse_float(row.get("% Coverage of reference", "0")),
            "identity": _parse_float(row.get("% Identity to reference", "0")),
            "scope": row.get("Scope", ""),
        }

    elif fmt == "abricate":
        gene = row.get("GENE", "")
        if not gene:
            return None
        sample = row.get("#FILE", "")
        if sample:
            sample = Path(sample).stem.replace(".result", "")
        else:
            sample = filepath.stem.replace("_abricate", "")
        return {
            "sample": sample,
            "gene": gene,
            "drug_class": row.get("RESISTANCE", row.get("PRODUCT", "NA")),
            "element_type": "AMR",
            "coverage": _parse_float(row.get("%COVERAGE", "0")),
            "identity": _parse_float(row.get("%IDENTITY", "0")),
            "scope": "",
        }

    elif fmt == "resfinder":
        gene = row.get("Resistance gene", "")
        if not gene:
            return None
        return {
            "sample": filepath.stem.replace("_resfinder", "").replace("ResFinder_results_tab", ""),
            "gene": gene,
            "drug_class": row.get("Phenotype", "NA"),
            "element_type": "AMR",
            "coverage": _parse_float(row.get("Coverage", "0").rstrip("%")),
            "identity": _parse_float(row.get("Identity", "0").rstrip("%")),
            "scope": "",
        }

    elif fmt == "rgi":
        gene = row.get("Best_Hit_ARO", "")
        if not gene:
            return None
        orf = row.get("ORF_ID", "")
        sample = orf.split(" ")[0] if orf else filepath.stem.replace("_rgi", "")
        # Clean up contig prefix from sample name
        if "_NODE_" in sample or "_contig" in sample:
            parts = sample.split("_NODE_")
            if len(parts) > 1:
                sample = parts[0]
            else:
                parts = sample.split("_contig")
                if len(parts) > 1:
                    sample = parts[0]
        return {
            "sample": sample,
            "gene": gene,
            "drug_class": row.get("Drug Class", "NA"),
            "element_type": row.get("Model_type", "AMR"),
            "coverage": 100.0,  # RGI doesn't report coverage directly
            "identity": _parse_float(row.get("Best_Hit_Identity", "0")),
            "scope": row.get("Cut_Off", ""),
        }

    return None


def _parse_float(val: str) -> float:
    """Safely parse a float, returning 0.0 on failure."""
    try:
        return float(val)
    except (ValueError, TypeError):
        return 0.0


def extract_sample_name(filepath: Path, rows: list[dict]) -> str:
    """Get sample name from the Name column or fall back to filename."""
    if rows:
        name = rows[0].get("Name", "")
        if name:
            return name
    return filepath.stem.replace(".result", "").replace("_amr", "")


def collect_inputs(input_paths: list[str]) -> list[Path]:
    """Resolve input paths - can be files or directories."""
    files = []
    for p in input_paths:
        path = Path(p)
        if path.is_dir():
            files.extend(sorted(path.glob("*.tsv")))
            files.extend(sorted(path.glob("*.txt")))
            files.extend(sorted(path.glob("*.csv")))
        elif path.is_file():
            files.append(path)
        else:
            print(f"Warning: {p} not found, skipping", file=sys.stderr)
    return files


def build_metadata(
    input_files: list[Path],
    amr_only: bool = False,
    min_coverage: float = 0.0,
    min_identity: float = 0.0,
    scope_filter: str | None = None,
    class_filter: set[str] | None = None,
):
    """Build the metadata table from AMR tool outputs.

    Args:
        input_files: List of input file paths.
        amr_only: Only include AMR-type genes (exclude STRESS, VIRULENCE).
        min_coverage: Minimum coverage % threshold.
        min_identity: Minimum identity % threshold.
        scope_filter: Only include genes with this scope (e.g. "core", "plus").
        class_filter: Only include these drug classes (set of strings).

    Returns (samples, drug_classes, all_genes, formats_detected) where:
      - samples: dict of sample data
      - drug_classes: sorted list of drug class names
      - all_genes: sorted list of all gene symbols
      - formats_detected: set of format names encountered
    """
    samples = {}
    all_classes = set()
    all_genes = set()
    formats_detected = set()

    for filepath in input_files:
        fmt = detect_format(filepath)
        if fmt == "unknown":
            print(f"Warning: Unknown format for {filepath}, skipping", file=sys.stderr)
            continue

        formats_detected.add(fmt)
        rows, _ = parse_tsv(filepath)
        if not rows:
            continue

        for row in rows:
            normalized = normalize_row(row, fmt, filepath)
            if normalized is None:
                continue

            sample_name = normalized["sample"]
            gene = normalized["gene"]
            drug_class = normalized["drug_class"]
            element_type = normalized["element_type"]
            coverage = normalized["coverage"]
            identity = normalized["identity"]
            scope = normalized["scope"]

            # Apply filters
            if amr_only and element_type not in ("AMR",):
                continue
            if coverage < min_coverage:
                continue
            if identity < min_identity:
                continue
            if scope_filter and scope and scope.lower() != scope_filter.lower():
                continue
            if class_filter and drug_class not in class_filter:
                continue

            if sample_name not in samples:
                samples[sample_name] = {
                    "classes": defaultdict(set),
                    "genes": set(),
                }

            samples[sample_name]["genes"].add(gene)
            all_genes.add(gene)

            if drug_class and drug_class != "NA":
                samples[sample_name]["classes"][drug_class].add(gene)
                all_classes.add(drug_class)

    drug_classes = sorted(all_classes)
    gene_list = sorted(all_genes)

    return samples, drug_classes, gene_list, formats_detected


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
    add_colours: bool = True,
):
    """Write the Microreact-compatible CSV."""
    header = ["id"]
    # Drug class summary columns + optional colour columns
    for dc in drug_classes:
        col = sanitize_column_name(dc)
        header.append(col)
        if add_colours:
            header.append(f"{col}__colour")
    # Individual gene columns + optional colour columns
    for gene in gene_list:
        col = sanitize_column_name(gene)
        header.append(col)
        if add_colours:
            header.append(f"{col}__colour")

    rows_out = []
    for sample_name in sorted(samples.keys()):
        sdata = samples[sample_name]
        row = [sample_name]

        for dc in drug_classes:
            genes_in_class = sdata["classes"].get(dc, set())
            summary = ",".join(sorted(genes_in_class)) if genes_in_class else "NA"
            row.append(summary)
            if add_colours:
                row.append(COLOUR_HAS_GENES if genes_in_class else COLOUR_NO_GENES)

        for gene in gene_list:
            present = gene in sdata["genes"]
            row.append("yes" if present else "no")
            if add_colours:
                row.append(COLOUR_PRESENT if present else COLOUR_ABSENT)

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
    if add_colours:
        print("  Microreact __colour columns included", file=sys.stderr)


def build_csv_string(
    samples: dict,
    drug_classes: list[str],
    gene_list: list[str],
    add_colours: bool = True,
) -> str:
    """Build CSV content as a string (for API upload)."""
    import io

    header = ["id"]
    for dc in drug_classes:
        col = sanitize_column_name(dc)
        header.append(col)
        if add_colours:
            header.append(f"{col}__colour")
    for gene in gene_list:
        col = sanitize_column_name(gene)
        header.append(col)
        if add_colours:
            header.append(f"{col}__colour")

    buf = io.StringIO()
    writer = csv.writer(buf)
    writer.writerow(header)

    for sample_name in sorted(samples.keys()):
        sdata = samples[sample_name]
        row = [sample_name]
        for dc in drug_classes:
            genes_in_class = sdata["classes"].get(dc, set())
            summary = ",".join(sorted(genes_in_class)) if genes_in_class else "NA"
            row.append(summary)
            if add_colours:
                row.append(COLOUR_HAS_GENES if genes_in_class else COLOUR_NO_GENES)
        for gene in gene_list:
            present = gene in sdata["genes"]
            row.append("yes" if present else "no")
            if add_colours:
                row.append(COLOUR_PRESENT if present else COLOUR_ABSENT)
        writer.writerow(row)

    return buf.getvalue()


# ---------------------------------------------------------------------------
# Microreact API integration
# ---------------------------------------------------------------------------

MICROREACT_API_URL = "https://microreact.org/api/projects/create"


def build_microreact_project_json(
    csv_data: str,
    tree_data: str | None = None,
    project_name: str = "AMR Microreact Project",
) -> dict:
    """Build a Microreact .microreact JSON document with embedded data."""
    csv_b64 = base64.b64encode(csv_data.encode()).decode()

    project = {
        "schema": "https://microreact.org/schema/v1.json",
        "meta": {"name": project_name},
        "datasets": {
            "dataset-1": {
                "file": "data-file-1",
                "idFieldName": "id",
            }
        },
        "files": {
            "data-file-1": {
                "id": "data-file-1",
                "format": "text/csv",
                "name": "metadata.csv",
                "url": f"data:text/csv;base64,{csv_b64}",
            }
        },
        "tables": {
            "table-1": {
                "dataset": "dataset-1",
                "title": "AMR Metadata",
            }
        },
    }

    if tree_data:
        tree_b64 = base64.b64encode(tree_data.encode()).decode()
        project["files"]["tree-file-1"] = {
            "id": "tree-file-1",
            "format": "text/x-nh",
            "name": "tree.nwk",
            "url": f"data:text/x-nh;base64,{tree_b64}",
        }
        project["trees"] = {
            "tree-1": {
                "file": "tree-file-1",
                "title": "Phylogenetic Tree",
            }
        }

    return project


def upload_to_microreact(
    csv_data: str,
    api_key: str,
    tree_data: str | None = None,
    project_name: str = "AMR Microreact Project",
) -> dict:
    """Upload data to Microreact and return the response with project URL.

    Returns dict with 'url' key on success, or raises an exception.
    """
    project_json = build_microreact_project_json(csv_data, tree_data, project_name)
    body = json.dumps(project_json).encode("utf-8")

    req = Request(
        MICROREACT_API_URL,
        data=body,
        headers={
            "Content-Type": "application/json; charset=utf-8",
            "Access-Token": api_key,
        },
        method="POST",
    )

    try:
        with urlopen(req) as resp:
            return json.loads(resp.read().decode())
    except HTTPError as e:
        error_body = e.read().decode() if e.fp else str(e)
        raise RuntimeError(f"Microreact API error ({e.code}): {error_body}") from e


def main():
    parser = argparse.ArgumentParser(
        description="Convert AMR tool output to Microreact metadata CSV. "
        "Supports AMRFinderPlus, ABRicate, ResFinder, and CARD RGI. "
        "Format is auto-detected."
    )
    parser.add_argument(
        "-i",
        "--input",
        nargs="+",
        required=True,
        help="Input file(s) or directory containing them",
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
    parser.add_argument(
        "--min-coverage",
        type=float,
        default=0.0,
        help="Minimum coverage %% threshold (default: 0)",
    )
    parser.add_argument(
        "--min-identity",
        type=float,
        default=0.0,
        help="Minimum identity %% threshold (default: 0)",
    )
    parser.add_argument(
        "--scope",
        choices=["core", "plus"],
        default=None,
        help="Only include genes with this scope (AMRFinderPlus only)",
    )
    parser.add_argument(
        "--classes",
        nargs="+",
        default=None,
        help="Only include these drug classes",
    )
    parser.add_argument(
        "--no-colours",
        action="store_true",
        help="Do not add __colour columns for Microreact",
    )
    parser.add_argument(
        "--tree",
        default=None,
        help="Newick tree file to include when uploading to Microreact",
    )
    parser.add_argument(
        "--microreact-api-key",
        default=None,
        help="Microreact API access token (or set MICROREACT_API_KEY env var)",
    )
    parser.add_argument(
        "--project-name",
        default="AMR Microreact Project",
        help="Name for the Microreact project (default: 'AMR Microreact Project')",
    )
    args = parser.parse_args()

    input_files = collect_inputs(args.input)
    if not input_files:
        print("Error: No input files found", file=sys.stderr)
        sys.exit(1)

    print(f"Processing {len(input_files)} files...", file=sys.stderr)

    class_filter = set(args.classes) if args.classes else None
    samples, drug_classes, gene_list, formats = build_metadata(
        input_files,
        amr_only=args.amr_only,
        min_coverage=args.min_coverage,
        min_identity=args.min_identity,
        scope_filter=args.scope,
        class_filter=class_filter,
    )

    print(f"  Formats detected: {', '.join(sorted(formats))}", file=sys.stderr)

    add_colours = not args.no_colours
    write_csv(samples, drug_classes, gene_list, Path(args.output), add_colours=add_colours)

    # Microreact API upload
    api_key = args.microreact_api_key or os.environ.get("MICROREACT_API_KEY")
    if api_key:
        csv_data = build_csv_string(samples, drug_classes, gene_list, add_colours=add_colours)

        tree_data = None
        if args.tree:
            tree_path = Path(args.tree)
            if tree_path.is_file():
                tree_data = tree_path.read_text()
            else:
                print(f"Warning: Tree file {args.tree} not found, uploading without tree", file=sys.stderr)

        print("Uploading to Microreact...", file=sys.stderr)
        try:
            result = upload_to_microreact(csv_data, api_key, tree_data, args.project_name)
            url = result.get("url", "")
            print(f"Microreact project created: {url}", file=sys.stderr)
        except RuntimeError as e:
            print(f"Error: {e}", file=sys.stderr)
            sys.exit(1)


if __name__ == "__main__":
    main()
