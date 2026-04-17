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


# ---------------------------------------------------------------------------
# Colour scheme for Microreact
# ---------------------------------------------------------------------------

# Default colours for drug class summary columns
COLOUR_HAS_GENES = "#E53935"   # red - resistance genes present
COLOUR_NO_GENES = "#43A047"    # green - no resistance genes
# Default colours for gene presence columns
COLOUR_PRESENT = "#E53935"     # red
COLOUR_ABSENT = "#EEEEEE"     # light grey

# Named colour palettes
COLOUR_PALETTES = {
    "default": {"present": "#E53935", "absent_class": "#43A047", "absent_gene": "#EEEEEE"},
    "blue-grey": {"present": "#1565C0", "absent_class": "#78909C", "absent_gene": "#ECEFF1"},
    "purple-teal": {"present": "#7B1FA2", "absent_class": "#00897B", "absent_gene": "#F3E5F5"},
    "orange-blue": {"present": "#E65100", "absent_class": "#0277BD", "absent_gene": "#FFF3E0"},
    "pink-green": {"present": "#C2185B", "absent_class": "#2E7D32", "absent_gene": "#FCE4EC"},
    "amber-indigo": {"present": "#FF8F00", "absent_class": "#283593", "absent_gene": "#FFF8E1"},
}


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


def _resolve_colours(palette: str | None = None) -> tuple[str, str, str]:
    """Return (present, absent_class, absent_gene) colour hex values."""
    if palette and palette in COLOUR_PALETTES:
        p = COLOUR_PALETTES[palette]
        return p["present"], p["absent_class"], p["absent_gene"]
    return COLOUR_HAS_GENES, COLOUR_NO_GENES, COLOUR_ABSENT


def write_csv(
    samples: dict,
    drug_classes: list[str],
    gene_list: list[str],
    output_path: Path,
    add_colours: bool = True,
    colour_palette: str | None = None,
):
    """Write the Microreact-compatible CSV."""
    col_present, col_absent_class, col_absent_gene = _resolve_colours(colour_palette)

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

    rows_out = []
    for sample_name in sorted(samples.keys()):
        sdata = samples[sample_name]
        row = [sample_name]

        for dc in drug_classes:
            genes_in_class = sdata["classes"].get(dc, set())
            summary = ",".join(sorted(genes_in_class)) if genes_in_class else "NA"
            row.append(summary)
            if add_colours:
                row.append(col_present if genes_in_class else col_absent_class)

        for gene in gene_list:
            present = gene in sdata["genes"]
            row.append("yes" if present else "no")
            if add_colours:
                row.append(col_present if present else col_absent_gene)

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
    colour_palette: str | None = None,
) -> str:
    """Build CSV content as a string (for API upload)."""
    import io

    col_present, col_absent_class, col_absent_gene = _resolve_colours(colour_palette)

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
                row.append(col_present if genes_in_class else col_absent_class)
        for gene in gene_list:
            present = gene in sdata["genes"]
            row.append("yes" if present else "no")
            if add_colours:
                row.append(col_present if present else col_absent_gene)
        writer.writerow(row)

    return buf.getvalue()


# ---------------------------------------------------------------------------
# Microreact API integration
# ---------------------------------------------------------------------------

MICROREACT_API_URL = "https://microreact.org/api/projects/create"


def _build_pane_layout(has_tree: bool) -> dict:
    """Build the Microreact pane layout model."""
    if has_tree:
        # Outer row has 1 child (a nested row), nested row stacks vertically
        layout_children = [
            {
                "type": "row",
                "children": [
                    {
                        "type": "tabset",
                        "weight": 60,
                        "children": [
                            {"type": "tab", "id": "tree-1", "name": "Tree", "component": "Tree"}
                        ],
                    },
                    {
                        "type": "tabset",
                        "weight": 40,
                        "children": [
                            {"type": "tab", "id": "table-1", "name": "Metadata", "component": "Table"}
                        ],
                    },
                ],
            },
        ]
    else:
        layout_children = [
            {
                "type": "tabset",
                "children": [
                    {"type": "tab", "id": "table-1", "name": "Metadata", "component": "Table"}
                ],
            },
        ]

    return {
        "model": {
            "global": {
                "splitterSize": 2,
                "tabEnableClose": False,
                "tabSetHeaderHeight": 1,
                "tabSetTabStripHeight": 1,
                "tabSetMinWidth": 160,
                "tabSetMinHeight": 160,
                "borderMinSize": 160,
                "borderBarSize": 20,
                "borderEnableDrop": False,
            },
            "borders": [
                {
                    "type": "border",
                    "size": 240,
                    "location": "right",
                    "children": [
                        {"type": "tab", "id": "--mr-legend-pane", "name": "Legend", "component": "Legend", "enableClose": False, "enableDrag": False},
                        {"type": "tab", "id": "--mr-selection-pane", "name": "Selection", "component": "Selection", "enableClose": False, "enableDrag": False},
                        {"type": "tab", "id": "--mr-history-pane", "name": "History", "component": "History", "enableClose": False, "enableDrag": False},
                        {"type": "tab", "id": "--mr-views-pane", "name": "Views", "component": "Views", "enableClose": False, "enableDrag": False},
                    ],
                }
            ],
            "layout": {
                "type": "row",
                "children": layout_children,
            },
        }
    }


def _build_tree_config(file_id: str) -> dict:
    """Build full Microreact tree panel configuration."""
    return {
        "alignLabels": True,
        "blockHeaderFontSize": 13,
        "blockPadding": 0,
        "blocks": [],
        "blockSize": 14,
        "branchLengthsDigits": 4,
        "controls": False,
        "fontSize": 16,
        "hideOrphanDataRows": False,
        "ids": None,
        "internalLabelsFilterRange": [0, 100],
        "internalLabelsFontSize": 13,
        "lasso": False,
        "nodeSize": 14,
        "path": None,
        "roundBranchLengths": True,
        "scaleLineAlpha": True,
        "showBlockHeaders": True,
        "showBlockLabels": False,
        "showBranchLengths": False,
        "showEdges": True,
        "showInternalLabels": False,
        "showLabels": True,
        "showLeafLabels": False,
        "showPiecharts": True,
        "showShapeBorders": True,
        "showShapes": True,
        "styleLeafLabels": False,
        "styleNodeEdges": False,
        "subtreeIds": None,
        "type": "rc",
        "title": "Tree",
        "labelField": "id",
        "file": file_id,
    }


def build_microreact_project_json(
    csv_data: str,
    tree_data: str | None = None,
    project_name: str = "AMR Microreact Project",
) -> dict:
    """Build a Microreact .microreact JSON document with embedded data."""
    csv_bytes = csv_data.encode()
    csv_b64 = base64.b64encode(csv_bytes).decode()

    # Parse CSV header to build table column definitions
    first_line = csv_data.split("\n")[0]
    columns = [{"field": col.strip(), "fixed": False} for col in first_line.split(",")]

    project = {
        "schema": "https://microreact.org/schema/v1.json",
        "charts": {},
        "datasets": {
            "dataset-1": {
                "id": "dataset-1",
                "file": "data-file-1",
                "idFieldName": "id",
            }
        },
        "files": {
            "data-file-1": {
                "id": "data-file-1",
                "format": "text/csv",
                "name": "metadata.csv",
                "type": "data",
                "size": len(csv_bytes),
                "blob": f"data:application/octet-stream;base64,{csv_b64}",
            }
        },
        "filters": {
            "paneFilters": [],
            "dataFilters": [],
            "chartFilters": [],
            "searchOperator": "includes",
            "searchValue": "",
            "selection": [],
            "selectionBreakdownField": None,
        },
        "maps": {},
        "matrices": {},
        "meta": {"name": project_name},
        "networks": {},
        "notes": {},
        "slicers": {},
        "styles": {
            "coloursField": None,
            "colourPalettes": [],
            "defaultColour": "transparent",
            "defaultShape": "circle",
            "colourSettings": {},
            "labelsField": None,
            "legendDirection": "row",
            "shapesField": None,
            "shapePalettes": [],
        },
        "tables": {
            "table-1": {
                "displayMode": "cosy",
                "hideUnselected": False,
                "dataset": "dataset-1",
                "file": "data-file-1",
                "title": "AMR Metadata",
                "paneId": "table-1",
                "columns": columns,
            }
        },
        "timelines": {},
        "views": [],
    }

    has_tree = tree_data is not None
    if has_tree:
        tree_bytes = tree_data.encode()
        tree_b64 = base64.b64encode(tree_bytes).decode()
        project["files"]["tree-file-1"] = {
            "id": "tree-file-1",
            "format": "text/x-nh",
            "name": "tree.nwk",
            "type": "tree",
            "size": len(tree_bytes),
            "blob": f"data:application/octet-stream;base64,{tree_b64}",
        }
        project["trees"] = {
            "tree-1": _build_tree_config("tree-file-1"),
        }

    project["panes"] = _build_pane_layout(has_tree)

    return project


def upload_to_microreact(
    csv_data: str,
    api_key: str,
    tree_data: str | None = None,
    project_name: str = "AMR Microreact Project",
) -> dict:
    """Upload data to Microreact and return the response with project URL.

    Returns dict with 'url' key on success, or raises an exception.
    Uses requests if available (better Cloudflare compat), falls back to urllib.
    """
    project_json = build_microreact_project_json(csv_data, tree_data, project_name)

    headers = {
        "Content-Type": "application/json; charset=utf-8",
        "Access-Token": api_key,
        "User-Agent": "Mozilla/5.0 (compatible; amr2microreact/1.0)",
        "Accept": "application/json",
    }

    try:
        import requests
        resp = requests.post(MICROREACT_API_URL, json=project_json, headers=headers)
        if resp.status_code != 200:
            raise RuntimeError(f"Microreact API error ({resp.status_code}): {resp.text}")
        return resp.json()
    except ImportError:
        from urllib.request import Request, urlopen
        from urllib.error import HTTPError

        body = json.dumps(project_json).encode("utf-8")
        req = Request(MICROREACT_API_URL, data=body, headers=headers, method="POST")
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
        "--colour-palette",
        choices=list(COLOUR_PALETTES.keys()),
        default="default",
        help="Colour palette for Microreact columns (default: red/green)",
    )
    parser.add_argument(
        "--tree",
        default=None,
        help="Newick tree file to include in Microreact project",
    )
    parser.add_argument(
        "--output-microreact",
        default=None,
        help="Save a .microreact file for drag-and-drop upload to microreact.org",
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
    palette = args.colour_palette if add_colours else None
    write_csv(samples, drug_classes, gene_list, Path(args.output),
              add_colours=add_colours, colour_palette=palette)

    # Read tree if provided
    tree_data = None
    if args.tree:
        tree_path = Path(args.tree)
        if tree_path.is_file():
            tree_data = tree_path.read_text()
        else:
            print(f"Warning: Tree file {args.tree} not found", file=sys.stderr)

    # Save .microreact file
    if args.output_microreact:
        csv_data = build_csv_string(samples, drug_classes, gene_list,
                                    add_colours=add_colours, colour_palette=palette)
        project = build_microreact_project_json(csv_data, tree_data, args.project_name)
        with open(args.output_microreact, "w") as f:
            json.dump(project, f)
        print(f"Saved Microreact project to {args.output_microreact}", file=sys.stderr)

    # Microreact API upload
    api_key = args.microreact_api_key or os.environ.get("MICROREACT_API_KEY")
    if api_key:
        csv_data = build_csv_string(samples, drug_classes, gene_list,
                                    add_colours=add_colours, colour_palette=palette)
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
