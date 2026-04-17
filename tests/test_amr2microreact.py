"""Tests for amr2microreact.py"""

import csv
import os
import tempfile
from pathlib import Path

import pytest

# Add parent to path so we can import the module
import sys

sys.path.insert(0, str(Path(__file__).parent.parent))

from amr2microreact import (
    build_metadata,
    collect_inputs,
    extract_sample_name,
    parse_amrfinder,
    sanitize_column_name,
    write_csv,
)


# ---------------------------------------------------------------------------
# Fixtures: create minimal AMRFinderPlus-style TSV files for testing
# ---------------------------------------------------------------------------

HEADER = (
    "Name\tProtein id\tContig id\tStart\tStop\tStrand\tElement symbol\t"
    "Element name\tScope\tType\tSubtype\tClass\tSubclass\tMethod\t"
    "Target length\tReference sequence length\t% Coverage of reference\t"
    "% Identity to reference\tAlignment length\tClosest reference accession\t"
    "Closest reference name\tHMM accession\tHMM description"
)


def _make_row(name, gene, drug_class, element_type="AMR", subclass="NA"):
    return (
        f"{name}\tNA\tcontig1\t100\t900\t+\t{gene}\tsome enzyme\tcore\t"
        f"{element_type}\t{element_type}\t{drug_class}\t{subclass}\tEXACTX\t"
        f"300\t300\t100\t100\t300\tWP_000.1\tref name\tNA\tNA"
    )


@pytest.fixture
def sample1_tsv(tmp_path):
    """Sample with 3 AMR genes across 2 drug classes."""
    p = tmp_path / "sample1_amr.tsv"
    lines = [
        HEADER,
        _make_row("sample1", "blaTEM-1", "BETA-LACTAM"),
        _make_row("sample1", "aadA1", "AMINOGLYCOSIDE"),
        _make_row("sample1", "sul1", "SULFONAMIDE"),
    ]
    p.write_text("\n".join(lines) + "\n")
    return p


@pytest.fixture
def sample2_tsv(tmp_path):
    """Sample with 2 AMR genes + 1 STRESS gene."""
    p = tmp_path / "sample2_amr.tsv"
    lines = [
        HEADER,
        _make_row("sample2", "blaTEM-1", "BETA-LACTAM"),
        _make_row("sample2", "tet(A)", "TETRACYCLINE"),
        _make_row("sample2", "silE", "SILVER", element_type="STRESS"),
    ]
    p.write_text("\n".join(lines) + "\n")
    return p


@pytest.fixture
def empty_tsv(tmp_path):
    """AMRFinderPlus output with header only (no hits)."""
    p = tmp_path / "empty_amr.tsv"
    p.write_text(HEADER + "\n")
    return p


@pytest.fixture
def sample_with_na_gene(tmp_path):
    """Sample where Element symbol is NA (should be skipped)."""
    p = tmp_path / "na_gene_amr.tsv"
    lines = [
        HEADER,
        _make_row("sampleNA", "NA", "BETA-LACTAM"),
        _make_row("sampleNA", "blaTEM-1", "BETA-LACTAM"),
    ]
    p.write_text("\n".join(lines) + "\n")
    return p


@pytest.fixture
def two_sample_dir(tmp_path, sample1_tsv, sample2_tsv):
    """Directory containing two sample TSVs."""
    return tmp_path


# ---------------------------------------------------------------------------
# Unit tests
# ---------------------------------------------------------------------------


class TestSanitizeColumnName:
    def test_parentheses(self):
        assert sanitize_column_name("aph(3'')-Ib") == "aph.3....Ib"

    def test_slash(self):
        assert sanitize_column_name("PHENICOL/QUINOLONE") == "PHENICOL.QUINOLONE"

    def test_dash(self):
        assert sanitize_column_name("BETA-LACTAM") == "BETA.LACTAM"

    def test_clean_name(self):
        assert sanitize_column_name("sul1") == "sul1"

    def test_quotes(self):
        assert sanitize_column_name("gene'name") == "gene.name"


class TestParseAmrfinder:
    def test_parses_rows(self, sample1_tsv):
        rows = parse_amrfinder(sample1_tsv)
        assert len(rows) == 3
        assert rows[0]["Element symbol"] == "blaTEM-1"
        assert rows[0]["Name"] == "sample1"

    def test_empty_file(self, empty_tsv):
        rows = parse_amrfinder(empty_tsv)
        assert rows == []


class TestExtractSampleName:
    def test_from_name_column(self, sample1_tsv):
        rows = parse_amrfinder(sample1_tsv)
        assert extract_sample_name(sample1_tsv, rows) == "sample1"

    def test_fallback_to_filename(self, tmp_path):
        p = tmp_path / "my_sample_amr.tsv"
        p.write_text(HEADER + "\n")
        rows = parse_amrfinder(p)
        assert extract_sample_name(p, rows) == "my_sample"

    def test_fallback_strips_result(self, tmp_path):
        p = tmp_path / "ESC_AC7680AA_AS.result.tsv"
        p.write_text(HEADER + "\n")
        rows = parse_amrfinder(p)
        assert extract_sample_name(p, rows) == "ESC_AC7680AA_AS"


class TestCollectInputs:
    def test_single_file(self, sample1_tsv):
        files = collect_inputs([str(sample1_tsv)])
        assert len(files) == 1

    def test_directory(self, two_sample_dir):
        files = collect_inputs([str(two_sample_dir)])
        assert len(files) == 2

    def test_missing_path(self, capsys):
        files = collect_inputs(["/nonexistent/path"])
        assert files == []

    def test_mixed_inputs(self, sample1_tsv, two_sample_dir):
        files = collect_inputs([str(sample1_tsv), str(two_sample_dir)])
        # sample1 appears both as direct file and in directory
        assert len(files) == 3


class TestBuildMetadata:
    def test_basic(self, sample1_tsv):
        samples, classes, genes = build_metadata([sample1_tsv])
        assert "sample1" in samples
        assert sorted(classes) == ["AMINOGLYCOSIDE", "BETA-LACTAM", "SULFONAMIDE"]
        assert sorted(genes) == ["aadA1", "blaTEM-1", "sul1"]

    def test_two_samples_share_gene(self, sample1_tsv, sample2_tsv):
        samples, classes, genes = build_metadata([sample1_tsv, sample2_tsv])
        assert len(samples) == 2
        assert "blaTEM-1" in genes  # shared gene
        # sample1 has blaTEM-1
        assert "blaTEM-1" in samples["sample1"]["genes"]
        # sample2 also has blaTEM-1
        assert "blaTEM-1" in samples["sample2"]["genes"]

    def test_amr_only_filter(self, sample2_tsv):
        samples, classes, genes = build_metadata([sample2_tsv], amr_only=True)
        # silE is STRESS, should be excluded
        assert "silE" not in genes
        assert "SILVER" not in classes
        assert "blaTEM-1" in genes
        assert "tet(A)" in genes

    def test_amr_only_false_includes_stress(self, sample2_tsv):
        samples, classes, genes = build_metadata([sample2_tsv], amr_only=False)
        assert "silE" in genes
        assert "SILVER" in classes

    def test_na_gene_skipped(self, sample_with_na_gene):
        samples, classes, genes = build_metadata([sample_with_na_gene])
        assert "NA" not in genes
        assert "blaTEM-1" in genes

    def test_empty_file_skipped(self, empty_tsv, sample1_tsv):
        samples, classes, genes = build_metadata([empty_tsv, sample1_tsv])
        assert len(samples) == 1  # only sample1


class TestWriteCsv:
    def test_output_structure(self, sample1_tsv, tmp_path):
        samples, classes, genes = build_metadata([sample1_tsv])
        out = tmp_path / "out.csv"
        write_csv(samples, classes, genes, out)

        with open(out) as f:
            reader = csv.reader(f)
            header = next(reader)
            rows = list(reader)

        # Header: id + drug classes + genes
        assert header[0] == "id"
        assert len(header) == 1 + len(classes) + len(genes)
        assert len(rows) == 1
        assert rows[0][0] == "sample1"

    def test_gene_presence_values(self, sample1_tsv, sample2_tsv, tmp_path):
        samples, classes, genes = build_metadata([sample1_tsv, sample2_tsv])
        out = tmp_path / "out.csv"
        write_csv(samples, classes, genes, out)

        with open(out) as f:
            reader = csv.DictReader(f)
            rows = {r["id"]: r for r in reader}

        # sample1 has sul1, sample2 does not
        assert rows["sample1"][sanitize_column_name("sul1")] == "yes"
        assert rows["sample2"][sanitize_column_name("sul1")] == "no"

        # Both have blaTEM-1
        assert rows["sample1"][sanitize_column_name("blaTEM-1")] == "yes"
        assert rows["sample2"][sanitize_column_name("blaTEM-1")] == "yes"

    def test_drug_class_summary(self, sample1_tsv, tmp_path):
        samples, classes, genes = build_metadata([sample1_tsv])
        out = tmp_path / "out.csv"
        write_csv(samples, classes, genes, out)

        with open(out) as f:
            reader = csv.DictReader(f)
            row = next(reader)

        assert row[sanitize_column_name("BETA-LACTAM")] == "blaTEM-1"
        assert row[sanitize_column_name("AMINOGLYCOSIDE")] == "aadA1"

    def test_na_for_missing_class(self, sample1_tsv, sample2_tsv, tmp_path):
        samples, classes, genes = build_metadata([sample1_tsv, sample2_tsv])
        out = tmp_path / "out.csv"
        write_csv(samples, classes, genes, out)

        with open(out) as f:
            reader = csv.DictReader(f)
            rows = {r["id"]: r for r in reader}

        # sample1 has no TETRACYCLINE genes
        assert rows["sample1"][sanitize_column_name("TETRACYCLINE")] == "NA"


# ---------------------------------------------------------------------------
# Acceptance / integration tests (using real test data if available)
# ---------------------------------------------------------------------------

REAL_AMR_DIR = Path(__file__).parent.parent / "test_data" / "amr_results"
REAL_TREE = Path(__file__).parent.parent / "test_data" / "test_tree.nwk"


@pytest.mark.skipif(
    not REAL_AMR_DIR.exists(), reason="Real test data not available"
)
class TestAcceptanceWithRealData:
    def test_processes_all_samples(self):
        files = collect_inputs([str(REAL_AMR_DIR)])
        samples, classes, genes = build_metadata(files)
        assert len(samples) == 20

    def test_output_csv_valid(self, tmp_path):
        files = collect_inputs([str(REAL_AMR_DIR)])
        samples, classes, genes = build_metadata(files)
        out = tmp_path / "acceptance.csv"
        write_csv(samples, classes, genes, out)

        with open(out) as f:
            reader = csv.DictReader(f)
            rows = list(reader)

        assert len(rows) == 20
        # Every row should have an id
        for row in rows:
            assert row["id"]
            assert row["id"].startswith("ESC_")

    def test_gene_values_only_yes_or_no(self, tmp_path):
        files = collect_inputs([str(REAL_AMR_DIR)])
        samples, classes, genes = build_metadata(files)
        out = tmp_path / "acceptance.csv"
        write_csv(samples, classes, genes, out)

        with open(out) as f:
            reader = csv.DictReader(f)
            header = reader.fieldnames
            gene_cols = [
                h for h in header if h not in ["id"] + [sanitize_column_name(c) for c in classes]
            ]
            for row in reader:
                for col in gene_cols:
                    assert row[col] in ("yes", "no"), f"Bad value in {col}: {row[col]}"

    def test_tree_tips_match_csv_ids(self, tmp_path):
        """The tree tip labels must match the id column in the CSV."""
        if not REAL_TREE.exists():
            pytest.skip("Test tree not available")

        files = collect_inputs([str(REAL_AMR_DIR)])
        samples, classes, genes = build_metadata(files)
        out = tmp_path / "acceptance.csv"
        write_csv(samples, classes, genes, out)

        # Extract CSV ids
        with open(out) as f:
            reader = csv.DictReader(f)
            csv_ids = {row["id"] for row in reader}

        # Extract tree tip labels
        tree_text = REAL_TREE.read_text()
        # Remove branch lengths and punctuation to get labels
        import re
        tips = set(re.findall(r"[A-Z_0-9]+(?=:)", tree_text))

        assert tips == csv_ids, f"Mismatch: tree-only={tips - csv_ids}, csv-only={csv_ids - tips}"

    def test_amr_only_has_fewer_genes(self):
        files = collect_inputs([str(REAL_AMR_DIR)])
        _, _, genes_all = build_metadata(files, amr_only=False)
        _, _, genes_amr = build_metadata(files, amr_only=True)
        assert len(genes_amr) <= len(genes_all)


# ---------------------------------------------------------------------------
# Regression tests
# ---------------------------------------------------------------------------


class TestRegressionMultipleGenesPerClass:
    """Ensure drug class summary correctly joins multiple genes."""

    def test_multiple_genes_comma_separated(self, tmp_path):
        p = tmp_path / "multi_amr.tsv"
        lines = [
            HEADER,
            _make_row("sampleM", "blaTEM-1", "BETA-LACTAM"),
            _make_row("sampleM", "blaCTX-M-15", "BETA-LACTAM"),
            _make_row("sampleM", "blaOXA-1", "BETA-LACTAM"),
        ]
        p.write_text("\n".join(lines) + "\n")

        samples, classes, genes = build_metadata([p])
        out = tmp_path / "out.csv"
        write_csv(samples, classes, genes, out)

        with open(out) as f:
            reader = csv.DictReader(f)
            row = next(reader)

        bl_val = row[sanitize_column_name("BETA-LACTAM")]
        gene_list = bl_val.split(",")
        assert len(gene_list) == 3
        assert "blaTEM-1" in gene_list
        assert "blaCTX-M-15" in gene_list
        assert "blaOXA-1" in gene_list


class TestRegressionDuplicateGenes:
    """Same gene appearing twice in one sample should not duplicate."""

    def test_no_duplicate_genes(self, tmp_path):
        p = tmp_path / "dup_amr.tsv"
        lines = [
            HEADER,
            _make_row("sampleD", "blaTEM-1", "BETA-LACTAM"),
            _make_row("sampleD", "blaTEM-1", "BETA-LACTAM"),
        ]
        p.write_text("\n".join(lines) + "\n")

        samples, classes, genes = build_metadata([p])
        assert genes.count("blaTEM-1") == 1

        out = tmp_path / "out.csv"
        write_csv(samples, classes, genes, out)

        with open(out) as f:
            reader = csv.DictReader(f)
            row = next(reader)

        assert row[sanitize_column_name("BETA-LACTAM")] == "blaTEM-1"


class TestRegressionSortedOutput:
    """Output should be deterministic - sorted samples and columns."""

    def test_samples_sorted(self, tmp_path):
        # Create files with names that would sort differently
        for name in ["sampleZ", "sampleA", "sampleM"]:
            p = tmp_path / f"{name}_amr.tsv"
            lines = [HEADER, _make_row(name, "blaTEM-1", "BETA-LACTAM")]
            p.write_text("\n".join(lines) + "\n")

        files = collect_inputs([str(tmp_path)])
        samples, classes, genes = build_metadata(files)
        out = tmp_path / "out.csv"
        write_csv(samples, classes, genes, out)

        with open(out) as f:
            reader = csv.DictReader(f)
            ids = [row["id"] for row in reader]

        assert ids == ["sampleA", "sampleM", "sampleZ"]
