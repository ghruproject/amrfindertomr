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
    COLOUR_ABSENT,
    COLOUR_HAS_GENES,
    COLOUR_NO_GENES,
    COLOUR_PRESENT,
    build_metadata,
    collect_inputs,
    detect_format,
    extract_sample_name,
    normalize_row,
    sanitize_column_name,
    write_csv,
)


# ---------------------------------------------------------------------------
# Fixtures: create minimal AMR tool-style files for testing
# ---------------------------------------------------------------------------

AMRFINDER_HEADER = (
    "Name\tProtein id\tContig id\tStart\tStop\tStrand\tElement symbol\t"
    "Element name\tScope\tType\tSubtype\tClass\tSubclass\tMethod\t"
    "Target length\tReference sequence length\t% Coverage of reference\t"
    "% Identity to reference\tAlignment length\tClosest reference accession\t"
    "Closest reference name\tHMM accession\tHMM description"
)

ABRICATE_HEADER = (
    "#FILE\tSEQUENCE\tSTART\tEND\tSTRAND\tGENE\tCOVERAGE\tCOVERAGE_MAP\t"
    "GAPS\t%COVERAGE\t%IDENTITY\tDATABASE\tACCESSION\tPRODUCT\tRESISTANCE"
)

RGI_HEADER = (
    "ORF_ID\tContig\tStart\tStop\tOrientation\tCut_Off\tPass_Bitscore\t"
    "Best_Hit_Bitscore\tBest_Hit_ARO\tBest_Hit_Identity\tARO\tModel_type\t"
    "SNPs_in_matching_region\tOther_SNPs\tDrug Class\tResistance Mechanism\t"
    "AMR Gene Family\tPredicted_DNA\tPredicted_Protein\t"
    "CARD_Protein_Sequence_Source\tHomologs"
)

RESFINDER_HEADER = (
    "Resistance gene\tIdentity\tAlignment Length/Gene Length\tCoverage\t"
    "Position in reference\tContig\tPosition in contig\tPhenotype\t"
    "Accession no."
)


def _make_amrfinder_row(name, gene, drug_class, element_type="AMR", subclass="NA",
                         coverage="100", identity="100", scope="core"):
    return (
        f"{name}\tNA\tcontig1\t100\t900\t+\t{gene}\tsome enzyme\t{scope}\t"
        f"{element_type}\t{element_type}\t{drug_class}\t{subclass}\tEXACTX\t"
        f"300\t300\t{coverage}\t{identity}\t300\tWP_000.1\tref name\tNA\tNA"
    )


def _make_abricate_row(filename, gene, resistance="AMR_DRUG", coverage="100.00",
                        identity="99.50"):
    return (
        f"{filename}\tcontig1\t100\t900\t+\t{gene}\t1-300/300\t========\t"
        f"0/0\t{coverage}\t{identity}\tncbi\tACC123\tsome product\t{resistance}"
    )


def _make_rgi_row(orf_id, gene, drug_class="beta-lactam", identity="99.5",
                   cut_off="Perfect"):
    return (
        f"{orf_id}\tcontig1\t100\t900\t+\t{cut_off}\t500\t600\t{gene}\t"
        f"{identity}\tARO:123\tProtein Homolog Model\tNA\tNA\t{drug_class}\t"
        f"antibiotic inactivation\tsome family\tATG\tMET\tNCBI\tNA"
    )


def _make_resfinder_row(gene, phenotype="Ampicillin", identity="100.00",
                          coverage="100"):
    return (
        f"{gene}\t{identity}\t300/300\t{coverage}\t1..900\tcontig1\t100..900\t"
        f"{phenotype}\tACC123"
    )


@pytest.fixture
def sample1_tsv(tmp_path):
    """AMRFinderPlus sample with 3 AMR genes across 3 drug classes."""
    p = tmp_path / "sample1_amr.tsv"
    lines = [
        AMRFINDER_HEADER,
        _make_amrfinder_row("sample1", "blaTEM-1", "BETA-LACTAM"),
        _make_amrfinder_row("sample1", "aadA1", "AMINOGLYCOSIDE"),
        _make_amrfinder_row("sample1", "sul1", "SULFONAMIDE"),
    ]
    p.write_text("\n".join(lines) + "\n")
    return p


@pytest.fixture
def sample2_tsv(tmp_path):
    """AMRFinderPlus sample with 2 AMR genes + 1 STRESS gene."""
    p = tmp_path / "sample2_amr.tsv"
    lines = [
        AMRFINDER_HEADER,
        _make_amrfinder_row("sample2", "blaTEM-1", "BETA-LACTAM"),
        _make_amrfinder_row("sample2", "tet(A)", "TETRACYCLINE"),
        _make_amrfinder_row("sample2", "silE", "SILVER", element_type="STRESS"),
    ]
    p.write_text("\n".join(lines) + "\n")
    return p


@pytest.fixture
def empty_tsv(tmp_path):
    """AMRFinderPlus output with header only (no hits)."""
    p = tmp_path / "empty_amr.tsv"
    p.write_text(AMRFINDER_HEADER + "\n")
    return p


@pytest.fixture
def sample_with_na_gene(tmp_path):
    """Sample where Element symbol is NA (should be skipped)."""
    p = tmp_path / "na_gene_amr.tsv"
    lines = [
        AMRFINDER_HEADER,
        _make_amrfinder_row("sampleNA", "NA", "BETA-LACTAM"),
        _make_amrfinder_row("sampleNA", "blaTEM-1", "BETA-LACTAM"),
    ]
    p.write_text("\n".join(lines) + "\n")
    return p


@pytest.fixture
def abricate_tsv(tmp_path):
    """ABRicate output file."""
    p = tmp_path / "sample_abricate.tsv"
    lines = [
        ABRICATE_HEADER,
        _make_abricate_row("sample_ab.fasta", "blaTEM-1", "BETA-LACTAM"),
        _make_abricate_row("sample_ab.fasta", "sul1", "SULFONAMIDE"),
    ]
    p.write_text("\n".join(lines) + "\n")
    return p


@pytest.fixture
def rgi_tsv(tmp_path):
    """CARD RGI output file."""
    p = tmp_path / "sample_rgi.txt"
    lines = [
        RGI_HEADER,
        _make_rgi_row("sample_rgi_NODE_1_len_5000 orf1", "TEM-1", "beta-lactam"),
        _make_rgi_row("sample_rgi_NODE_2_len_3000 orf2", "OXA-1", "beta-lactam"),
    ]
    p.write_text("\n".join(lines) + "\n")
    return p


@pytest.fixture
def resfinder_tsv(tmp_path):
    """ResFinder output file."""
    p = tmp_path / "sample_resfinder.tsv"
    lines = [
        RESFINDER_HEADER,
        _make_resfinder_row("blaTEM-1", "Ampicillin"),
        _make_resfinder_row("sul1", "Sulfamethoxazole"),
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


class TestDetectFormat:
    def test_amrfinderplus(self, sample1_tsv):
        assert detect_format(sample1_tsv) == "amrfinderplus"

    def test_abricate(self, abricate_tsv):
        assert detect_format(abricate_tsv) == "abricate"

    def test_rgi(self, rgi_tsv):
        assert detect_format(rgi_tsv) == "rgi"

    def test_resfinder(self, resfinder_tsv):
        assert detect_format(resfinder_tsv) == "resfinder"

    def test_unknown(self, tmp_path):
        p = tmp_path / "random.tsv"
        p.write_text("col1\tcol2\tcol3\n")
        assert detect_format(p) == "unknown"


class TestExtractSampleName:
    def test_from_name_column(self, sample1_tsv):
        from amr2microreact import parse_tsv
        rows, _ = parse_tsv(sample1_tsv)
        assert extract_sample_name(sample1_tsv, rows) == "sample1"

    def test_fallback_to_filename(self, tmp_path):
        p = tmp_path / "my_sample_amr.tsv"
        p.write_text(AMRFINDER_HEADER + "\n")
        from amr2microreact import parse_tsv
        rows, _ = parse_tsv(p)
        assert extract_sample_name(p, rows) == "my_sample"

    def test_fallback_strips_result(self, tmp_path):
        p = tmp_path / "ESC_AC7680AA_AS.result.tsv"
        p.write_text(AMRFINDER_HEADER + "\n")
        from amr2microreact import parse_tsv
        rows, _ = parse_tsv(p)
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
        assert len(files) == 3


class TestBuildMetadata:
    def test_basic(self, sample1_tsv):
        samples, classes, genes, fmts = build_metadata([sample1_tsv])
        assert "sample1" in samples
        assert sorted(classes) == ["AMINOGLYCOSIDE", "BETA-LACTAM", "SULFONAMIDE"]
        assert sorted(genes) == ["aadA1", "blaTEM-1", "sul1"]
        assert "amrfinderplus" in fmts

    def test_two_samples_share_gene(self, sample1_tsv, sample2_tsv):
        samples, classes, genes, _ = build_metadata([sample1_tsv, sample2_tsv])
        assert len(samples) == 2
        assert "blaTEM-1" in genes
        assert "blaTEM-1" in samples["sample1"]["genes"]
        assert "blaTEM-1" in samples["sample2"]["genes"]

    def test_amr_only_filter(self, sample2_tsv):
        samples, classes, genes, _ = build_metadata([sample2_tsv], amr_only=True)
        assert "silE" not in genes
        assert "SILVER" not in classes
        assert "blaTEM-1" in genes
        assert "tet(A)" in genes

    def test_amr_only_false_includes_stress(self, sample2_tsv):
        samples, classes, genes, _ = build_metadata([sample2_tsv], amr_only=False)
        assert "silE" in genes
        assert "SILVER" in classes

    def test_na_gene_skipped(self, sample_with_na_gene):
        samples, classes, genes, _ = build_metadata([sample_with_na_gene])
        assert "NA" not in genes
        assert "blaTEM-1" in genes

    def test_empty_file_skipped(self, empty_tsv, sample1_tsv):
        samples, classes, genes, _ = build_metadata([empty_tsv, sample1_tsv])
        assert len(samples) == 1


# ---------------------------------------------------------------------------
# Tests for new features: colours, filtering, multi-format
# ---------------------------------------------------------------------------


class TestColourColumns:
    def test_colour_columns_present(self, sample1_tsv, tmp_path):
        samples, classes, genes, _ = build_metadata([sample1_tsv])
        out = tmp_path / "out.csv"
        write_csv(samples, classes, genes, out, add_colours=True)

        with open(out) as f:
            header = f.readline().strip().split(",")

        colour_cols = [h for h in header if "__colour" in h]
        # Should have one colour col per drug class + one per gene
        assert len(colour_cols) == len(classes) + len(genes)

    def test_colour_values_correct(self, sample1_tsv, sample2_tsv, tmp_path):
        samples, classes, genes, _ = build_metadata([sample1_tsv, sample2_tsv])
        out = tmp_path / "out.csv"
        write_csv(samples, classes, genes, out, add_colours=True)

        with open(out) as f:
            reader = csv.DictReader(f)
            rows = {r["id"]: r for r in reader}

        # sample1 has BETA-LACTAM genes -> red
        bl_col = sanitize_column_name("BETA-LACTAM") + "__colour"
        assert rows["sample1"][bl_col] == COLOUR_HAS_GENES

        # sample1 has no TETRACYCLINE -> green
        tet_col = sanitize_column_name("TETRACYCLINE") + "__colour"
        assert rows["sample1"][tet_col] == COLOUR_NO_GENES

        # sample1 has blaTEM-1 -> red
        gene_col = sanitize_column_name("blaTEM-1") + "__colour"
        assert rows["sample1"][gene_col] == COLOUR_PRESENT

        # sample1 has no tet(A) -> grey
        teta_col = sanitize_column_name("tet(A)") + "__colour"
        assert rows["sample1"][teta_col] == COLOUR_ABSENT

    def test_no_colours_flag(self, sample1_tsv, tmp_path):
        samples, classes, genes, _ = build_metadata([sample1_tsv])
        out = tmp_path / "out.csv"
        write_csv(samples, classes, genes, out, add_colours=False)

        with open(out) as f:
            header = f.readline().strip().split(",")

        colour_cols = [h for h in header if "__colour" in h]
        assert len(colour_cols) == 0


class TestFilteringOptions:
    def test_min_coverage_filter(self, tmp_path):
        p = tmp_path / "cov_amr.tsv"
        lines = [
            AMRFINDER_HEADER,
            _make_amrfinder_row("s1", "blaTEM-1", "BETA-LACTAM", coverage="95"),
            _make_amrfinder_row("s1", "aadA1", "AMINOGLYCOSIDE", coverage="50"),
        ]
        p.write_text("\n".join(lines) + "\n")

        samples, _, genes, _ = build_metadata([p], min_coverage=80.0)
        assert "blaTEM-1" in genes
        assert "aadA1" not in genes

    def test_min_identity_filter(self, tmp_path):
        p = tmp_path / "id_amr.tsv"
        lines = [
            AMRFINDER_HEADER,
            _make_amrfinder_row("s1", "blaTEM-1", "BETA-LACTAM", identity="99"),
            _make_amrfinder_row("s1", "aadA1", "AMINOGLYCOSIDE", identity="70"),
        ]
        p.write_text("\n".join(lines) + "\n")

        samples, _, genes, _ = build_metadata([p], min_identity=90.0)
        assert "blaTEM-1" in genes
        assert "aadA1" not in genes

    def test_scope_filter_core(self, tmp_path):
        p = tmp_path / "scope_amr.tsv"
        lines = [
            AMRFINDER_HEADER,
            _make_amrfinder_row("s1", "blaTEM-1", "BETA-LACTAM", scope="core"),
            _make_amrfinder_row("s1", "mexA", "EFFLUX", scope="plus"),
        ]
        p.write_text("\n".join(lines) + "\n")

        samples, _, genes, _ = build_metadata([p], scope_filter="core")
        assert "blaTEM-1" in genes
        assert "mexA" not in genes

    def test_class_filter(self, sample1_tsv):
        samples, classes, genes, _ = build_metadata(
            [sample1_tsv], class_filter={"BETA-LACTAM"}
        )
        assert "BETA-LACTAM" in classes
        assert "AMINOGLYCOSIDE" not in classes
        # aadA1 has class AMINOGLYCOSIDE, filtered out
        assert "aadA1" not in genes


class TestMultiFormatInput:
    def test_abricate_parsing(self, abricate_tsv):
        samples, classes, genes, fmts = build_metadata([abricate_tsv])
        assert "abricate" in fmts
        assert "blaTEM-1" in genes
        assert "sul1" in genes
        assert len(samples) == 1

    def test_rgi_parsing(self, rgi_tsv):
        samples, classes, genes, fmts = build_metadata([rgi_tsv])
        assert "rgi" in fmts
        assert "TEM-1" in genes
        assert "OXA-1" in genes

    def test_resfinder_parsing(self, resfinder_tsv):
        samples, classes, genes, fmts = build_metadata([resfinder_tsv])
        assert "resfinder" in fmts
        assert "blaTEM-1" in genes
        assert "sul1" in genes

    def test_mixed_formats(self, sample1_tsv, abricate_tsv):
        samples, _, _, fmts = build_metadata([sample1_tsv, abricate_tsv])
        assert "amrfinderplus" in fmts
        assert "abricate" in fmts
        assert len(samples) >= 2

    def test_unknown_format_skipped(self, tmp_path, sample1_tsv):
        unknown = tmp_path / "random.tsv"
        unknown.write_text("col1\tcol2\tcol3\nval1\tval2\tval3\n")
        samples, _, _, fmts = build_metadata([unknown, sample1_tsv])
        assert "unknown" not in fmts
        assert len(samples) == 1


class TestWriteCsv:
    def test_output_structure_no_colours(self, sample1_tsv, tmp_path):
        samples, classes, genes, _ = build_metadata([sample1_tsv])
        out = tmp_path / "out.csv"
        write_csv(samples, classes, genes, out, add_colours=False)

        with open(out) as f:
            reader = csv.reader(f)
            header = next(reader)
            rows = list(reader)

        assert header[0] == "id"
        assert len(header) == 1 + len(classes) + len(genes)
        assert len(rows) == 1
        assert rows[0][0] == "sample1"

    def test_gene_presence_values(self, sample1_tsv, sample2_tsv, tmp_path):
        samples, classes, genes, _ = build_metadata([sample1_tsv, sample2_tsv])
        out = tmp_path / "out.csv"
        write_csv(samples, classes, genes, out, add_colours=False)

        with open(out) as f:
            reader = csv.DictReader(f)
            rows = {r["id"]: r for r in reader}

        assert rows["sample1"][sanitize_column_name("sul1")] == "yes"
        assert rows["sample2"][sanitize_column_name("sul1")] == "no"
        assert rows["sample1"][sanitize_column_name("blaTEM-1")] == "yes"
        assert rows["sample2"][sanitize_column_name("blaTEM-1")] == "yes"

    def test_drug_class_summary(self, sample1_tsv, tmp_path):
        samples, classes, genes, _ = build_metadata([sample1_tsv])
        out = tmp_path / "out.csv"
        write_csv(samples, classes, genes, out, add_colours=False)

        with open(out) as f:
            reader = csv.DictReader(f)
            row = next(reader)

        assert row[sanitize_column_name("BETA-LACTAM")] == "blaTEM-1"
        assert row[sanitize_column_name("AMINOGLYCOSIDE")] == "aadA1"

    def test_na_for_missing_class(self, sample1_tsv, sample2_tsv, tmp_path):
        samples, classes, genes, _ = build_metadata([sample1_tsv, sample2_tsv])
        out = tmp_path / "out.csv"
        write_csv(samples, classes, genes, out, add_colours=False)

        with open(out) as f:
            reader = csv.DictReader(f)
            rows = {r["id"]: r for r in reader}

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
        samples, classes, genes, fmts = build_metadata(files)
        assert len(samples) == 20
        assert "amrfinderplus" in fmts

    def test_output_csv_valid(self, tmp_path):
        files = collect_inputs([str(REAL_AMR_DIR)])
        samples, classes, genes, _ = build_metadata(files)
        out = tmp_path / "acceptance.csv"
        write_csv(samples, classes, genes, out, add_colours=False)

        with open(out) as f:
            reader = csv.DictReader(f)
            rows = list(reader)

        assert len(rows) == 20
        for row in rows:
            assert row["id"]
            assert row["id"].startswith("ESC_")

    def test_gene_values_only_yes_or_no(self, tmp_path):
        files = collect_inputs([str(REAL_AMR_DIR)])
        samples, classes, genes, _ = build_metadata(files)
        out = tmp_path / "acceptance.csv"
        write_csv(samples, classes, genes, out, add_colours=False)

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
        if not REAL_TREE.exists():
            pytest.skip("Test tree not available")

        files = collect_inputs([str(REAL_AMR_DIR)])
        samples, classes, genes, _ = build_metadata(files)
        out = tmp_path / "acceptance.csv"
        write_csv(samples, classes, genes, out, add_colours=False)

        with open(out) as f:
            reader = csv.DictReader(f)
            csv_ids = {row["id"] for row in reader}

        tree_text = REAL_TREE.read_text()
        import re
        tips = set(re.findall(r"[A-Z_0-9]+(?=:)", tree_text))

        assert tips == csv_ids

    def test_amr_only_has_fewer_genes(self):
        files = collect_inputs([str(REAL_AMR_DIR)])
        _, _, genes_all, _ = build_metadata(files, amr_only=False)
        _, _, genes_amr, _ = build_metadata(files, amr_only=True)
        assert len(genes_amr) <= len(genes_all)

    def test_with_colours(self, tmp_path):
        files = collect_inputs([str(REAL_AMR_DIR)])
        samples, classes, genes, _ = build_metadata(files)
        out = tmp_path / "acceptance_colour.csv"
        write_csv(samples, classes, genes, out, add_colours=True)

        with open(out) as f:
            header = f.readline().strip().split(",")

        colour_cols = [h for h in header if "__colour" in h]
        assert len(colour_cols) == len(classes) + len(genes)


# ---------------------------------------------------------------------------
# Regression tests
# ---------------------------------------------------------------------------


class TestRegressionMultipleGenesPerClass:
    def test_multiple_genes_comma_separated(self, tmp_path):
        p = tmp_path / "multi_amr.tsv"
        lines = [
            AMRFINDER_HEADER,
            _make_amrfinder_row("sampleM", "blaTEM-1", "BETA-LACTAM"),
            _make_amrfinder_row("sampleM", "blaCTX-M-15", "BETA-LACTAM"),
            _make_amrfinder_row("sampleM", "blaOXA-1", "BETA-LACTAM"),
        ]
        p.write_text("\n".join(lines) + "\n")

        samples, classes, genes, _ = build_metadata([p])
        out = tmp_path / "out.csv"
        write_csv(samples, classes, genes, out, add_colours=False)

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
    def test_no_duplicate_genes(self, tmp_path):
        p = tmp_path / "dup_amr.tsv"
        lines = [
            AMRFINDER_HEADER,
            _make_amrfinder_row("sampleD", "blaTEM-1", "BETA-LACTAM"),
            _make_amrfinder_row("sampleD", "blaTEM-1", "BETA-LACTAM"),
        ]
        p.write_text("\n".join(lines) + "\n")

        samples, classes, genes, _ = build_metadata([p])
        assert genes.count("blaTEM-1") == 1

        out = tmp_path / "out.csv"
        write_csv(samples, classes, genes, out, add_colours=False)

        with open(out) as f:
            reader = csv.DictReader(f)
            row = next(reader)

        assert row[sanitize_column_name("BETA-LACTAM")] == "blaTEM-1"


class TestRegressionSortedOutput:
    def test_samples_sorted(self, tmp_path):
        for name in ["sampleZ", "sampleA", "sampleM"]:
            p = tmp_path / f"{name}_amr.tsv"
            lines = [AMRFINDER_HEADER, _make_amrfinder_row(name, "blaTEM-1", "BETA-LACTAM")]
            p.write_text("\n".join(lines) + "\n")

        files = collect_inputs([str(tmp_path)])
        samples, classes, genes, _ = build_metadata(files)
        out = tmp_path / "out.csv"
        write_csv(samples, classes, genes, out, add_colours=False)

        with open(out) as f:
            reader = csv.DictReader(f)
            ids = [row["id"] for row in reader]

        assert ids == ["sampleA", "sampleM", "sampleZ"]
