"""
Streamlit web app for converting AMR tool output to Microreact metadata.
Supports AMRFinderPlus, ABRicate, ResFinder, and CARD RGI (auto-detected).
"""

import csv
import io
import tempfile
from pathlib import Path

import pandas as pd
import streamlit as st

from amr2microreact import (
    COLOUR_ABSENT,
    COLOUR_HAS_GENES,
    COLOUR_NO_GENES,
    COLOUR_PRESENT,
    build_metadata,
    sanitize_column_name,
)

st.set_page_config(
    page_title="AMR → Microreact",
    page_icon="🧬",
    layout="wide",
)

# Logo + title
col_logo, col_title = st.columns([1, 5])
with col_logo:
    logo_path = Path(__file__).parent / "logo.svg"
    if logo_path.exists():
        st.image(str(logo_path), width=80)
with col_title:
    st.title("AMR → Microreact Converter")

st.markdown(
    "Upload AMR tool output files and get a "
    "[Microreact](https://microreact.org)-compatible metadata CSV. "
    "Supports **AMRFinderPlus**, **ABRicate**, **ResFinder**, and **CARD RGI** "
    "(format auto-detected)."
)

# --- File upload ---
uploaded_files = st.file_uploader(
    "Upload AMR output files",
    type=["tsv", "txt", "csv"],
    accept_multiple_files=True,
    help="Select one or more AMR tool output files. Format is auto-detected from column headers.",
)

# --- Sidebar filters ---
with st.sidebar:
    st.header("Filters")

    amr_only = st.checkbox(
        "AMR genes only",
        value=False,
        help="Exclude stress and virulence genes.",
    )

    min_coverage = st.slider(
        "Min coverage (%)",
        min_value=0.0,
        max_value=100.0,
        value=0.0,
        step=5.0,
        help="Minimum coverage threshold.",
    )

    min_identity = st.slider(
        "Min identity (%)",
        min_value=0.0,
        max_value=100.0,
        value=0.0,
        step=5.0,
        help="Minimum identity threshold.",
    )

    scope_filter = st.selectbox(
        "Scope (AMRFinderPlus only)",
        options=["all", "core", "plus"],
        index=0,
    )

    add_colours = st.checkbox(
        "Add __colour columns",
        value=True,
        help="Include Microreact colour columns (red=present, green=absent).",
    )


if uploaded_files:
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_paths = []
        for uf in uploaded_files:
            p = Path(tmpdir) / uf.name
            p.write_bytes(uf.getvalue())
            tmp_paths.append(p)

        samples, drug_classes, gene_list, formats = build_metadata(
            tmp_paths,
            amr_only=amr_only,
            min_coverage=min_coverage,
            min_identity=min_identity,
            scope_filter=scope_filter if scope_filter != "all" else None,
        )

    if not samples:
        st.warning("No AMR genes found in the uploaded files.")
    else:
        fmt_str = ", ".join(sorted(formats))
        st.success(
            f"Processed **{len(samples)}** samples from **{fmt_str}** — "
            f"**{len(drug_classes)}** drug classes, **{len(gene_list)}** genes detected."
        )

        # Build CSV in memory
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
                    row.append(COLOUR_HAS_GENES if genes_in_class else COLOUR_NO_GENES)
            for gene in gene_list:
                present = gene in sdata["genes"]
                row.append("yes" if present else "no")
                if add_colours:
                    row.append(COLOUR_PRESENT if present else COLOUR_ABSENT)
            rows_out.append(row)

        buf = io.StringIO()
        writer = csv.writer(buf)
        writer.writerow(header)
        writer.writerows(rows_out)
        csv_data = buf.getvalue()

        # Summary metrics
        col1, col2, col3 = st.columns(3)
        with col1:
            st.metric("Samples", len(samples))
        with col2:
            st.metric("Drug Classes", len(drug_classes))
        with col3:
            st.metric("Genes", len(gene_list))

        # Preview (hide colour columns for readability)
        st.subheader("Preview")
        df = pd.read_csv(io.StringIO(csv_data))
        display_cols = [c for c in df.columns if "__colour" not in c]
        st.dataframe(df[display_cols], use_container_width=True, height=400)

        # Download
        st.download_button(
            label="Download Microreact CSV",
            data=csv_data,
            file_name="microreact_metadata.csv",
            mime="text/csv",
        )

        # Drug class breakdown
        with st.expander("Drug class breakdown"):
            for dc in drug_classes:
                genes_in_dc = set()
                for s in samples.values():
                    genes_in_dc.update(s["classes"].get(dc, set()))
                st.markdown(f"**{dc}**: {', '.join(sorted(genes_in_dc))}")

else:
    st.info("Upload one or more AMR tool output files to get started.")
    with st.expander("Supported formats"):
        st.markdown("""
- **AMRFinderPlus** — `amrfinder -o sample.tsv`
- **ABRicate** — `abricate sample.fasta > sample.tsv`
- **ResFinder** — `ResFinder_results_tab.txt`
- **CARD RGI** — `rgi main -o sample` → `sample.txt`

Format is **auto-detected** from column headers. You can mix formats in a single upload.
        """)
