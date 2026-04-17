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
import streamlit.components.v1 as components

from amr2microreact import (
    COLOUR_ABSENT,
    COLOUR_HAS_GENES,
    COLOUR_NO_GENES,
    COLOUR_PRESENT,
    build_csv_string,
    build_metadata,
    sanitize_column_name,
    upload_to_microreact,
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

        # Build CSV
        csv_data = build_csv_string(samples, drug_classes, gene_list, add_colours=add_colours)

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

        # Download CSV
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

        # -----------------------------------------------------------
        # Export to Microreact
        # -----------------------------------------------------------
        st.divider()
        st.subheader("Export to Microreact")
        st.markdown(
            "Create a Microreact project directly. You need an API access token from "
            "[your Microreact account settings](https://microreact.org/my-account)."
        )

        # Load API key from localStorage via JS
        # This injects a hidden iframe that reads localStorage and writes to a query param
        _LS_KEY = "amr2microreact_api_key"

        # Read stored key from query params (set by JS below)
        stored_key = st.query_params.get("_mr_key", "")

        api_key = st.text_input(
            "Microreact API Access Token",
            value=stored_key,
            type="password",
            help="Your token is stored in your browser's local storage only — never sent to our server.",
        )

        # JS component to persist API key to localStorage and read it back
        components.html(
            f"""
            <script>
            const KEY = "{_LS_KEY}";
            const currentInput = "{api_key}";

            // On load: if we have a stored key and no query param, set it
            const stored = localStorage.getItem(KEY) || "";
            const urlParams = new URLSearchParams(window.location.search);
            const paramKey = urlParams.get("_mr_key") || "";

            if (stored && !paramKey) {{
                // Push stored key into query params so Streamlit can read it
                urlParams.set("_mr_key", stored);
                const newUrl = window.location.pathname + "?" + urlParams.toString();
                window.history.replaceState(null, "", newUrl);
                window.location.reload();
            }}

            // If user entered a key, save it
            if (currentInput && currentInput !== stored) {{
                localStorage.setItem(KEY, currentInput);
            }}
            </script>
            """,
            height=0,
        )

        # Tree upload
        tree_file = st.file_uploader(
            "Upload Newick tree (optional)",
            type=["nwk", "newick", "tre", "tree", "nhx", "treefile"],
            help="Upload a phylogenetic tree. Tip labels must match the id column.",
        )

        project_name = st.text_input(
            "Project name",
            value="AMR Microreact Project",
        )

        if st.button("Create Microreact Project", type="primary", disabled=not api_key):
            tree_data = None
            if tree_file:
                tree_data = tree_file.getvalue().decode("utf-8")

            with st.spinner("Uploading to Microreact..."):
                try:
                    result = upload_to_microreact(
                        csv_data, api_key, tree_data, project_name
                    )
                    url = result.get("url", "")
                    if url:
                        st.success(f"Project created! [Open in Microreact]({url})")
                        st.balloons()
                    else:
                        st.success(f"Project created! Response: {result}")
                except RuntimeError as e:
                    st.error(f"Failed: {e}")

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
