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

        # ---------------------------------------------------------------
        # Tabs: Data | Visualisation | Export to Microreact
        # ---------------------------------------------------------------
        tab_data, tab_viz, tab_export = st.tabs(
            ["📋 Data", "📊 Visualisation", "🚀 Export to Microreact"]
        )

        # --- Tab 1: Data ---
        with tab_data:
            st.subheader("Metadata Preview")
            df = pd.read_csv(io.StringIO(csv_data))
            display_cols = [c for c in df.columns if "__colour" not in c]
            st.dataframe(df[display_cols], use_container_width=True, height=400)

            st.download_button(
                label="Download Microreact CSV",
                data=csv_data,
                file_name="microreact_metadata.csv",
                mime="text/csv",
            )

            with st.expander("Drug class breakdown"):
                for dc in drug_classes:
                    genes_in_dc = set()
                    for s in samples.values():
                        genes_in_dc.update(s["classes"].get(dc, set()))
                    st.markdown(f"**{dc}**: {', '.join(sorted(genes_in_dc))}")

        # --- Tab 2: Visualisation ---
        with tab_viz:
            # Build a clean dataframe for viz (no colour columns)
            df_clean = df[display_cols].set_index("id")

            # Separate drug class summary cols from gene presence cols
            dc_sanitized = [sanitize_column_name(dc) for dc in drug_classes]
            gene_sanitized = [sanitize_column_name(g) for g in gene_list]
            gene_cols_in_df = [c for c in gene_sanitized if c in df_clean.columns]
            dc_cols_in_df = [c for c in dc_sanitized if c in df_clean.columns]

            # --- Heatmap: Gene presence/absence ---
            if gene_cols_in_df:
                st.subheader("Gene Presence / Absence Heatmap")
                heatmap_df = df_clean[gene_cols_in_df].replace({"yes": 1, "no": 0})

                # Use Streamlit's built-in charting - show as a styled dataframe
                styled = heatmap_df.style.map(
                    lambda v: f"background-color: {COLOUR_PRESENT}; color: white"
                    if v == 1
                    else f"background-color: {COLOUR_ABSENT}; color: #999"
                )
                st.dataframe(styled, use_container_width=True, height=400)

            # --- Bar chart: Genes per drug class ---
            if drug_classes:
                st.subheader("Resistance Genes per Drug Class")
                class_counts = {}
                for dc in drug_classes:
                    genes_in_dc = set()
                    for s in samples.values():
                        genes_in_dc.update(s["classes"].get(dc, set()))
                    class_counts[dc] = len(genes_in_dc)

                class_df = pd.DataFrame(
                    {"Drug Class": list(class_counts.keys()),
                     "Unique Genes": list(class_counts.values())}
                ).sort_values("Unique Genes", ascending=False)
                st.bar_chart(class_df, x="Drug Class", y="Unique Genes")

            # --- Bar chart: Genes per sample ---
            st.subheader("Resistance Genes per Sample")
            sample_counts = {
                name: len(sdata["genes"]) for name, sdata in samples.items()
            }
            sample_df = pd.DataFrame(
                {"Sample": list(sample_counts.keys()),
                 "Total Genes": list(sample_counts.values())}
            ).sort_values("Total Genes", ascending=False)
            st.bar_chart(sample_df, x="Sample", y="Total Genes")

            # --- Drug class coverage across samples ---
            if dc_cols_in_df:
                st.subheader("Drug Class Resistance Across Samples")
                # Count how many samples have at least one gene per class
                class_sample_counts = {}
                for dc in drug_classes:
                    count = sum(
                        1 for s in samples.values()
                        if s["classes"].get(dc, set())
                    )
                    class_sample_counts[dc] = count

                coverage_df = pd.DataFrame(
                    {"Drug Class": list(class_sample_counts.keys()),
                     "Samples with Resistance": list(class_sample_counts.values())}
                ).sort_values("Samples with Resistance", ascending=False)
                st.bar_chart(coverage_df, x="Drug Class", y="Samples with Resistance")

        # --- Tab 3: Export to Microreact ---
        with tab_export:
            st.subheader("Create a Microreact Project")
            st.warning(
                "**Experimental** — API project creation may not work reliably due to "
                "Cloudflare restrictions on Streamlit Cloud. If it fails, download the CSV "
                "and tree and upload them manually at [microreact.org](https://microreact.org)."
            )
            st.markdown(
                "Get your API access token from "
                "[your Microreact account settings](https://microreact.org/my-account/settings)."
            )

            # localStorage persistence for API key
            _LS_KEY = "amr2microreact_api_key"
            stored_key = st.query_params.get("_mr_key", "")

            api_key = st.text_input(
                "Microreact API Access Token",
                value=stored_key,
                type="password",
                help="Stored in your browser only — never sent to our server.",
            )

            components.html(
                f"""
                <script>
                const KEY = "{_LS_KEY}";
                const currentInput = "{api_key}";
                const stored = localStorage.getItem(KEY) || "";
                const urlParams = new URLSearchParams(window.location.search);
                const paramKey = urlParams.get("_mr_key") || "";
                if (stored && !paramKey) {{
                    urlParams.set("_mr_key", stored);
                    const newUrl = window.location.pathname + "?" + urlParams.toString();
                    window.history.replaceState(null, "", newUrl);
                    window.location.reload();
                }}
                if (currentInput && currentInput !== stored) {{
                    localStorage.setItem(KEY, currentInput);
                }}
                </script>
                """,
                height=0,
            )

            tree_file = st.file_uploader(
                "Upload Newick tree (optional)",
                type=["nwk", "newick", "tre", "tree", "nhx", "treefile"],
                help="Tip labels must match the id column in the metadata.",
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
