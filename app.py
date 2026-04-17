"""
Streamlit web app for converting AMRFinderPlus output to Microreact metadata.
"""

import csv
import io

import streamlit as st

from amr2microreact import build_metadata, parse_amrfinder, extract_sample_name, sanitize_column_name

st.set_page_config(
    page_title="AMRFinder → Microreact",
    page_icon="🧬",
    layout="wide",
)

st.title("AMRFinder → Microreact Converter")
st.markdown(
    "Upload [AMRFinderPlus](https://github.com/ncbi/amr) output files and get a "
    "[Microreact](https://microreact.org)-compatible metadata CSV."
)

uploaded_files = st.file_uploader(
    "Upload AMRFinderPlus TSV output files",
    type=["tsv", "txt"],
    accept_multiple_files=True,
    help="Select one or more AMRFinderPlus output files (tab-separated).",
)

amr_only = st.checkbox(
    "AMR genes only (exclude stress/virulence)",
    value=False,
    help="If checked, only AMR-type genes are included. Stress and virulence genes are excluded.",
)

if uploaded_files:
    # Parse all uploaded files using the existing logic
    from collections import defaultdict
    from pathlib import Path
    import tempfile
    import os

    # Write uploaded files to temp dir so we can reuse existing functions
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_paths = []
        for uf in uploaded_files:
            p = Path(tmpdir) / uf.name
            p.write_bytes(uf.getvalue())
            tmp_paths.append(p)

        samples, drug_classes, gene_list = build_metadata(tmp_paths, amr_only=amr_only)

    if not samples:
        st.warning("No AMR genes found in the uploaded files.")
    else:
        st.success(f"Processed **{len(samples)}** samples — "
                   f"**{len(drug_classes)}** drug classes, **{len(gene_list)}** genes detected.")

        # Build the CSV in memory
        header = ["id"]
        for dc in drug_classes:
            header.append(sanitize_column_name(dc))
        for gene in gene_list:
            header.append(sanitize_column_name(gene))

        rows_out = []
        for sample_name in sorted(samples.keys()):
            sdata = samples[sample_name]
            row = [sample_name]
            for dc in drug_classes:
                genes_in_class = sdata["classes"].get(dc, set())
                row.append(",".join(sorted(genes_in_class)) if genes_in_class else "NA")
            for gene in gene_list:
                row.append("yes" if gene in sdata["genes"] else "no")
            rows_out.append(row)

        # Write to string buffer
        buf = io.StringIO()
        writer = csv.writer(buf)
        writer.writerow(header)
        writer.writerows(rows_out)
        csv_data = buf.getvalue()

        # Preview
        st.subheader("Preview")
        import pandas as pd
        df = pd.read_csv(io.StringIO(csv_data))
        st.dataframe(df, use_container_width=True, height=400)

        # Download button
        st.download_button(
            label="Download Microreact CSV",
            data=csv_data,
            file_name="microreact_metadata.csv",
            mime="text/csv",
        )

        # Show summary stats
        col1, col2, col3 = st.columns(3)
        with col1:
            st.metric("Samples", len(samples))
        with col2:
            st.metric("Drug Classes", len(drug_classes))
        with col3:
            st.metric("Genes", len(gene_list))

        # Optional: show drug class breakdown
        with st.expander("Drug class breakdown"):
            for dc in drug_classes:
                genes_in_dc = set()
                for s in samples.values():
                    genes_in_dc.update(s["classes"].get(dc, set()))
                st.markdown(f"**{dc}**: {', '.join(sorted(genes_in_dc))}")
else:
    st.info("Upload one or more AMRFinderPlus output files to get started.")
    with st.expander("How to generate AMRFinderPlus output"):
        st.code(
            "amrfinder --nucleotide assembly.fasta --name sample_id -o sample_amr.tsv",
            language="bash",
        )
        st.markdown(
            "Run this for each of your assemblies, then upload all the `.tsv` files here."
        )
