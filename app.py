import streamlit as st
import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import plotly.express as px
import io

st.title("🧬 BioSeq Toolkit")
st.info(
    """
    👋 Welcome to the BioSeq Toolkit!

    This app allows you to:
    - 📂 Upload DNA sequences in FASTA format
    - 🔍 Calculate GC content
    - 🧬 Translate DNA to protein
    - 📊 Visualize base composition as a bar chart
    """
)

# Sidebar options
st.sidebar.header("Display Options")
show_gc = st.sidebar.checkbox("Show GC Content", value=True)
show_translation = st.sidebar.checkbox("Show Translated Protein", value=True)
show_base_plot = st.sidebar.checkbox("Show Base Composition Plot", value=True)

# Upload section
uploaded_file = st.file_uploader("📁 Upload a DNA FASTA file", type=["fasta", "fa"])

def analyze_sequence(seq_record):
    seq = seq_record.seq
    gc = gc_fraction(seq) * 100

    # Trim the sequence if it's not a multiple of 3
    trimmed_length = len(seq) - (len(seq) % 3)
    trimmed_seq = seq[:trimmed_length]
    protein = trimmed_seq.translate(to_stop=True)

    base_counts = pd.Series(list(seq)).value_counts().sort_index()
    base_df = base_counts.rename_axis("Base").reset_index(name="Count")
    return gc, protein, base_df

# If file is uploaded
if uploaded_file:
    stringio = io.StringIO(uploaded_file.getvalue().decode("utf-8"))
    records = list(SeqIO.parse(stringio, "fasta"))

    if not records:
        st.error("❌ No sequences found in file. Make sure it's a proper FASTA format.")
    else:
        st.success(f"✅ {len(records)} sequence(s) loaded!")

        for record in records:
            st.subheader(f"🧬 Sequence ID: {record.id}")
            st.code(str(record.seq), language="text")

            gc_content, protein, base_df = analyze_sequence(record)

            if show_gc:
                st.write(f"**GC Content:** `{gc_content:.2f}%`")

            if show_translation:
                st.write("**Translated Protein (First ORF):**")
                st.code(str(protein), language="text")

            if show_base_plot:
                fig = px.bar(base_df, x="Base", y="Count", title="Base Composition", color="Base")
                st.plotly_chart(fig, use_container_width=True)
else:
    st.warning("👆 Please upload a FASTA file to begin analysis.")