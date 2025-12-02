# app.py
import streamlit as st
import pandas as pd
from utils import parse_23andme_csv, parse_vcf
from model import compute_drug_scores, RSID_TO_GENE, DRUGS

st.set_page_config(page_title="Genomic Drug Effectiveness AI", layout="centered")

st.title("🧬 Genomic Drug Effectiveness AI (Prototype)")
st.markdown(
    "Upload a **23andMe-style CSV** or a **VCF** file. The app reads pharmacogenomic markers "
    "and returns a simple predicted effectiveness score for a few common drugs. **Non-clinical demo only.**"
)

uploaded = st.file_uploader("Upload genomic file (CSV or VCF)", type=["csv", "txt", "vcf", "gz"])
if uploaded is None:
    st.info("Try the sample files in `sample_data/` or upload your own 23andMe-style CSV / VCF.")
    st.markdown("**What to upload:** a 23andMe raw data file (rsid,chromosome,position,genotype) or a VCF containing rsIDs.")
else:
    st.write("File uploaded:", uploaded.name)
    # detect type
    fname = uploaded.name.lower()
    try:
        if fname.endswith(".vcf") or fname.endswith(".vcf.gz") or uploaded.name.lower().endswith(".gz"):
            # assume VCF (gz handler: Streamlit gives file-like; parse raw text)
            snp_dict = parse_vcf(uploaded)
            st.success(f"Parsed VCF: {len(snp_dict)} rsIDs found (showing up to 10).")
        else:
            # CSV or txt
            snp_dict = parse_23andme_csv(uploaded)
            st.success(f"Parsed CSV: {len(snp_dict)} markers found (showing up to 10).")
    except Exception as e:
        st.error("Failed to parse file: " + str(e))
        st.stop()

    # show top markers
    if len(snp_dict) == 0:
        st.warning("No recognizable markers found in the file. Make sure the file contains rsIDs or a genotype table.")
    else:
        df_preview = pd.DataFrame(list(snp_dict.items()), columns=["marker", "genotype"]).head(25)
        st.write("Sample of parsed markers:")
        st.dataframe(df_preview)

    # compute drug scores
    with st.spinner("Analyzing pharmacogenomic markers and computing drug-effectiveness..."):
        results = compute_drug_scores(snp_dict)

    st.markdown("## Results")
    for drug, res in results.items():
        score = res["score"]
        phenotype = res["phenotype"]
        explanation = res["explanation"]
        col1, col2 = st.columns([1,3])
        with col1:
            # simple visual score
            st.metric(label=drug, value=f"{score} %")
        with col2:
            st.write(f"**Phenotype:** {phenotype.capitalize()}")
            st.write(explanation)
            # advice-like message (prototype)
            if drug == "Clopidogrel" and phenotype == "poor":
                st.warning("⚠️ Warning: Poor CYP2C19 function may reduce Clopidogrel activation. This prototype is not clinical advice; consult a clinician.")
            if drug == "Codeine" and phenotype == "poor":
                st.info("Note: Codeine may be ineffective for poor CYP2D6 metabolizers.")
    st.markdown("---")
    st.caption("This is a simplified prototype. Real clinical decision requires validated lab reports and medical professionals.")
