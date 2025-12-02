import streamlit as st
import pandas as pd
import re

st.title("Genomic Drug Response AI (Prototype)")

st.write("Upload your genomic data (CSV or VCF) to predict drug effectiveness.")

# -----------------------------
# FILE PARSER FUNCTIONS INSIDE SAME FILE
# -----------------------------

def parse_csv(file):
    df = pd.read_csv(file)
    return df

def parse_vcf(file):
    rows = []
    for line in file:
        line = line.decode("utf-8")
        if line.startswith("#"):
            continue
        parts = line.strip().split("\t")
        if len(parts) > 4:
            chrom, pos, rsid, ref, alt = parts[:5]
            rows.append([rsid, chrom, pos, ref + "/" + alt])
    df = pd.DataFrame(rows, columns=["rsid", "chromosome", "position", "genotype"])
    return df

# -----------------------------
# DRUG–GENE DATABASE (VERY BASIC PROTOTYPE)
# -----------------------------
drug_genes = {
    "Paracetamol": {
        "UGT1A1": {"rs8175347": {"TA6/TA7": "Slow", "TA7/TA7": "Very Slow"}}
    },
    "Ibuprofen": {
        "CYP2C9": {"rs1057910": {"AA": "Fast", "AC": "Normal", "CC": "Slow"}}
    },
    "Caffeine": {
        "CYP1A2": {"rs762551": {"AA": "Fast", "AC": "Medium", "CC": "Slow"}}
    }
}

# -----------------------------
# EFFECTIVENESS SCORE
# -----------------------------
def score_effectiveness(speed):
    if speed == "Fast":
        return 85
    if speed == "Medium" or speed == "Normal":
        return 70
    if speed == "Slow":
        return 40
    if speed == "Very Slow":
        return 20
    return 50

# -----------------------------
# ANALYSIS FUNCTION
# -----------------------------
def analyze(df, drug):
    result = {}
    gene_info = drug_genes[drug]

    for gene, snps in gene_info.items():
        for rsid, variants in snps.items():
            row = df[df["rsid"] == rsid]
            if not row.empty:
                gt = row.iloc[0]["genotype"]
                # Genotype formatting correction
                gt = gt.replace("/", "").replace("-", "")
                for v, metabolism in variants.items():
                    v_clean = v.replace("/", "")
                    if gt == v_clean:
                        result[gene] = metabolism

    if result:
        metabolism = list(result.values())[0]
        effectiveness = score_effectiveness(metabolism)
        return metabolism, effectiveness
    else:
        return "Unknown", 50

# -----------------------------
# FILE UPLOAD
# -----------------------------
file = st.file_uploader("Upload CSV or VCF file", type=["csv", "txt", "vcf"])

if file:
    if file.name.endswith(".csv"):
        df = parse_csv(file)
    else:
        df = parse_vcf(file)

    st.success("File uploaded & processed!")
    st.dataframe(df.head())

    drug = st.selectbox("Select a drug to evaluate:", list(drug_genes.keys()))

    if st.button("Predict Drug Effectiveness"):
        metabolism, effectiveness = analyze(df, drug)

        st.subheader("Result")
        st.write(f"**Metabolism speed:** {metabolism}")
        st.write(f"**Predicted drug effectiveness:** {effectiveness} / 100")
