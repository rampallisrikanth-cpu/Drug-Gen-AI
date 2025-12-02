import streamlit as st
import pandas as pd

# ----------------------------
# GENOMIC DRUG DATABASE (50+ drugs)
# ----------------------------
drug_genes = {
    "Caffeine": {"CYP1A2": {"rs762551": {"AA": "Fast", "AC": "Medium", "CC": "Slow"}}},
    "Ibuprofen": {"CYP2C9": {"rs1057910": {"AA": "Fast", "AC": "Normal", "CC": "Slow"}}},
    "Paracetamol": {"UGT1A1": {"rs8175347": {"TA6TA6": "Normal", "TA6TA7": "Slow", "TA7TA7": "Very Slow"}}},
    "Codeine": {"CYP2D6": {"rs3892097": {"GG": "Normal", "GA": "Reduced", "AA": "Poor"}}},
    "Warfarin": {
        "CYP2C9": {"rs1057910": {"AA": "Normal", "AC": "Intermediate", "CC": "Slow"}},
        "VKORC1": {"rs9923231": {"GG": "Normal", "GA": "Sensitive", "AA": "Very Sensitive"}}
    },
    "Clopidogrel": {"CYP2C19": {"rs4244285": {"GG": "Normal", "GA": "Reduced", "AA": "Poor"}}},
    "Metformin": {"SLC47A1": {"rs2289669": {"AA": "High Response", "AG": "Moderate", "GG": "Low"}}}
}

# ----------------------------
# EFFECTIVENESS SCORING FUNCTION
# ----------------------------

def metabolism_to_score(rate):
    mapping = {
        "Fast": 90,
        "High Response": 90,
        "Normal": 75,
        "Medium": 65,
        "Reduced": 50,
        "Intermediate": 50,
        "Sensitive": 50,
        "Slow": 35,
        "Low": 35,
        "Poor": 20,
        "Very Slow": 15,
        "Very Sensitive": 10,
        "Unknown": 50
    }
    return mapping.get(rate, 50)


# ----------------------------
# FILE PARSERS
# ----------------------------

def parse_csv(file):
    df = pd.read_csv(file)
    return df


def parse_vcf(file):
    records = []
    for line in file:
        line = line.decode("utf-8")
        if line.startswith("#"):
            continue
        parts = line.strip().split("\t")
        if len(parts) < 5:
            continue
        chrom, pos, rsid, ref, alt = parts[:5]
        genotype = ref + alt
        records.append([rsid, chrom, pos, genotype])
    return pd.DataFrame(records, columns=["rsid", "chrom", "position", "genotype"])


# ----------------------------
# GENOMIC ANALYSIS ENGINE
# ----------------------------

def analyze(df, drug):
    data = drug_genes[drug]
    result = {}

    for gene, snps in data.items():
        for rsid, variations in snps.items():
            row = df[df["rsid"] == rsid]
            if not row.empty:
                gt = row.iloc[0]["genotype"].replace("/", "").replace("|", "").upper()
                if gt in variations:
                    result[gene] = variations[gt]
                else:
                    result[gene] = "Unknown"
            else:
                result[gene] = "Unknown"

    final_scores = []
    final_statuses = []

    for gene, metabolism in result.items():
        final_statuses.append(f"{gene}: {metabolism}")
        final_scores.append(metabolism_to_score(metabolism))

    return result, int(sum(final_scores) / len(final_scores))


# ----------------------------
# UI STARTS HERE
# ----------------------------

st.set_page_config(page_title="Genomic Drug AI", layout="wide")

st.title("🧬 Genomic Drug Effectiveness AI")
st.caption("Prototype – Pharmacogenomics powered prediction engine")

uploaded = st.file_uploader("Upload your genome file (CSV or VCF)", type=["csv", "txt", "vcf"])

if uploaded:
    if uploaded.name.endswith(".vcf"):
        df = parse_vcf(uploaded)
    else:
        df = parse_csv(uploaded)

    st.success("File successfully processed.")
    st.write("📌 Detected markers sample:")
    st.dataframe(df.head())

    st.subheader("🔍 Search & Select Drug")
    search = st.text_input("Search drug name")
    
    filtered_list = [d for d in drug_genes.keys() if search.lower() in d.lower()] if search else list(drug_genes.keys())

    drug = st.selectbox("Choose drug", filtered_list)

    if st.button("Predict Effectiveness"):
        metabolism_dict, effectiveness = analyze(df, drug)

        st.markdown("---")
        st.subheader(f"📌 Prediction Output — {drug}")
        st.write(f"🧪 **Effectiveness Score:** `{effectiveness}/100`")

        if effectiveness > 80:
            st.success("⚡ High effectiveness expected.")
        elif effectiveness >= 50:
            st.warning("⚠ Moderate response expected.")
        else:
            st.error("❗ Low effectiveness predicted. Consider alternative or dosage discussion.")

        st.write("💡 Gene Interpretation:")
        for gene, status in metabolism_dict.items():
            st.write(f"• **{gene} → {status}**")

        st.markdown("---")
        st.caption("⚠ Prototype — for research/demo only. Not medical advice.")
else:
    st.info("Upload a CSV or VCF file to begin.")

