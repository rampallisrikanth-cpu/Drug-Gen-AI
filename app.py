import streamlit as st
import pandas as pd
from model import predict_effectiveness

# App Title
st.title("🧬 Genomic Drug Effectiveness Predictor")
st.subheader("Personalized Medicine Prototype")

# Sidebar
st.sidebar.title("📁 Upload Section")
uploaded_file = st.sidebar.file_uploader("Upload Genome CSV", type=["csv"])

st.markdown("""
<style>
    .main {
        background-color: #F0F2F6;
    }
</style>
""", unsafe_allow_html=True)

if uploaded_file:
    st.success("File uploaded successfully!")

    df = pd.read_csv(uploaded_file)
    st.write("### 📊 Uploaded Genomic Data")
    st.dataframe(df)

    st.write("### 🔍 Prediction")
    genes = df.iloc[0].to_dict()
    score = predict_effectiveness(genes)

    st.markdown(f"""
    <div style='padding:15px;background:white;border-radius:12px;box-shadow:0px 3px 10px rgba(0,0,0,0.1);'>
    <h2>Drug Effectiveness: {score}%</h2>
    </div>
    """, unsafe_allow_html=True)
else:
    st.info("Please upload a genomic CSV using the sidebar.")

