import streamlit as st
import pandas as pd
from model import predict_effectiveness

st.title("Genomic Drug Effectiveness AI Prototype")
st.write("Upload your genomic data CSV to predict how effective the drug will be.")

uploaded_file = st.file_uploader("Upload Genome CSV", type=["csv"])

if uploaded_file:
    df = pd.read_csv(uploaded_file)
    st.write("Uploaded Genome Data:")
    st.dataframe(df)

    # Convert first row to dictionary for prototype
    genes = df.iloc[0].to_dict()

    score = predict_effectiveness(genes)

    st.subheader("Predicted Drug Effectiveness:")
    st.success(f"{score} % Effective")
