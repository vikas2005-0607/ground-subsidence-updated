
import streamlit as st

st.set_page_config(
    page_title="Ground Subsidence Prediction",
    page_icon="🌍",
    layout="wide"
)

st.title("Ground Subsidence Prediction System")
st.caption("Mining Engineering | IIT (ISM) Dhanbad")

st.info(
    "This application will implement your MATLAB "
    "subsidence prediction model."
)

st.header("Observation Point")

col1, col2 = st.columns(2)

with col1:
    x = st.number_input("X coordinate (m)", value=0.0)

with col2:
    y = st.number_input("Y coordinate (m)", value=0.0)

if st.button("Run Prediction"):
    st.write(f"Observation point: X = {x}, Y = {y}")
    st.warning(
        "The MATLAB calculation engine has not been "
        "integrated yet. No subsidence prediction was calculated."
    )
