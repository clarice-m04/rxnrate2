import streamlit as st
import base64
from pathlib import Path
from rxnrate2.Interface_rxnrate.functions import set_background


#from rxnrate2.Interface_rxnrate.pages import SimpleRxn, ComplexRxn

# Times new roman font
st.markdown("""
    <style>
    * {
        font-family: 'Times New Roman', Times, serif !important;
    }
    html, body, [class*="css"] {
        font-family: 'Times New Roman', Times, serif !important;
    }
    .stText, .stMarkdown, .stDataFrame, .stTable, .stButton, .stHeader, .stSubheader, .stCaption {
        font-family: 'Times New Roman', Times, serif !important;
    }s
    </style>
""", unsafe_allow_html=True)

# backround of the application

set_background("rxnrate.jpg")

# Initialisation of the active page
if "page" not in st.session_state:
    st.session_state.page = "home"

def home_page():
    st.header("Welcome to RxnRate!📈✨")
    st.subheader("Are you curious to see how your reaction advances ?")
    st.write("This program provides you the progression of your reaction over time")
    st.write("Now let's see what kind of reaction you have! " \
    "To do this, read the two categories below carefully and select the option that best describes your system. " \
    "Click on the button and you'll be taken straight to the perfect page for your model!")
    


if st.session_state.page == "home":
   home_page()

st.divider()

col1, col2 = st.columns(2)
with col1:
    st.subheader("**Main characteristics:**")
    st.write("- Reaction of type A -> B -> C")
    st.write("- One reactant gives one product")
    st.page_link("pages/SimpleRxn.py", label="Simple reaction", icon="1️⃣")

with col2:
    st.subheader("**Main characteristics:**")
    st.write("- Reaction of type A + B -> C")
    st.write("- Two or more reactants give one product")
    st.page_link("pages/ComplexRxn2.py", label="More complex reaction", icon="2️⃣")
    
