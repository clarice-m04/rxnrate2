import streamlit as st
from rxnrate2.ODE_nonlinear import solve_reactions, plot_solution
from rdkit import Chem
from PIL import Image, ImageDraw, ImageFont
import pubchempy as pcp
import base64
from pathlib import Path
from rxnrate2.Interface_rxnrate.functions import set_background
from SimpleRxn import get_smiles, draw_reaction


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
    }
    </style>
""", unsafe_allow_html=True)

# Set of the backround
set_background("./rxnratepic.jpg")

# Title of the page; definition of the nonlinear part:
st.title("Welcome in non-linear reaction part")


# same structure of the chemical reaction collector
st.title("Chemical Reaction Collector")

# Data structures
if 'reagents' not in st.session_state:
    st.session_state.reagents = []
    st.session_state.products = []
    st.session_state.kf = []
    st.session_state.kb = []
    st.session_state.init_conc = {}
    st.session_state.reactions = []

if "fixed_reagent" not in st.session_state:
    st.session_state.fixed_reagent = None

if "species_list" not in st.session_state:
    st.session_state.species_list = []

if "reaction_tuples" not in st.session_state:
    st.session_state.reaction_tuples = []

# Imput of the users 

col1, col2 = st.columns(2)

with col1:
    st.subheader("Reagents")
    reagent = st.text_input("Reagent (Formula or Name, e.g. H2O)")
   
    init_conc_reagent = st.number_input("Initial concentration of Reagent", min_value=0.0, value=1.0, format="%.3f")

    # Button to add the reagent and the associated concentration in two different lists; 
        #   one for the reagents and one for the concentration
    if st.button("Add this reagent "):
        if reagent:
            st.session_state.reagents.append(reagent)
            st.success(f"'{reagent}' is added to your reagents list !")
        else:
            st.write("Please insert a reagent")
        
        if init_conc_reagent:
            st.session_state.init_conc.append(init_conc_reagent)
            st.success(f"'{init_conc_reagent}' is added to your concentration of reagents list !")

    st.write("Please insert the velocity constant of your reaction")
    k_forward = st.number_input("k_forward", min_value=0.0, value=1.0, format="%.6f")
