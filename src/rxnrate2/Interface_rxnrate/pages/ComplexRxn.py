import streamlit as st

from rxnrate2.ODE_nonlinear import solve_reactions, plot_solution
from rdkit import Chem
from PIL import Image, ImageDraw, ImageFont
import pubchempy as pcp

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

# Title of the page; definition of the nonlinear part:
st.title("Welcome in non-linear reaction part")


