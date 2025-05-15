import streamlit as st
import sys
import os
import pubchempy as pcp

from rxnrate2.ODE_linearrxn import solve_reaction, plot_solution
from rdkit import Chem
from rdkit.Chem import Draw
from PIL import Image, ImageDraw, ImageFont
from rxnrate2.Interface_rxnrate.pages.SimpleRxn import get_smiles, draw_reaction


# Define interactive page with buttons, in the same way as the one made for linear

## Times new roman font
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

## Title
st.title("Welcome in nonlinear reaction part")

## Subtitle
st.title("Chemical Reaction Collector")

## Data structures
if 'reagents' not in st.session_state:
    st.session_state.reagents = []
    st.session_state.products = []
    st.session_state.kf = []
    st.session_state.kb = []
    st.session_state.init_conc = {}
    st.session_state.reactions = []

if "species_list" not in st.session_state:
    st.session_state.species_list = []

if "reaction_tuples" not in st.session_state:
    st.session_state.reaction_tuples = []

#if st.session_state.get("last_action_message"):
 #   st.success(st.session_state["last_action_message"])
  #  if os.path.exists(st.session_state["new_image_to_show"]):
   #     st.image(
    #        st.session_state["new_image_to_show"],
     #       caption="Reaction Rate Picture",
      #      use_container_width=True)
    # Clear the message and image info so it's not shown again
    #del st.session_state["last_action_message"]
    #del st.session_state["new_image_to_show"]

## Form input
col1, col2 = st.columns(2)
with col1:
    st.subheader("Reagents")

    reagent = st.text_input("Reagent (Formula or Name, e.g. H2O)")
    init_conc_reagent = st.number_input("Initial concentration of Reagent", min_value=0.0, value=1.0, format="%.3f")


    #for i in range(len(st.session_state.reactants)):
     #   st.session_state.reactants[i] = st.text_input(f"Reactant {i+1}", key=f"reactant_{i}")
    #kf = st.number_input("Forward Rate Constant (kf)", min_value=0.0, value=1.0, format="%.5f")

    #if st.button("Add Reactant Field"):
     #   st.session_state.reactants.append([])


    # Button to add the reagent and the associated concentration in two different lists; 
        #   one for the reagents and one for the concentration
    if st.button("Add this reagent"):
        if reagent:
            reagents = st.session_state.reagents.append(reagent)
            st.success(f"'{reagent}' is added to your reagents list !")
            st.write(f"Your list is then : {reagents}")
        else:
            st.write("Please insert a reagent")

        if init_conc_reagent:
            for specie, conc in [(reagent, init_conc_reagent)]:
                st.session_state.init_conc[specie] = conc

            st.success(f"'{init_conc_reagent}' is added to your concentration of reagents list !")
        else:
            st.write("Please insert a concentration for this reagent")

    k_forward = st.number_input("k_forward", min_value=0.0, value=1.0, format="%.6f")
    


with col2:

    st.subheader("Products")

    product = st.text_input("Product (Formula or Name, e.g. CH4)")
    init_conc_product = st.number_input("Initial concentration of Product", min_value=0.0, value=1.0, format="%.3f")

    # Button to add the product and the associated concentration in two different lists; 
        #   one for the product and one for the concentration
    if st.button("Add this product"):
        if product:
            products = st.session_state.products.append(product)
            st.success(f"'{product}' is added to your products list !")
            st.write(f"Your list is then : {products}")
        else:
            st.write("Please insert a product")
        
        if init_conc_product:
            for specie, conc in  [(product, init_conc_product)]:
                st.session_state.init_conc[specie] = conc

            st.success(f"'{init_conc_product}' is added to your concentration of products list !")
        else:
            st.write("Please insert a concentration for this product")

    k_backward = st.number_input("k_backward", min_value=0.0, value=1.0, format="%.6f")
  

   