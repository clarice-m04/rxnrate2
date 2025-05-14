import streamlit as st
import sys
import os
import pubchempy as pcp

from rxnrate2.ODE_linearrxn import solve_reaction, plot_solution
from rdkit import Chem
from rdkit.Chem import Draw
from PIL import Image, ImageDraw, ImageFont



### Definition of functions to draw the reactions using SMILES ###

## Helper: fallback SMILESdef get_smiles(query):
def get_smiles(query): 
    fallback_smiles = {
        'H2O': 'O',
        'CO2': 'O=C=O',
        'O2': 'O=O',
        'H2': '[H][H]',
        'N2': 'N#N',
        'CH4': 'C',
        'NH3': 'N',#
    }
    try:#
        compound = pcp.get_compounds(query, 'name')
        if compound:
            return compound[0].isomeric_smiles
    except:
        pass
    return fallback_smiles.get(query.strip(), None)

## Helper: drawing function
def draw_reaction(reagent_smiles, product_smiles, reagent_label, product_label, conc_reagent, conc_product, kf, kb):
    # Try to use Times New Roman; fallback to default
    try:
        font = ImageFont.truetype("Times New Roman.ttf", 14)
    except:
        font = ImageFont.load_default()


    mol1 = Chem.MolFromSmiles(reagent_smiles) if reagent_smiles else None
    mol2 = Chem.MolFromSmiles(product_smiles) if product_smiles else None

    if not mol1 or not mol2:
        st.warning("Could not generate one or both molecule images.")
        return None

    img1 = Draw.MolToImage(mol1, size=(200, 200))
    img2 = Draw.MolToImage(mol2, size=(200, 200))

    canvas = Image.new('RGB', (500, 250), 'white')
    draw = ImageDraw.Draw(canvas)

    # Paste molecule images
    canvas.paste(img1, (10, 10))
    canvas.paste(img2, (290, 10))

    # Draw double arrows
    arrow_y = 100
    draw.line((220, arrow_y - 10, 280, arrow_y - 10), fill='black', width=2)
    draw.line((280, arrow_y + 10, 220, arrow_y + 10), fill='black', width=2)
    draw.polygon([(275, arrow_y - 13), (285, arrow_y - 10), (275, arrow_y - 7)], fill='black')
    draw.polygon([(225, arrow_y + 7), (215, arrow_y + 10), (225, arrow_y + 13)], fill='black')

    # Decimal format for kinetic constants
    draw.text((230, arrow_y - 30), f"kf = {kf:.2f}", fill='black', font=font)
    draw.text((230, arrow_y + 15), f"kb = {kb:.2f}", fill='black', font=font)

    # Labels and concentrations
    draw.text((10, 210), f"{reagent_label}", fill='black', font=font)
    draw.text((10, 225), f"[{reagent_label}]i = {conc_reagent:.2f}M", fill='black', font=font)

    draw.text((290, 210), f"{product_label}", fill='black', font=font)
    draw.text((290, 225), f"[{product_label}]i = {conc_product:.2f}M", fill='black', font=font)

    return canvas



### Define interactive page with buttons ###

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
st.title("Welcome in linear reaction part")

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

if st.session_state.get("last_action_message"):
    st.success(st.session_state["last_action_message"])
    if os.path.exists(st.session_state["new_image_to_show"]):
        st.image(
            st.session_state["new_image_to_show"],
            caption="Reaction Rate Picture",
            use_container_width=True
        )
    # Clear the message and image info so it's not shown again
    del st.session_state["last_action_message"]
    del st.session_state["new_image_to_show"]

## Form input
col1, col2 = st.columns(2)
with col1:
    st.subheader("Reagents")

    reagent = st.text_input("Reagent (Formula or Name, e.g. H2O)")
    init_conc_reagent = st.number_input("Initial concentration of Reagent", min_value=0.0, value=1.0, format="%.3f")

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
  

   