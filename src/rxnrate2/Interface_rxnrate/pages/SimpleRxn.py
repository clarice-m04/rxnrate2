import streamlit as st

import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../../')))

from rxnrate2.ODE_linearrxn import solve_reaction, plot_solution

from PIL import Image

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

#def linear_rxn():
st.title("Welcome in linear reaction part")

###Inputs from the user###
import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw
from PIL import Image, ImageDraw, ImageFont
import pubchempy as pcp

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

# Helper: fallback SMILES
def get_smiles(query):
    fallback_smiles = {
        'H2O': 'O',
        'CO2': 'O=C=O',
        'O2': 'O=O',
        'H2': '[H][H]',
        'N2': 'N#N',
        'CH4': 'C',
        'NH3': 'N',
    }
    try:
        compound = pcp.get_compounds(query, 'name')
        if compound:
            return compound[0].isomeric_smiles
    except:
        pass
    return fallback_smiles.get(query.strip(), None)

# Helper: drawing function
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


# Form input
col1, col2 = st.columns(2)
with col1:
    # If previous reactions exist, use the last product as the default reagent
    if st.session_state.fixed_reagent:
        reagent = st.session_state.fixed_reagent
        st.markdown(f"**Reagent:** {reagent} *(auto-filled from previous product)*")
    else:
        reagent = st.text_input("Reagent (Formula or Name, e.g. H2O)")

    k_forward = st.number_input("k_forward", min_value=0.0, value=1.0, format="%.6f")
    init_conc_reagent = st.number_input("Initial concentration of Reagent", min_value=0.0, value=1.0, format="%.3f")
with col2:
    product = st.text_input("Product (Formula or Name, e.g. CO2)")
    k_backward = st.number_input("k_backward", min_value=0.0, value=0.5, format="%.6f")
    init_conc_product = st.number_input("Initial concentration of Product", min_value=0.0, value=0.0, format="%.3f")

# Submit button
if st.button("Add Reaction"):
    if reagent and product:
        reagent_smiles = get_smiles(reagent)
        product_smiles = get_smiles(product)

        if reagent_smiles and product_smiles:
            st.session_state.reagents.append(reagent)
            st.session_state.products.append(product)
            st.session_state.kf.append(k_forward)
            st.session_state.kb.append(k_backward)

            for specie, conc in [(reagent, init_conc_reagent), (product, init_conc_product)]:
                st.session_state.init_conc[specie] = conc

            st.session_state.reactions.append({
                'reagent': reagent,
                'product': product,
                'reagent_smiles': reagent_smiles,
                'product_smiles': product_smiles,
                'kf': k_forward,
                'kb': k_backward
            })

            # ✅ Append tuple with kb=None if it's 0    

            kb_value = k_backward if k_backward != 0 else None
            reaction_tuple = (reagent, product, k_forward, kb_value)
            st.session_state.reaction_tuples.append(reaction_tuple)

            # After adding the reaction
            if not st.session_state.fixed_reagent:
            # Set reagent to the first manually entered one
                st.session_state.fixed_reagent = product
            else:
            # Update reagent to next product to continue chain
                st.session_state.fixed_reagent = product

            # Append reagent if not already present
            if reagent not in st.session_state.species_list:
                st.session_state.species_list.append(reagent)

            # Append product if not already present
            if product not in st.session_state.species_list:
                st.session_state.species_list.append(product)

            st.success(f"Added reaction: {reagent} ⇌ {product}")
        else:
            st.error("Could not resolve SMILES for reagent or product.")
    if st.session_state.reactions:
        i_conc_list = list(st.session_state.init_conc.values())

        s,m = solve_reaction(st.session_state.species_list, st.session_state.reaction_tuples, i_conc_list)
        plot_solution(s,st.session_state.species_list)
        current_dir = os.path.dirname(__file__)        # points to pages/
        parent_dir = os.path.dirname(current_dir)       # points to Interface_rxnrate/
        image_path = os.path.join(parent_dir, "reaction_plot.jpg")
        st.image(image_path, caption="Reaction Rate Picture", use_container_width=True)

    
    else:
        st.error("Please enter both reagent and product.")


if st.button("Remove Last Reaction"):
    if st.session_state.reactions:
        # Remove the last elements
        st.session_state.reactions.pop()
        st.session_state.reagents.pop()
        st.session_state.products.pop()
        st.session_state.kf.pop()
        st.session_state.kb.pop()
        st.session_state.reaction_tuples.pop()

        # Remove species if no longer used
        used_species = set()
        for r in st.session_state.reactions:
            used_species.add(r['reagent'])
            used_species.add(r['product'])

        st.session_state.init_conc = {
            k: v for k, v in st.session_state.init_conc.items() if k in used_species
        }

        # Rebuild species_list preserving order
        species_seen = set()
        st.session_state.species_list = []
        for r in st.session_state.reactions:
            if r["reagent"] not in species_seen:
                st.session_state.species_list.append(r["reagent"])
                species_seen.add(r["reagent"])
            if r["product"] not in species_seen:
                st.session_state.species_list.append(r["product"])
                species_seen.add(r["product"])

        # Update fixed reagent to new last product (or None if empty)
        if st.session_state.reactions:
            st.session_state.fixed_reagent = st.session_state.reactions[-1]["product"]
        else:
            st.session_state.fixed_reagent = None

        st.success("✅ Last reaction removed.")
        st.rerun()
    else:
        st.warning("⚠️ No reactions to remove.")

if st.button("Clear All Reactions"):
    st.session_state.reactions = []
    st.session_state.reagents = []
    st.session_state.products = []
    st.session_state.kf = []
    st.session_state.kb = []
    st.session_state.init_conc = {}
    st.session_state.fixed_reagent = None
    st.session_state.species_list = []
    st.session_state.reaction_tuples = []

    st.success("All reactions have been cleared.")

# Visualization
st.header("Reaction Visualizations")

for idx, rxn in enumerate(st.session_state.reactions):
    st.markdown(f"### Reaction {idx+1}: {rxn['reagent']} ⇌ {rxn['product']}")
    
    image = draw_reaction(
        rxn['reagent_smiles'],
        rxn['product_smiles'],
        rxn['reagent'],
        rxn['product'],
        st.session_state.init_conc[rxn['reagent']],
        st.session_state.init_conc[rxn['product']],
        rxn['kf'],
        rxn['kb']
    )

    if image:
        st.image(image)

st.subheader("Stored Reaction Tuples")
for r in st.session_state.reaction_tuples:
    st.write(r)

st.markdown("### All Species (Ordered, No Duplicates):")
st.write(st.session_state.species_list)
