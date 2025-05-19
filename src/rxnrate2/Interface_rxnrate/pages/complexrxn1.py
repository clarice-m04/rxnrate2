import streamlit as st
import os
import pubchempy as pcp
from rxnrate2.ODE_linearrxn import solve_reaction, plot_solution
from rdkit import Chem
from rdkit.Chem import Draw
from PIL import Image, ImageDraw, ImageFont

# Custom Times New Roman font styling
st.markdown("""
    <style>
    * {
        font-family: 'Times New Roman', Times, serif !important;
    }
    html, body, [class*="css"] {
        font-family: 'Times New Roman', Times, serif !important;
    }
    </style>
""", unsafe_allow_html=True)

st.title("Linear Chemical Reaction Simulator")
st.subheader("Define, Visualize, and Simulate Reversible Reactions")

# Initialize session state
for key in ["reagents", "products", "kf", "kb", "init_conc", "reactions", "species_list", "reaction_tuples"]:
    if key not in st.session_state:
        st.session_state[key] = [] if key != "init_conc" else {}

# Helper to get SMILES from compound name
@st.cache_data
def get_smiles(query):
    fallback_smiles = {
        'H2O': 'O', 'CO2': 'O=C=O', 'O2': 'O=O', 'H2': '[H][H]',
        'N2': 'N#N', 'CH4': 'C', 'NH3': 'N'
    }
    try:
        compound = pcp.get_compounds(query, 'name')
        if compound:
            return compound[0].isomeric_smiles
    except:
        pass
    return fallback_smiles.get(query.strip(), None)

# Helper to draw reaction scheme
def draw_reaction(reagent_smiles, product_smiles, reagent_label, product_label, conc_reagent, conc_product, kf, kb):
    try:
        font = ImageFont.truetype("Times New Roman.ttf", 14)
    except:
        font = ImageFont.load_default()

    mol1 = Chem.MolFromSmiles(reagent_smiles) if reagent_smiles else None
    mol2 = Chem.MolFromSmiles(product_smiles) if product_smiles else None

    if not mol1 or not mol2:
        return None

    img1 = Draw.MolToImage(mol1, size=(200, 200))
    img2 = Draw.MolToImage(mol2, size=(200, 200))

    canvas = Image.new('RGB', (500, 250), 'white')
    draw = ImageDraw.Draw(canvas)
    canvas.paste(img1, (10, 10))
    canvas.paste(img2, (290, 10))

    # Draw reversible arrow and rate constants
    arrow_y = 100
    draw.line((220, arrow_y - 10, 280, arrow_y - 10), fill='black', width=2)
    draw.line((280, arrow_y + 10, 220, arrow_y + 10), fill='black', width=2)
    draw.polygon([(275, arrow_y - 13), (285, arrow_y - 10), (275, arrow_y - 7)], fill='black')
    draw.polygon([(225, arrow_y + 7), (215, arrow_y + 10), (225, arrow_y + 13)], fill='black')

    draw.text((230, arrow_y - 30), f"kf = {kf:.2f}", fill='black', font=font)
    draw.text((230, arrow_y + 15), f"kb = {kb:.2f}", fill='black', font=font)
    draw.text((10, 210), f"{reagent_label}", fill='black', font=font)
    draw.text((10, 225), f"[{reagent_label}]i = {conc_reagent:.2f}M", fill='black', font=font)
    draw.text((290, 210), f"{product_label}", fill='black', font=font)
    draw.text((290, 225), f"[{product_label}]i = {conc_product:.2f}M", fill='black', font=font)

    return canvas

# User inputs for reaction definition
col1, col2 = st.columns(2)
with col1:
    reagent = st.text_input("Reagent (e.g. H2O)")
    kf = st.number_input("Forward rate constant (kf)", min_value=0.0, value=1.0)
    conc_r = st.number_input("Initial concentration of reagent", min_value=0.0, value=1.0)
with col2:
    product = st.text_input("Product (e.g. CO2)")
    kb = st.number_input("Backward rate constant (kb)", min_value=0.0, value=0.5)
    conc_p = st.number_input("Initial concentration of product", min_value=0.0, value=0.0)

# Add reaction
if st.button("Add Reaction"):
    if reagent and product:
        sm_r = get_smiles(reagent)
        sm_p = get_smiles(product)
        if sm_r and sm_p:
            st.session_state.reagents.append(reagent)
            st.session_state.products.append(product)
            st.session_state.kf.append(kf)
            st.session_state.kb.append(kb)
            st.session_state.init_conc[reagent] = conc_r
            st.session_state.init_conc[product] = conc_p
            reaction_tuple = (reagent, product, kf, None if kb == 0 else kb)
            st.session_state.reaction_tuples.append(reaction_tuple)
            for s in [reagent, product]:
                if s not in st.session_state.species_list:
                    st.session_state.species_list.append(s)
            st.success(f"Reaction added: {reagent} ⇌ {product}")
        else:
            st.error("Could not resolve SMILES for one or more inputs.")

# Solve and plot reactions
if st.button("Simulate Reactions"):
    if st.session_state.reaction_tuples:
        try:
            os.makedirs("figures", exist_ok=True)
            filename = f"figures/{','.join(st.session_state.reagents)}_to_{','.join(st.session_state.products)}.jpg"
            conc_list = [st.session_state.init_conc[s] for s in st.session_state.species_list]
            sol, _ = solve_reaction(st.session_state.species_list, st.session_state.reaction_tuples, conc_list)
            plot_solution(sol, st.session_state.species_list, filename=filename)
            st.image(filename, caption="Reaction Progress", use_container_width=True)
        except Exception as e:
            st.error(f"Simulation failed: {e}")
    else:
        st.warning("No reactions defined.")

# Display all reaction visuals
st.header("Defined Reactions")
for idx, (r, p, kf, kb) in enumerate(st.session_state.reaction_tuples):
    sm_r = get_smiles(r)
    sm_p = get_smiles(p)
    image = draw_reaction(sm_r, sm_p, r, p, st.session_state.init_conc[r], st.session_state.init_conc[p], kf, kb or 0)
    if image:
        st.image(image, caption=f"Reaction {idx+1}: {r} ⇌ {p}")

# Show all species and reaction tuples
st.subheader("All Species")
st.write(st.session_state.species_list)

st.subheader("Stored Reaction Tuples")
for rt in st.session_state.reaction_tuples:
    st.write(rt)

# Option to clear all
if st.button("Clear All Reactions"):
    for key in ["reagents", "products", "kf", "kb", "init_conc", "reactions", "species_list", "reaction_tuples"]:
        st.session_state[key] = [] if key != "init_conc" else {}
    st.success("Cleared all reaction data.")