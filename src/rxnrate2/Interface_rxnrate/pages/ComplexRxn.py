import streamlit as st
import numpy as np
import pubchempy as pcp
import rdkit as rd

st.set_page_config(page_title="Chemical Reaction Simulator", layout="centered")

from rdkit import Chem
from rdkit.Chem import Draw
import matplotlib.pyplot as plt
from rxnrate2.ODE_nonlinear import solve_reactions, plot_solution
from rxnrate2.Interface_rxnrate.pages.SimpleRxn import get_smiles
from PIL import Image, ImageDraw, ImageFont

st.title("Nonlinear Chemical Reaction Simulator")

# Species input (reagnts and products)
species_input = st.text_input("Enter species (comma-separated)", "H2, O2, H2O")
species = [s.strip() for s in species_input.split(",") if s.strip()]
num_species = len(species)

# Initial concentrations
st.subheader("Initial Concentrations")
initial_conc = []
for s in species:
    conc = st.number_input(f"[{s}]₀", min_value=0.0, value=1.0, key=f"conc_{s}")
    initial_conc.append(conc)

# Reaction input
st.subheader("Define Reactions")
reaction_list = []
num_reactions = st.number_input("Number of reactions", min_value=1, value=1, step=1)

for i in range(num_reactions):
    with st.expander(f"Reaction {i+1}"):
        reactants_str = st.text_input(f"Reactants (comma-separated) - R{i+1}", key=f"reactants_{i}")
        products_str = st.text_input(f"Products (comma-separated) - R{i+1}", key=f"products_{i}")
        kf = st.number_input(f"Forward rate constant kf - R{i+1}", min_value=0.0, value=1.0, key=f"kf_{i}")
        kr = st.number_input(f"Reverse rate constant kr (0 for irreversible) - R{i+1}", min_value=0.0, value=0.0, key=f"kr_{i}")

        reactants = [r.strip() for r in reactants_str.split(",") if r.strip()]
        products = [p.strip() for p in products_str.split(",") if p.strip()]
        kr_val = kr if kr > 0 else None

        reaction_list.append((reactants, products, kf, kr_val))

if reactants and products:
    for reactant in reactants:
        reagent_smile = get_smiles(reactant)
    for product in products:
        product_smile = get_smiles(product)

def rxn_diagram(reagent_smiles, product_smiles):
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

    # Labels and concentrations
    #draw.text((10, 210), f"{reagent_labels}", fill='black', font=font)
    #draw.text((290, 210), f"{product_labels}", fill='black', font=font)
 

    return canvas


rxn_diagram(reagent_smile, product_smile)

# Simulation time
st.subheader("Simulation Time")
t_max = st.number_input("Maximum time", min_value=1.0, value=10.0)
t_eval = np.linspace(0, t_max, 300)

# Run simulation
if st.button("Run Simulation"):
    try:
        sol = solve_reactions(species, reaction_list, initial_conc, t_span=(0, t_max), t_eval=t_eval)

        # Use your custom plot function
        fig = plt.figure()
        plot_solution(sol, species)
        st.pyplot(fig)

    except Exception as e:
        st.error(f"An error occurred: {e}")
