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


def rxn_diagram_multi(reagent_smiles_list, product_smiles_list):
    try:
        font = ImageFont.truetype("Times New Roman.ttf", 18)
    except:
        font = ImageFont.load_default()

    # Convert SMILES to RDKit molecules
    reagent_mols = [Chem.MolFromSmiles(smi) for smi in reagent_smiles_list if smi]
    product_mols = [Chem.MolFromSmiles(smi2) for smi2 in product_smiles_list if smi2]

    reagent_mols = [mol for mol in reagent_mols if mol]
    product_mols = [mol2 for mol2 in product_mols if mol2]

    if not reagent_mols or not product_mols:
        st.warning("Could not generate one or more molecule images.")
        return None

    mol_size = (150, 150)
    plus_width = 20
    arrow_space = 60
    spacing = 10

    # Create images for molecules
    reactant_imgs = [Draw.MolToImage(mol, size=mol_size) for mol in reagent_mols]
    product_imgs = [Draw.MolToImage(mol, size=mol_size) for mol in product_mols]

    # Total width calculation
    num_reactants = len(reactant_imgs)
    num_products = len(product_imgs)

    reactant_width = num_reactants * mol_size[0] + (num_reactants - 1) * plus_width
    product_width = num_products * mol_size[0] + (num_products - 1) * plus_width
    total_width = reactant_width + product_width + arrow_space + 4 * spacing
    canvas_height = mol_size[1] + 40

    canvas = Image.new("RGB", (total_width, canvas_height), "white")
    draw = ImageDraw.Draw(canvas)

    # Draw reactants with plus signs
    x_offset = spacing
    for i, img in enumerate(reactant_imgs):
        canvas.paste(img, (x_offset, 20))
        x_offset += mol_size[0]
        if i < num_reactants - 1:
            draw.text((x_offset + 5, canvas_height // 2 - 10), "+", fill="black", font=font)
            x_offset += plus_width

    # Draw reaction arrow
    arrow_x_start = x_offset + spacing
    arrow_x_end = arrow_x_start + arrow_space - 10
    arrow_y = canvas_height // 2
    draw.line((arrow_x_start, arrow_y, arrow_x_end, arrow_y), fill="black", width=3)
    draw.polygon([
        (arrow_x_end, arrow_y),
        (arrow_x_end - 10, arrow_y - 5),
        (arrow_x_end - 10, arrow_y + 5)
    ], fill="black")
    x_offset = arrow_x_end + spacing

    # Draw products with plus signs
    for i, img in enumerate(product_imgs):
        canvas.paste(img, (x_offset, 20))
        x_offset += mol_size[0]
        if i < num_products - 1:
            draw.text((x_offset + 5, canvas_height // 2 - 10), "+", fill="black", font=font)
            x_offset += plus_width

    return canvas

if reactants and products:
    reagents_smiles_list = []
    products_smiles_list = []

    for reactant in reactants:
        reagent_smile = get_smiles(reactant)
        reagents_smiles_list.append(reagent_smile)
    for product in products:
        product_smile = get_smiles(product)
        products_smiles_list.append(product)
    st.write(f"The reagents smiles are: {reagents_smiles_list}")
    st.write(f"The products smiles are: {products_smiles_list}")

    rxn_diagram_multi(reagents_smiles_list, products_smiles_list)

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
