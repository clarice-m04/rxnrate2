import streamlit as st
import numpy as np
import pubchempy as pcp
import os
import pandas as pd

st.set_page_config(page_title="Chemical Reaction Simulator", layout="centered")

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

from rxnrate2.ODE_nonlinear import plot_solution_nl_save, solve_reactions
from rdkit import Chem
from rdkit.Chem import Draw
import matplotlib.pyplot as plt
from PIL import Image, ImageDraw, ImageFont
from rxnrate2.Interface_rxnrate.__init__ import set_background

set_background("rxnrate.jpg")

st.title("Nonlinear Chemical Reaction Simulator")


######################################### Functions #################################################

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
            return compound[0].canonical_smiles
    except:
        pass
    return fallback_smiles.get(query.strip(), None)


# Check if your reactants are in the database
path_database = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))),"database/Data_projet.csv")

data = pd.read_csv(path_database, sep = ";")
df = pd.DataFrame(data)

def check(reactant1, reactant2) -> bool | int:
    pb1 = get_smiles(reactant1)
    pb2 = get_smiles(reactant2)

    mol1 = Chem.MolFromSmiles(pb1)
    mol2 = Chem.MolFromSmiles(pb2)

    smiles_can1 = Chem.MolToSmiles(mol1)
    smiles_can2 = Chem.MolToSmiles(mol2)

    if smiles_can1 in df['smiles r1'].values :
        index1 = []
        for a in range(0, len(df), 1):
                if df.loc[a, 'smiles r1'] == smiles_can1:
                    index1.append(a)
        if smiles_can2 in df['smiles r2'].values:
            index2 = []
            for b in range(0, len(df), 1):
                if df.loc[b, 'smiles r2'] == smiles_can2:
                    index2.append(b)
            common_index = list(set(index1) & set(index2))
            if len(common_index) == 1:                
                return True, int(''.join(map(str, common_index)))
            else:
                return False, 0
        else:
            return False, 0
            
    elif smiles_can1 in df['smiles r2'].values :
        index3 = []
        for c in range(0, len(df), 1):
                if df.loc[c, 'smiles r2'] == smiles_can1:
                    index3.append(c)
        if smiles_can2 in df['smiles r1'].values:
            index4 = []
            for d in range(0, len(df), 1):
                if df.loc[d, 'smiles r1'] == smiles_can2:
                    index4.append(d)
            common_index_bis = list(set(index3) & set(index4))
            if len(common_index_bis) == 1:
                return True, int(''.join(map(str, common_index_bis)))
            else:
                return False, 0
        else:
            return False, 0

    else:
        return False, 0

# Calculate k at the chosen temperature

def calc_temperature(temp: float,E: float, A: float) -> float:
    return (A * (np.exp(-E/temp)))

# Draw the reaction 

def rxn_diagram_multi(reagents_list, products_list, kf, kb):
    try:
        font = ImageFont.truetype("Times New Roman.ttf", 18)
    except:
        font = ImageFont.load_default()

    # Create new lists with the canonical smiles of reactants and products
    reagents_smiles_list: list = []
    products_smiles_list: list = []

    for r in reagents_list:
        reagents_smiles_list.append(get_smiles(r))
    for p in products_list:
        products_smiles_list.append(get_smiles(p))
    st.write(f"The reagents SMILES are : {reagents_smiles_list}")
    st.write(f"The products SMILES are : {products_smiles_list}")


    # Convert SMILES to RDKit Molecules

    reagents_molecule_list = []
    products_molecule_list = []

    for r_smile in reagents_smiles_list:
        reagents_molecule_list.append(Chem.MolFromSmiles(r_smile))
    
    for p_smile in products_smiles_list:
        products_molecule_list.append(Chem.MolFromSmiles(p_smile))

    if (not reagents_molecule_list) or (not products_molecule_list):
        st.warning("Could not generate one or more molecule images.")
        return None

    # Creation of the molecules images list
    
    reagents_images_list: list = []
    products_images_list: list = []

    for r_molecule in reagents_molecule_list:
        reagents_images_list.append(Draw.MolToImage(r_molecule, size=(150,150)))
    for p_molecule in products_molecule_list:
        products_images_list.append(Draw.MolToImage(p_molecule, size=(150,150)))

    ## Drawing on the interface

    # Total width calculation
    num_reagents = len(reagents_images_list)
    num_products = len(products_images_list)

    # Basis of the canevas dimension
    plus_width = 20
    arrow_space = 100
    spacing = 10
    mol_size = (150,150)

    reagent_width = num_reagents * mol_size[0] + (num_reagents - 1) * plus_width
    product_width = num_products * mol_size[0] + (num_products - 1) * plus_width
    total_width = reagent_width + product_width + arrow_space + 4 * spacing
    canvas_height = mol_size[1] + 40

    canvas = Image.new("RGB", (total_width, canvas_height), "white")
    draw = ImageDraw.Draw(canvas)

    # Draw reagents with plus signs
    x_offset = spacing
    for i, img in enumerate(reagents_images_list):
        canvas.paste(img, (x_offset, 20))
        x_offset += mol_size[0]
        if i < num_reagents - 1:
            draw.text((x_offset + 6, canvas_height // 2 - 16), "+", fill="black", font_size=24)
            x_offset += plus_width

    # Draw reaction arrow

    if kb == None:
        arrow_x_start = x_offset + spacing
        arrow_x_end = arrow_x_start + arrow_space - 10
        arrow_y = canvas_height // 2
        draw.line((arrow_x_start, arrow_y, arrow_x_end, arrow_y), fill="black", width=3)
        draw.polygon([
            (arrow_x_end + 5, arrow_y),
            (arrow_x_end - 10, arrow_y - 5),
            (arrow_x_end - 10, arrow_y + 5)
        ], fill="black")

        draw.text((arrow_x_start, arrow_y - 26), f"kf = {kf}", fill="black", font_size=16)

        x_offset = arrow_x_end + spacing
    else:
        arrow_up_x_start = x_offset + spacing
        arrow_up_x_end = arrow_up_x_start + arrow_space - 10
        arrow_up_y = canvas_height // 2
        draw.line((arrow_up_x_start, arrow_up_y, arrow_up_x_end - 1, arrow_up_y), fill="black", width=2)
        draw.polygon([
            (arrow_up_x_end, arrow_up_y + 1),
            (arrow_up_x_end - 10, arrow_up_y - 6),
            (arrow_up_x_end - 10, arrow_up_y)
        ], fill="black")

        draw.text((arrow_up_x_start, arrow_up_y - 28), f"kf = {kf}", fill="black", font_size=16)

        arrow_down_x_start = x_offset + spacing
        arrow_down_x_end = arrow_down_x_start + arrow_space - 10
        arrow_down_y = (canvas_height // 2)  + 5
        draw.line((arrow_down_x_start, arrow_down_y, arrow_down_x_end, arrow_down_y), fill="black", width=2)
        draw.polygon([
            (arrow_down_x_start, arrow_down_y + 1),
            (arrow_down_x_start + 10, arrow_down_y),
            (arrow_down_x_start + 10, arrow_down_y + 6)
        ], fill="black")

        draw.text((arrow_down_x_start, arrow_down_y + 12), f"kb = {kb}", fill="black", font_size=16)

        x_offset = arrow_down_x_end + spacing

        

    # Draw products with plus signs
    for i, img in enumerate(products_images_list):
        canvas.paste(img, (x_offset, 20))
        x_offset += mol_size[0]
        if i < num_products - 1:
            draw.text((x_offset + 5, canvas_height // 2 - 10), "+", fill="black", font=font)
            x_offset += plus_width

    return canvas
    

#####################################################################################################

# Species input (reagents and products)
species_input = st.text_input("Enter species (comma-separated)", "H2, N2, NH3")
species = [s.strip() for s in species_input.split(",") if s.strip()]
num_species = len(species)

# Initial concentrations
st.subheader("Initial Concentrations")
initial_conc = []
for s in species:
    conc = st.number_input(f"[{s}]₀", min_value=0.0, value=1.0, format="%.5f", key=f"conc_{s}")
    initial_conc.append(conc)

# Reaction input
st.subheader("Define Reactions")
reaction_list = []
num_reactions = st.number_input("Number of reactions", min_value=1, value=1, step=1)

for i in range(num_reactions):
    with st.expander(f"Reaction {i+1}"):
        reactants_str = st.text_input(f"Reactants (comma-separated) - Reaction{i+1}", key=f"reactants_{i}")
        products_str = st.text_input(f"Products (comma-separated) - Reaction{i+1}", key=f"products_{i}")
        
        temperature_reaction = st.number_input(f"Temperature of the reaction {i+1} [°C]")

        reactants = [r.strip() for r in reactants_str.split(",") if r.strip()]
        products = [p.strip() for p in products_str.split(",") if p.strip()]
        
        if len(reactants) == 2:
            is_in_data = check(reactants[0], reactants[1])
            print(is_in_data)
            if is_in_data[0] == True:
                kf_value = calc_temperature(temperature_reaction, df.loc[is_in_data[1], 'E'], df.loc[is_in_data[1], 'Arrhenius factor'])
                kf_written = "%.2E" %kf_value
                kf_scheme = kf_written
                st.write(f"The value of kf was found using the database, kf = {kf_scheme}")
            else:
                kf_value = st.number_input(f"Forward rate constant kf - Reaction{i+1}", min_value=0.0, value=1.0, format="%.3f", key=f"kf_{i}")
                kf_scheme = kf_value
                st.warning("No reaction found")
        else:
            kf_value = st.number_input(f"Forward rate constant kf - Reaction{i+1}", min_value=0.0, value=1.0, format="%.3f", key=f"kf_{i}")
            kf_scheme = kf_value
            st.warning("Lenght != 2")

        kb = st.number_input(f"Backward rate constant kb (0 for irreversible) - Reaction{i+1}", min_value=0.0, value=0.0, format="%.3f", key=f"kr_{i}")
        kb_val = kb if kb > 0 else None

        reaction_list.append((reactants, products, kf_value, kb_val))


    if reactants and products:
        image = rxn_diagram_multi(reactants, products, kf_scheme, kb_val)
        st.image(image)
    else:
        st.warning("Please insert both reactants and products for your reaction")


# Simulation time
st.subheader("Simulation Time")
t_max = st.number_input("Maximum time", min_value=1.0, value=10.0)
t_eval = np.linspace(0, t_max, 300)

# Run simulation
if st.button("Run Simulation"):
    try:
        sol = solve_reactions(species, reaction_list, initial_conc, t_span=(0, t_max), t_eval=t_eval)

        #look for the folder to put the graph in
        try:
            os.mkdir("figures")
        except FileExistsError:  # directory already exists
            print(f"Directory for figures already exists")
        except Exception as error:  # other error
            print(f"An Error occured: {error}")

        filename = f"./figures/complex_{reaction_list}.jpg"

        # Use your custom plot function
        fig = plt.figure()
        plot_solution_nl_save(sol, species, filename)
        st.subheader("Concentration Plot")
        st.pyplot(fig)

    except Exception as e:
        st.error(f"An error occurred: {e}")


### Ajouter un bouton "Clear all" qui remet tout à 0
### Faire en sorte que ca affiche toujours la premiere reaction meme si in en ajoute une deuxieme