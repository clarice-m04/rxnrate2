import streamlit as st
import numpy as np
from rxnrate2.ODE_nonlinear import solve_reactions, plot_solution
import matplotlib.pyplot as plt

st.set_page_config(page_title="Chemical Reaction Simulator", layout="centered")

st.title("Nonlinear Chemical Reaction Simulator")

# Species input
species_input = st.text_input("Enter species (comma-separated)", "A, B, C, D")
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
