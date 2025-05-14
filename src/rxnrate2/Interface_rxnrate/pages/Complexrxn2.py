import streamlit as st
import numpy as np
import os
from rxnrate2.ODE_nonlinear import solve_reactions, plot_solution  # Assuming you saved the nonlinear logic here


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

# Title and intro
st.title("Nonlinear Reaction Simulator")
st.markdown("Define nonlinear reactions with multiple reactants and products.")

# Initialize session state
if 'reaction_inputs' not in st.session_state:
    st.session_state.reaction_inputs = []
if 'reactants' not in st.session_state:
    st.session_state.reactants = [[]]  # list of lists
if 'products' not in st.session_state:
    st.session_state.products = [[]]  # list of lists


## Subtitle
st.title("Chemical Reaction Collector")

col1, col2 = st.columns(2)

with col1:
    # Reactants
    st.markdown("#### Reactants")
    for i in range(len(st.session_state.reactants)):
        st.session_state.reactants[i] = st.text_input(f"Reactant {i+1}", key=f"reactant_{i}")
    kf = st.number_input("Forward Rate Constant (kf)", min_value=0.0, value=1.0, format="%.5f")

    if st.button("Add Reactant Field"):
        st.session_state.reactants.append([])


with col2:
# Products
    st.markdown("#### Products")
    for i in range(len(st.session_state.products)):
        st.session_state.products[i] = st.text_input(f"Product {i+1}", key=f"product_{i}")
    kr = st.number_input("Reverse Rate Constant (kr, optional)", min_value=0.0, value=0.0, format="%.5f")
    kr = kr if kr != 0.0 else None

    if st.button("Add Product Field"):
        st.session_state.products.append([])



# Add reaction button
if st.button("✅ Add Reaction"):
    # Filter out empty entries
    reactants = [r.strip() for r in st.session_state.reactants if r.strip()]
    products = [p.strip() for p in st.session_state.products if p.strip()]

    if reactants and products:
        st.session_state.reaction_inputs.append((reactants, products, kf, kr))
        st.success(f"Reaction added: {' + '.join(reactants)} ⇌ {' + '.join(products)}")

        # Reset input fields
        st.session_state.reactants = [[]]
        st.session_state.products = [[]]
    else:
        st.error("Please enter at least one reactant and one product.")

# Clear all
if st.button("🗑️ Clear All Reactions"):
    st.session_state.reaction_inputs = []
    st.session_state.reactants = [[]]
    st.session_state.products = [[]]
    st.success("All reactions cleared.")

# Show added reactions
st.subheader("Stored Reactions")
for idx, (reactants, products, kf, kr) in enumerate(st.session_state.reaction_inputs):
    arrow = "⇌" if kr is not None else "→"
    st.write(f"{idx + 1}. {' + '.join(reactants)} {arrow} {' + '.join(products)} (kf={kf}, kr={kr if kr else 'N/A'})")

# Compute and plot
st.subheader("Initial Concentrations")
all_species = sorted(set(s for r in st.session_state.reaction_inputs for lst in (r[0], r[1]) for s in lst))
initial_concs = {}
for specie in all_species:
    initial_concs[specie] = st.number_input(f"[{specie}]₀", min_value=0.0, value=1.0, format="%.3f")

if st.button("🚀 Simulate Reaction"):
    if not st.session_state.reaction_inputs:
        st.error("No reactions defined.")
    else:
        y0 = [initial_concs[sp] for sp in all_species]
        try:
            sol = solve_reactions(all_species, st.session_state.reaction_inputs, y0)
            plot_solution(sol, all_species)
        except Exception as e:
            st.error(f"Simulation failed: {e}")

