import streamlit as st

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



import pandas as pd

# Default table structure
default_data = pd.DataFrame({
    "Species": [""] * 5,
    "k foward": [0.0] * 5,
    "k backward": [0.0] * 5,
    "Initial Concentration": [None] * 5
})

# Load or initialize editable table
table = st.data_editor(
    data=default_data,
    num_rows="dynamic",
    use_container_width=True,
    key="reaction_table"
)

# Convert blank strings to None for concentration
def clean_concentration(value):
    try:
        return float(value)
    except (ValueError, TypeError):
        return None

# Extract cleaned lists
species_list = table["Species"].tolist()
k_list = table["k"].tolist()
concentration_list = [clean_concentration(val) for val in table["Initial Concentration"]]

# Display extracted data
st.subheader("Extracted Lists")
st.write("**Species:**", species_list)
st.write("**k Values:**", k_list)
st.write("**Initial Concentrations:**", concentration_list)

# reaction_model.py




#Draw scheme
def generate_reaction_graphviz(species_list, k_list):
    """
    Generate a Graphviz DOT string for the reaction scheme.
    Assumes a linear chain: A --k1--> B --k2--> C ...
    """
    if len(species_list) < 2 or len(k_list) < 1:
        return "digraph G {\nlabel=\"Not enough data\"\n}"

    lines = ["digraph G {", "rankdir=LR;"]  # left-to-right orientation
    for i in range(len(species_list) - 1):
        source = species_list[i]
        target = species_list[i + 1]
        k = k_list[i] if i < len(k_list) else ""
        lines.append(f'"{source}" -> "{target}" [label="k={k}"];')
    lines.append("}")
    return "\n".join(lines)

