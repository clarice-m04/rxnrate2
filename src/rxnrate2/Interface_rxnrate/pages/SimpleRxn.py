import streamlit as st
import pandas as pd

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

def process_reaction_data(species_list, k_forward, k_backward, concentration_list):
    return {
        "Species": species_list,
        "k Forward": k_forward,
        "k Backward": k_backward,
        "Initial Concentrations": concentration_list
    }

def generate_reaction_graphviz(species_list, k_forward, k_backward, concentration_list):
    """
    Draws a reaction graph with forward/backward arrows and shows initial concentration under each species.
    """
    if len(species_list) < 2:
        return "digraph G {\nlabel=\"Not enough species\"\n}"

    lines = ["digraph G {", "rankdir=LR;", "node [shape=box];"]

    # Define nodes with concentrations
    for i, species in enumerate(species_list):
        conc = concentration_list[i]
        conc_str = f"initial concentration={conc}"
        label = f"<<TABLE BORDER='0' CELLBORDER='1' CELLSPACING='0'><TR><TD>{species}</TD></TR><TR><TD>{conc_str}</TD></TR></TABLE>>"
        lines.append(f'"{species}" [label={label}];')

    # Define edges
    for i in range(len(species_list) - 1):
        a = species_list[i]
        b = species_list[i + 1]
        kf = k_forward[i] if i < len(k_forward) else None
        kb = k_backward[i] if i < len(k_backward) else None

        if kf is not None and kb is not None:
            lines.append(f'"{a}" -> "{b}" [label="kf={kf}"];')
            lines.append(f'"{b}" -> "{a}" [label="kb={kb}"];')
        elif kf is not None:
            lines.append(f'"{a}" -> "{b}" [label="kf={kf}"];')
        elif kb is not None:
            lines.append(f'"{b}" -> "{a}" [label="kb={kb}"];')

    lines.append("}")
    return "\n".join(lines)


#from reaction_model import process_reaction_data, generate_reaction_graphviz

st.title("Editable Reaction Table (Reversible Reactions)")

# Default data table
default_data = pd.DataFrame({
    "Species": [None] * 5,
    "k_forward": [None] * 5,
    "k_backward": [None] * 5,
    "Initial Concentration": [0.0] * 5
})

# Display editable table
table = st.data_editor(
    data=default_data,
    num_rows="dynamic",
    use_container_width=True,
    key="reaction_table"
)

# Clean concentration and k values
def parse_float(value):
    try:
        return float(value)
    except (ValueError, TypeError):
        return None

# Extract lists
species_list = table["Species"].tolist()
k_forward_list = [parse_float(val) for val in table["k_forward"]]
k_backward_list = [parse_float(val) for val in table["k_backward"]]
concentration_list = [parse_float(val) for val in table["Initial Concentration"]]

# Show lists
st.subheader("Extracted Lists")
st.write("**Species:**", species_list)
st.write("**k Forward:**", k_forward_list)
st.write("**k Backward:**", k_backward_list)
st.write("**Initial Concentrations:**", concentration_list)

# Run external function
"""if st.button("Run Model"):
    result = process_reaction_data(species_list, k_forward_list, k_backward_list, concentration_list)
    st.subheader("Model Output")
    st.json(result)"""

# Draw reaction scheme
if st.button("Draw Reaction Scheme"):
    dot_string = generate_reaction_graphviz(species_list, k_forward_list, k_backward_list, concentration_list)
    st.graphviz_chart(dot_string)

# reaction_model.py

