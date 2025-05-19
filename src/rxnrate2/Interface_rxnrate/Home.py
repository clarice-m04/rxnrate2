import streamlit as st
import base64
from pathlib import Path


#from rxnrate2.Interface_rxnrate.pages import SimpleRxn, ComplexRxn

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
    }s
    </style>
""", unsafe_allow_html=True)

# backround of the application

def set_background(jpg_file):
    current_dir = Path(__file__).parent
    image_path = current_dir / jpg_file

    with image_path.open("rb") as image_file:
        encoded = base64.b64encode(image_file.read()).decode()

    st.markdown(
        f"""
        <style>
        .stApp {{
            background: linear-gradient(rgba(255, 255, 255, 0.7), rgba(255, 255, 255, 0.7)),
                        url("data:image/jpg;base64,{encoded}");
            background-size: cover;
            background-repeat: no-repeat;
            background-attachment: fixed;
        }}
        </style>
        """,
        unsafe_allow_html=True
    )

set_background("rxnratepic.jpg")

# Initialisation of the active page
if "page" not in st.session_state:
    st.session_state.page = "home"

def home_page():
    st.header("Welcome to RxnRate!📈✨")
    st.subheader("Are you curious to see how your reaction advances ?")
    st.write("This program provides you the progression of your reaction over time")
    st.write("Now let's see what kind of reaction you have! " \
    "To do this, read the two categories below carefully and select the option that best describes your system. " \
    "Click on the button and you'll be taken straight to the perfect page for your model!")
    


if st.session_state.page == "home":
   home_page()

st.divider()

col1, col2 = st.columns(2)
with col1:
    st.subheader("**Main characteristics:**")
    st.write("- Reaction of type A -> B -> C")
    st.write("- One reactant gives one product")
    st.page_link("pages/SimpleRxn.py", label="Simple reaction", icon="1️⃣")

with col2:
    st.subheader("**Main characteristics:**")
    st.write("- Reaction of type A + B -> C")
    st.write("- Two or more reactants give one product")
    st.page_link("pages/ComplexRxn2.py", label="More complex reaction", icon="2️⃣")
    
# col1, col2 = st.columns(2)
#with col1:
     #if st.button("Single reaction"):
     #   st.session_state.page = "ODE_singlerxn_int"
   
#with col2:
    #if st.button("Composit reaction"):
        #st.session_state.page = "ODE_linearrxn_int"

# Navigation

#if st.session_state.page == "ODE_singlerxn_int":
    #ODE_singlerxn_int.single_rxn()  
#elif st.session_state.page == "ODE_linearrxn_int":
    #ODE_linearrxn_int.linear_rxn()
