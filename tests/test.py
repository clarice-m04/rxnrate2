"test pour voir si j'arrive à push and pull"

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

def linear_rxn():
    st.title("Welcome in linear reaction part")