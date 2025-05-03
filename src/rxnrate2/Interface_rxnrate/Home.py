import streamlit as st
from pages import Page1, Page2

st.markdown(""" <style> html, body, [class^="css"] {font-family : 'Times New Roman', Times, serif !important;} h1, h2, h3, h4, h5, h6 
            {font-family : 'Times New Roman', Times, serif !important;} </style>""", unsafe_allow_html=True)

# Initialisation of the active page
if "page" not in st.session_state:
    st.session_state.page = "home"

def main():
    st.header("Welcome to Rxnrate!")
    st.subheader("Are you curious to see how your reaction advances ?")
    st.write("This program provides you the progression of your reaction over time")


if st.session_state.page == "home":
   main()
    
col1, col2 = st.columns(2)
with col1:
     if st.button("Simple reaction"):
        st.session_state.page = "page1"
with col2:
    if st.button("Composite reaction"):
        st.session_state.page = "page2"

    elif st.session_state.page == "page1":
       
    
    elif st.session_state.page == "page2":
      
