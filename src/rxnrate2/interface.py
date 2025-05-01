import streamlit as st

# Name of our application for rxnrate
st.title("RxnRate")

# Type of reaction; simple or composed
st.write("Select the type of your reaction")
click1 = st.button("Single reaction")
click2 = st.button("Composite reaction") 

if click1:
    
# Reaction the user wants to monitor
rxn = st.text_input("What is your reaction ? ")

# Type of reaction the user has (A + B -> C, A -> B, ...)
options = ["A + B -> C", " A -> B", "A + B -> C + D", "A -> C + D"]
type_rxn1 = st.selectbox("Select the type of reaction: ", options)

st.write("You have chosen the mode : Reaction of type ")

