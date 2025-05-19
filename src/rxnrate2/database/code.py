# Se mettre dans le dossier git\rxnrate2\src\rxnrate2\database pour run ce fichier
#(sinon Data_projet.csv est pas trouvé par le code)

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
import sys
import os

path_data = os.path.join(os.path.dirname(__file__), "Data_projet.csv")

data = pd.read_csv(path_data, sep = ";")
df = pd.DataFrame(data)

# Check if your reactants are in the database

def check(smiles1: str, smiles2: str) -> bool | int:
    if smiles1 in df['smiles r1'].values :
        index1 = []
        for a in range(0, len(df), 1):
                if df.loc[a, 'smiles r1'] == smiles1:
                    index1.append(a)
        if smiles2 in df['smiles r2'].values:
            index2 = []
            for b in range(0, len(df), 1):
                if df.loc[b, 'smiles r2'] == smiles2:
                    index2.append(b)
            common_index = list(set(index1) & set(index2))
            if len(common_index) == 1:
                print("Yeah ! Your reaction is in our database ! Let's check its evolution with time !")
                return True, int(''.join(map(str, common_index)))
            else:
                print("There are multiple reactions corresponding to your research, please be more precise.")
                return False, 0
        else:
            print("Reactant 1 is in our database, but there is no reaction with reactant 2 found.")
            return False, 0
            
    elif smiles1 in df['smiles r2'].values :
        index3 = []
        for c in range(0, len(df), 1):
                if df.loc[c, 'smiles r2'] == smiles1:
                    index3.append(c)
        if smiles2 in df['smiles r1'].values:
            index4 = []
            for d in range(0, len(df), 1):
                if df.loc[d, 'smiles r1'] == smiles2:
                    index4.append(d)
            common_index_bis = list(set(index3) & set(index4))
            if len(common_index_bis) == 1:
                print("Yeah ! Your reaction is in our database ! Let's check its evolution with time !")
                return True, int(''.join(map(str, common_index_bis)))
            else:
                print("There are multiple reactions corresponding to your research, please be more precise.")
                return False, 0
        else:
            print("Reactant 2 is not present in this database.")
            return False, 0

    else:
        print("Sorry, the reaction you are looking for in not available in our database.")
        return False, 0

# Calculate k at the chosen temperature

def calc_temperature(temp: float,E: float, A: float) -> float:
    return (A * (np.exp(-E/temp)))


# Canonicalize the smiles of the reactants

def canonicalize_smiles(smiles: str) -> str:
    if not isinstance(smiles, str):
        raise TypeError(f"Invalid type {type(smiles)}: smiles must be a string")
    
    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        raise ValueError(f"Could not convert smiles to mol")

    return Chem.MolToSmiles(mol)

# Show the images to make sure it is the right molecule

def check_molecule(canonicalized_smile: str) -> bool :
    print("Have a look at the molecule and check if it's the one you entered in this program: ")
    reactant = Chem.MolFromSmiles(canonicalized_smile)
    drawing = Draw.MolToImage(reactant)
    drawing.show()
    response = str(input('Please enter "yes" if it is the right molecule, and "no" if it is not.'))

    if response == "Yes" or response == "yes":
        return True
    else:
        return False
    


# Code

compound1 = str(input('Please enter the smiles of the first reactant: '))
compound2 = str(input('Please enter the smiles of the second reactant: '))

def link_database(compound1, compound2, temperature):
    smiles1_can = canonicalize_smiles(compound1)
    smiles2_can = canonicalize_smiles(compound2)

    drawing1 = check_molecule(smiles1_can)
    drawing2 = check_molecule(smiles2_can)

    if drawing1 and drawing2 :
        is_ok = check(smiles1_can, smiles2_can)
        print(is_ok)
    
        if is_ok[0] == True:
            #temperature = float(input('Please enter the temperature of the reaction you will perform in [°C]: '))
            if temperature < 273 :
                temperature += 273     # Conversion en Kelvin si la temperature a bien été donnée en °C (supposant que peu de réaction vont avoir lieu à plus de 273 °C)

            value_k_corrected = calc_temperature(temperature, df.loc[is_ok[1], 'E'], df.loc[is_ok[1], 'Arrhenius factor'])
            return value_k_corrected
            #("%.2E" %value_k_corrected, "L/(mol s)")