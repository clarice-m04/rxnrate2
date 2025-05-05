import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import sympy as sp

def build_M_matrix(species, reactions):
    """Constructs matrix M from the rxn definitions
    species: list like ['A', 'B', 'C']
    reactions: tuple like [('A', 'B', 1.0, None) or ('B', 'C', 0.5, 0.3)]
                first of tuple is the reagent/from species, second is product/to species, 3rd is k forward, last is a k reverse or none"""
    n = len(species)
    species_idx = {s: i for i, s in enumerate(species)}
    M = np.zeros((n, n))

    for from_s, to_s, kf, kr in reactions:
        i = species_idx[from_s]
        j = species_idx[to_s]
        M[i, i] -= kf
        M[j, i] += kf
        if kr is not None:
            M[j, j] -= kr
            M[i, j] += kr

    return M

def ode_system(t, y, M):
    """define the derivative of vector y: dy/dt= My"""
    return M @ y

def solve_reaction(species, reactions, y0_vals, t_span=(0, 10), t_eval=None):
    """Builds M from inputs
    sets up initial concentrations y_0
    uses RK45 with solve_ivp"""
    M = build_M_matrix(species, reactions)
    y0 = np.array(y0_vals)
    if t_eval is None:
        t_eval = np.linspace(t_span[0], t_span[1], 200)
    sol = solve_ivp(ode_system, t_span, y0, args=(M,), method='RK45', t_eval=t_eval)
    return sol, M

def plot_solution(sol, species):
    """plots each speacies' concentration"""
    for i, s in enumerate(species):
        plt.plot(sol.t, sol.y[i], label=s)
    plt.xlabel('Time')
    plt.ylabel('Concentration')
    plt.title('Concentration over time of each species')
    plt.legend()
    #plt.grid(True)
    plt.show()

#def symbolic_solution(M_numeric, y0_numeric):
#    """gives symbolic function for each concentration"""
#     t = sp.Symbol('t')
#     M = sp.Matrix(M_numeric)
#    y0 = sp.Matrix(y0_numeric)
#    yt = (M * t).exp() * y0
#    return yt, t

# === Example Usage ===
species_test = ['A', 'B', 'C', 'D']
reactions_test = [
    ('A', 'B', 1.0, None),
    ('B', 'C', 0.5, 0.3),
    ('C', 'A', 0.2, None),
    ('C', 'D', 2.0, 0.4)
]
y0_vals = [1.0, 0.0, 0.0, 0.0]  # Initial concentrations for A, B, C, D

sol, M = solve_reaction(species_test, reactions_test, y0_vals)


# Get symbolic expression
#yt, t_sym = symbolic_solution(M, y0_vals)
#print("Symbolic expression for concentrations:")
#for i, s in enumerate(species):
#    print(f"{s}(t) = {yt[i]}")

plot_solution(sol, species_test)

#NB: je pourrais rajouter l'evaluation de la concentration a un temps t donné et qu'il soit marqué sur le graphe