import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp


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

def plot_solution(sol, species, filename="reaction_plot.jpg"):
    """plots each speacies' concentration"""
    for i, s in enumerate(species):
        plt.plot(sol.t, sol.y[i], label=s)
    plt.xlabel('Time', fontname='Times New Roman', fontsize=14)
    plt.ylabel('Concentration', fontname='Times New Roman', fontsize=14)
    plt.xticks(fontname= 'Times New Roman')
    plt.yticks(fontname= 'Times New Roman') 
    plt.title('Concentration over time of each species', fontname='Times New Roman', fontsize=14)
    plt.legend(prop={'family': 'Times New Roman'})
    #plt.show()
    #plt.grid(True)
    plt.tight_layout()  # Ensures labels aren't cut off
    plt.savefig(filename, dpi=300)
    plt.close() # Close the plot to avoid display if running in batch mode

def plot_solution_forjup(sol, species):
    """plots each speacies' concentration"""
    for i, s in enumerate(species):
        plt.plot(sol.t, sol.y[i], label=s)
    plt.xlabel('Time', fontname='Times New Roman', fontsize=14)
    plt.ylabel('Concentration', fontname='Times New Roman', fontsize=14)
    plt.xticks(fontname= 'Times New Roman')
    plt.yticks(fontname= 'Times New Roman') 
    plt.title('Concentration over time of each species', fontname='Times New Roman', fontsize=14)
    plt.legend(prop={'family': 'Times New Roman'})
    #plt.show()
    #plt.grid(True)
    plt.tight_layout()  # Ensures labels aren't cut off
    plt.show() 