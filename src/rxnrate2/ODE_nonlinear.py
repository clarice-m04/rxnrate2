import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

def build_RHS(species, reactions):
    """Constructs the right-hand side of the ODE dy/dt = f(t,y).
    input:  species: list like ['A', 'B', 'C']
            reactions: tuple like [['A','B'], 1.0, None) or ('B', 'C', 0.5, 0.3)]
                first of tuple is the reagent/from species, second is product/to species, 3rd is k forward, last is a k reverse or none
    returns: function f(t,y) that computes the rate of change of each species"""
   
    species_idx = {s: i for i, s in enumerate(species)}

    def f(t,y):
        dy= np.zeros(len(species))

        for reactants, products, kf, kr in reactions:
            #forward rate
            rate_fwd = kf
            for r in reactants:
                rate_fwd *= y[species_idx[r]]
            
            #apply forward stoichiometry
            for r in reactants:
                dy[species_idx[r]] -= rate_fwd
            for p in products:
                dy[species_idx[p]] += rate_fwd
            
            #if reversible
            if kr is not None:
                rate_rev = kr
                for p in products:
                    rate_rev *= y[species_idx[p]]
                for p in products:
                    dy[species_idx[p]] -= rate_rev
                for r in reactants:
                    dy[species_idx[r]] += rate_rev
        return dy
    return f


def solve_reactions(species, reactions, y0_vals, t_span=(0, 10), t_eval=None):
    """solves the system of ODEs over a time interval with solve_ivp
    inputs: species: list of species names ['A', 'B', 'C']
            reactions: list of reactions [['A','B'], 1.0, None) or ('B', 'C', 0.5, 0.3)]
            y0_vals: initial concentrations [1.0, 0.5, 0.0]
            t_span: time interval (start, end)
            t_eval: time points to evaluate the solution at (optional)
    returns: solution object with time points and concentrations
    """
    rhs = build_RHS(species, reactions)
    if t_eval is None:
        t_eval = np.linspace(t_span[0], t_span[1], 300)
    sol = solve_ivp(rhs, t_span, y0_vals, method='RK45', t_eval=t_eval)
    return sol

def plot_solution(sol, species):
    """plots each species' concentration over time"""
    for i, s in enumerate(species):
        plt.plot(sol.t, sol.y[i], label=s)
    plt.xlabel('Time', fontname= 'Times New Roman')
    plt.ylabel('Concentration', fontname= 'Times New Roman')
    plt.xticks(fontname= 'Times New Roman')
    plt.yticks(fontname= 'Times New Roman')
    plt.title('Species concentration over time', fontname= 'Times New Roman')
    plt.legend(prop={'family': 'Times New Roman'})
    plt.show()


# Example usage
species_test = ['A', 'B', 'C', 'D']
reactions_test = [
    (['A', 'B'], ['C'], 1.0, 0.5),  # reversible
    (['C'], ['A', 'B'], 0.1, None), #irreversible
    (['A', 'C', 'B'], ['D', 'B'], 0.5, 0.5)  #multiple 
]
y0_vals = [1.0, 1.0, 0.0, 0.0]

sol = solve_reactions(species_test, reactions_test, y0_vals)
plot_solution(sol, species_test)
