"test pour voir si j'arrive à push and pull"
import numpy as np
from rxnrate2 import ODE_linearrxn, ODE_nonlinear

def test_linear_solution_dimensions():
    species = ['A', 'B']
    reactions = [('A', 'B', 1.0, None)]
    y0 = [1.0, 0.0]
    sol, M = ODE_linearrxn.solve_reaction(species, reactions, y0)
    assert sol.y.shape[0] == len(species)
    assert sol.t.shape[0] == sol.y.shape[1]
    assert np.allclose(sol.y[:, 0], y0, atol=1e-8)

def test_linear_mass_conservation():
    species = ['A', 'B']
    reactions = [('A', 'B', 1.0, None)]
    y0 = [1.0, 0.0]
    sol, _ = ODE_linearrxn.solve_reaction(species, reactions, y0)
    total = sol.y[0] + sol.y[1]
    np.testing.assert_allclose(total, np.full_like(total, 1.0), atol=1e-6)

def test_nonlinear_solution_dimensions():
    species = ['A', 'B', 'C']
    reactions = [ (['A', 'B'], ['C'], 1.0, None) ]
    y0 = [1.0, 1.0, 0.0]
    sol = ODE_nonlinear.solve_reactions(species, reactions, y0)
    assert sol.y.shape[0] == len(species)
    assert sol.t.shape[0] == sol.y.shape[1]
    assert np.allclose(sol.y[:, 0], y0, atol=1e-8)

def test_nonlinear_conservation_like_behavior():
    species = ['A', 'B', 'C']
    reactions = [ (['A', 'B'], ['C'], 1.0, None) ]
    y0 = [1.0, 1.0, 0.0]
    sol = ODE_nonlinear.solve_reactions(species, reactions, y0)
    total = sol.y[0] + sol.y[1] + sol.y[2]
    np.testing.assert_allclose(total, np.full_like(total, 2.0), atol=1e-6)
