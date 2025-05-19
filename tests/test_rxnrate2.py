

import numpy as np
from rxnrate2.ODE_linearrxn import build_M_matrix, ode_system, solve_reaction, plot_solution, plot_solution_forjupyter
from rxnrate2.ODE_nonlinear import build_RHS, solve_reactions, plot_solution_nl


def test_build_M_matrix_simple_irreversible():
    species = ['A', 'B']
    reactions = [('A', 'B', 1.0, None)]
    M = build_M_matrix(species, reactions)
    expected = np.array([[-1.0, 0.0],
                         [1.0, 0.0]])
    np.testing.assert_allclose(M, expected)


def test_build_M_matrix_reversible():
    species = ['A', 'B']
    reactions = [('A', 'B', 2.0, 1.0)]
    M = build_M_matrix(species, reactions)
    expected = np.array([[-2.0, 1.0],
                         [2.0, -1.0]])
    np.testing.assert_allclose(M, expected)


def test_ode_system_linear():
    M = np.array([[-1, 0], [1, 0]])
    y = np.array([1.0, 0.0])
    dy = ode_system(0, y, M)
    expected = np.array([-1.0, 1.0])
    np.testing.assert_allclose(dy, expected)


def test_solve_reaction_basic():
    species = ['A', 'B']
    reactions = [('A', 'B', 1.0, None)]
    y0 = [1.0, 0.0]
    sol, M = solve_reaction(species, reactions, y0, t_span=(0, 5), t_eval=np.linspace(0, 5, 50))
    assert sol.success
    assert sol.y.shape == (2, 50)
    # Concentration of A should decrease, B should increase
    assert sol.y[0, -1] < y0[0]
    assert sol.y[1, -1] > y0[1]


def test_build_RHS_irreversible():
    species = ['A', 'B']
    # reaction A -> B with k=1.0 irreversible
    reactions = [(['A'], ['B'], 1.0, None)]
    rhs = build_RHS(species, reactions)
    y = np.array([1.0, 0.0])
    dy = rhs(0, y)
    expected = np.array([-1.0, 1.0])
    np.testing.assert_allclose(dy, expected)


def test_build_RHS_reversible():
    species = ['A', 'B']
    # reaction A <-> B with kf=2.0 and kr=1.0
    reactions = [(['A'], ['B'], 2.0, 1.0)]
    rhs = build_RHS(species, reactions)
    y = np.array([1.0, 1.0])
    dy = rhs(0, y)
    # forward rate = 2 * [A] = 2.0, reverse rate = 1 * [B] = 1.0
    # dA/dt = -forward + reverse = -2 + 1 = -1
    # dB/dt = +forward - reverse = +2 - 1 = +1
    expected = np.array([-1.0, 1.0])
    np.testing.assert_allclose(dy, expected)


def test_solve_reactions_basic():
    species = ['A', 'B']
    reactions = [(['A'], ['B'], 1.0, None)]
    y0 = [1.0, 0.0]
    sol = solve_reactions(species, reactions, y0, t_span=(0, 5), t_eval=np.linspace(0, 5, 50))
    assert sol.success
    assert sol.y.shape == (2, 50)
    # A concentration should decrease, B increase
    assert sol.y[0, -1] < y0[0]
    assert sol.y[1, -1] > y0[1]


def test_plot_solution_creates_file(tmp_path):
    species = ['A', 'B']
    reactions = [('A', 'B', 1.0, None)]
    y0 = [1.0, 0.0]
    sol, _ = solve_reaction(species, reactions, y0)
    filename = tmp_path / "test_plot.jpg"
    plot_solution(sol, species, filename=str(filename))
    assert filename.exists()


def test_plot_solution_forjupyter_runs_without_error():
    species = ['A', 'B']
    reactions = [('A', 'B', 1.0, None)]
    y0 = [1.0, 0.0]
    sol, _ = solve_reaction(species, reactions, y0)
    # Should not raise any exceptions (plot shows and closes)
    plot_solution_forjupyter(sol, species)


def test_plot_solution_nl_runs_without_error():
    species = ['A', 'B']
    reactions = [(['A'], ['B'], 1.0, None)]
    y0 = [1.0, 0.0]
    sol = solve_reactions(species, reactions, y0)
    # Should not raise any exceptions when plotting nonlinear solution
    plot_solution_nl(sol, species)
