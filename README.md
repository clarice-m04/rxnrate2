![Project Logo](assets/banner.png)

![Coverage Status](assets/coverage-badge.svg)

<h1 align="center">
rxnrate2
</h1>

<br>



Calculates and graphs the way concentrations evolve in a chemical reaction system using ODEs.  
Supports both linear (matrix-based) and nonlinear (mass-action) kinetics, with options for interactive plotting and a Streamlit GUI.

---

##  Usage

### Linear (first-order) system example

```python
from rxnrate2.linear import build_M_matrix, solve_reaction, plot_solution

M = build_M_matrix([
    ("A", "B", 1.0),
    ("B", "C", 0.5),
    ("C", "B", 0.2),
])
y0 = [0.5, 0.0, 0.0]  # initial concentrations
t_span = (0, 20)

sol = solve_reaction(M, y0, t_span)
plot_solution(sol)
```

### Nonlinear (mass-action) system example

```python
from rxnrate2.nonlinear import solve_reactions, plot_solution

reactions = [
    ("A + B", "C", 1.0),
    ("C", "A + B", 0.5),
    ("C", "A", 0.2),
    ("A + B + C", "D + E", 0.5),
    ("D + E", "A + B + C", 0.5),
]
initial_conc = {"A": 1.0, "B": 1.0, "C": 0.0, "D": 0.0, "E": 0.2}
t_span = (0, 30)

sol = solve_reactions(reactions, initial_conc, t_span)
plot_solution(sol)
```

### GUI

You can also launch the interactive interface:

```bash
streamlit run app.py
```

---

##  Installation

Create a new environment (recommended), you may give it any name:

```bash
conda create -n rxnrate2 python=3.10
conda activate rxnrate2
```

Install the package locally:

```bash
pip install .
```

If you're working in a notebook:

```bash
pip install jupyterlab
```

To install full dependencies including the GUI:

```bash
pip install ".[full]"
```

---

## 📁 Project Structure

```
rxnrate2/
├── linear/         # Matrix-based linear ODE solver
├── nonlinear/      # Mass-action nonlinear ODE solver
├── app.py          # Streamlit GUI
├── examples/       # Jupyter notebooks and tests
└── ...
```

---

## 🧪 Examples

See the `examples/` folder for real usage of:

- Reversible and irreversible reactions
- Linear and nonlinear systems
- Chained and complex multi-step mechanisms

