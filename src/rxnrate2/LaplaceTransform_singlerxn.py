import sympy as sp
from functools import reduce
from operator import mul

# Define symbols
s, t = sp.symbols('s t')

def calculate_laplace_transforms(reagents, products, k, initial_concentrations):
    """
    Calculates the Laplace transforms for the reactants and products.

    Parameters:
    - reagents: List of reactant names (e.g., ['A', 'B'])
    - products: List of product names (e.g., ['C', 'D'])
    - k: Rate constant (numeric or symbolic)
    - initial_concentrations: Dictionary with initial concentrations for each species

    Returns:
    - A dictionary containing the Laplace transforms for all reactants and products.
    """
    
    # Laplace transforms for reactants: L{A(t)} = A_0 / (s + k)
    reactant_laplace = {
        reagent: initial_concentrations[reagent] / (s + k)
        for reagent in reagents
    }

    # Product of all reactant Laplace transforms: L{P(t)} ∝ ∏ L{R(t)}
    total_laplace_product = reduce(mul, [reactant_laplace[r] for r in reagents])

    # Laplace transforms for products: assume same formation for each (∏ L{R(t)} / s)
    product_laplace = {
        product: total_laplace_product / s
        for product in products
    }

    # Merge both dictionaries
    return {**reactant_laplace, **product_laplace}


def simplify_inverse_laplace(expr):
    """
    Simplifies the inverse Laplace transform expression by replacing Heaviside(t) with 1.

    Parameters:
    - expr: A sympy expression (typically from inverse_laplace_transform)

    Returns:
    - A simplified sympy expression with Heaviside(t) replaced by 1.
    """
    return sp.simplify(expr.subs(sp.Heaviside(t), 1))


def inverse_laplace_transform(laplace_results):
    """
    Calculates the inverse Laplace transform for each Laplace domain expression.

    Parameters:
    - laplace_results: Dictionary of Laplace-transformed species

    Returns:
    - A dictionary with time-domain expressions for each species.
    """
    inverse_results = {}

    for species, laplace_expr in laplace_results.items():
        time_expr = sp.inverse_laplace_transform(laplace_expr, s, t)
        inverse_results[species] = simplify_inverse_laplace(time_expr)

    return inverse_results

reagents = ['water', 'potassium']
products = ['hydrogen peroxide', 'oxygen']
k = 2
initial_concentrations = {'water': 1, 'potassium': 2, 'hydrogen peroxide': 0, 'oxygen': 0}

laplace_res = calculate_laplace_transforms(reagents, products, k, initial_concentrations)
inverse_res = inverse_laplace_transform(laplace_res)

print("Laplace Transforms:")
for name, expr in laplace_res.items():
    print(f"{name}(s) = {expr}")

print("\nInverse Laplace Transforms:")
for name, expr in inverse_res.items():
    print(f"{name}(t) = {expr}")
