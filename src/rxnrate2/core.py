from singlerxn1 import calculate_laplace_transforms

reagents = ['water', 'potassium']
products = ['hydrogen peroxide', 'oxygen']
k = 2
initial_concentrations = {'water': 1, 'potassium': 2, 'hydrogen peroxide': 0, 'oxygen': 0}

laplace_res = calculate_laplace_transforms(reagents, products, k, initial_concentrations)
print(laplace_res)