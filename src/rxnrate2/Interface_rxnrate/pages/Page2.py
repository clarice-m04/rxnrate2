from LaplaceTransform_singlerxn import calculate_laplace_transforms 

reagents = ['A', 'B', 'C', 'F']
products = ['D', 'E']
k = 2
initial_concentrations = {'A': 1, 'B': 2, 'C': 1, 'D': 0, 'E': 0, 'F': 1.0}

laplace_res = calculate_laplace_transforms(reagents, products, k, initial_concentrations)