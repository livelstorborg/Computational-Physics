import numpy as np
import matplotlib.pyplot as plt

n = 100
N = n + 1
h = 1/n
a = -1/(h**2)
d = 2/(h**2)
x_hat_analytical = np.linspace(0, 1, n + 1)

numerical_eigvecs = np.loadtxt("numerical_smallest_lambda.txt", usecols=(0, 1, 2), skiprows=1)
N_numerical = numerical_eigvecs.shape[0]
x_hat_numerical = np.linspace(0, 1, N_numerical)

#normalizing vectors
num_eigvec_norm = np.zeros(shape=numerical_eigvecs.shape)
for i in range(3):
    num_eigvec_norm[:, i] = numerical_eigvecs[:, i] / np.max(np.abs(numerical_eigvecs[:, i]))

#multiplying eigenvector 3 with -1 to get the same shape as numerical plot.
for i in range(len(numerical_eigvecs)):
	num_eigvec_norm[i][2] = -1 * num_eigvec_norm[i][2]

for i in range(3):
	#setting boundary points
	num_eigvec_norm[0, i] = 0
	num_eigvec_norm[-1, i] = 0

fig1, ax1 = plt.subplots(3, 1, figsize=(8, 6))
for i in range(3):
	analytical_solution = np.sin((i + 1) * np.pi * x_hat_analytical) if i != 2 else -np.sin((i + 1) * np.pi * x_hat_analytical)
	ax1[i].plot(x_hat_analytical, analytical_solution, label="Analytical", color="red")
	ax1[i].plot(x_hat_numerical, num_eigvec_norm[:, i], label=f"Numerical Eigenvector {i+1}",linestyle="--", marker="x", color='blue')
	ax1[i].legend()
	ax1[i].grid(True)
	ax1[i].set_xlabel("x_hat")
	ax1[i].set_ylabel(f"Eigenvector {i + 1}")
plt.tight_layout()
fig1.savefig("numerical_vs_analytical_" + str(len(numerical_eigvecs)) + ".pdf")