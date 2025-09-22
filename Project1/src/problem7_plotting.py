import numpy as np
import matplotlib.pyplot as plt

def u(x):
	func = (1 - ((1 - np.exp(-10)) * x) - np.exp((-10) * x))
	return func

def plot_problem7(filename_list, n_step_list):
	fig, ax = plt.subplots()

	for i in range(len(n_step_list)):
		v_values = np.zeros(n_step_list[i] + 1)
		x_values = np.zeros(n_step_list[i] + 1)
		ux_values = np.zeros(n_step_list[i] + 1)

		with open(filename_list[i], "r") as infile:
			v_list = []
			x_list = []
			infile.readline()
			for line in infile:
				v_list.append(line.split()[0])
				x_list.append(line.split()[1])

			for j in range(n_step_list[i] + 1):
				v_values[j] = float(v_list[j])
				x_values[j] = float(x_list[j])
				ux_values[j] = u(x_values[j])

			ax.plot(x_values, v_values, "--", label=f"numerical solution n_step = {n_step_list[i]}")

			if i == len(n_step_list) - 1: #only plotting exact solution for n_step = 1000
				ax.plot(x_values, ux_values, label="exact solution", alpha=0.5)

	ax.set_xlabel('x')
	ax.set_ylabel('u(x)')
	ax.set_title(f'Comparison of Numerical Solution and Exact Solution with different n_step')
	ax.legend(loc="upper right", fontsize="small")
	ax.grid(True)
	fig.savefig("plot_comparing_exact_numerical.pdf")
	

filename_list = ["problem7_v_x_10steps.txt", "problem7_v_x_100steps.txt", "problem7_v_x_1000steps.txt", "problem7_v_x_10000steps.txt", "problem7_v_x_100000steps.txt"]
n_step_list = [10, 100, 1000, 10000, 100000]
plot_problem7(filename_list, n_step_list)



