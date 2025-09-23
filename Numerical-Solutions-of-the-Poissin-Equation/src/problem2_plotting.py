import matplotlib.pyplot as plt
import numpy as np

def plot(path):
    x_values = np.zeros(100)
    ux_values = np.zeros(100)


    with open(path, "r") as infile:
        x_list = []
        ux_list = []
        infile.readline()
        for line in infile:
            x_list.append(line.split()[0])
            ux_list.append(line.split()[1])

        for i in range(0, 100):
            x_values[i] = float(x_list[i])
            ux_values[i] = float(ux_list[i])



    fig, ax = plt.subplots()
    ax.plot(x_values, ux_values)
    ax.set_xlabel('x-values')
    ax.set_ylabel('u(x) values')
    ax.set_title('u(x) = 100*e^(-10x) for x in [0, 1]')
    fig.savefig("problem2_ux_plot.pdf")

plot("x_and_ux_values_output.txt")