import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.animation import FuncAnimation


""" Collectiong data from CSV files """
potential_file = 'Potential_2.csv'
P_file = 'P_2.csv'
U_real_file = 'U_real_2.csv'
U_imag_file = 'U_imag_2.csv'
prob_file = 'Prob_vec_2.csv'





pot_data = pd.read_csv(potential_file, header=None)
P = pd.read_csv(P_file, header = None)                # P is 2D representing 3D
U_real = pd.read_csv(U_real_file, header = None)     
U_imag = pd.read_csv(U_imag_file, header = None)      
prob = pd.read_csv(prob_file, header = None)

P = P.values
U_real = U_real.values
U_imag = U_imag.values
prob = prob.values[0]

rows, cols = P.shape
M_2 = int(np.sqrt(rows))
M = M_2 + 2
timesteps = cols

# Reshape the 2D array into a 3D array
P_sim = P.reshape((M_2, M_2, timesteps)) 
U_real_sim = U_real.reshape((M_2, M_2, timesteps)) 
U_imag_sim = U_imag.reshape((M_2, M_2, timesteps)) 



"""

All the matrices are within the boundary conditions, and in order for the axes 
to be correct when plotting, we need to add a border of 0's

"""

P = np.zeros((M, M, timesteps))
U_real = np.zeros((M, M, timesteps))
U_imag = np.zeros((M, M, timesteps))

P[1:M_2+1, 1:M_2+1, :] = P_sim
U_real[1:M_2+1, 1:M_2+1, :] = U_real_sim
U_imag[1:M_2+1, 1:M_2+1, :] = U_imag_sim

num_labels = 6  
ticks = np.linspace(0, M, num_labels)  
x_labels = np.linspace(0, 1, num_labels) 
y_labels = np.linspace(1, 0, num_labels) #The y-axis' are inverted





def plot_probabilities(T):

    time = np.linspace(0, T * 1000, timesteps)
    y_min = np.min(prob)
    y_max = np.max(prob)
    plt.plot(time, prob)
    plt.ylim(y_min, y_max)
    plt.xlabel(r'time $\cdot 10^{3}$', fontsize = 16)
    plt.ylabel('probability', fontsize = 16)
    plt.xticks(fontsize = 16)
    plt.yticks(fontsize = 16)
    plt.grid()
    plt.tight_layout()
    plt.savefig("probabilities_with_pot.pdf")
    plt.show()





def plot_potential():

    sns.heatmap(pot_data, cmap='viridis', square=True, xticklabels=ticks, yticklabels=ticks)
    plt.xticks(ticks, [f"{label:.1f}" for label in x_labels], fontsize = 16)  
    plt.yticks(ticks, [f"{label:.1f}" for label in y_labels], fontsize = 16) 
    cbar = plt.gca().collections[-1].colorbar
    cbar.ax.tick_params(labelsize=16)  
    plt.xlabel('x', fontsize = 16)
    plt.ylabel('y', fontsize = 16)
    plt.tight_layout()
    plt.savefig("two_slits.pdf")
    plt.show()



def plot_three_times(data, label):

    'data is either P, U_real, or U_imag'
    t_3_id = [0, int(timesteps / 2), timesteps - 1]
    for i in range(len(t_3_id)):
        sns.heatmap(data[:,:,t_3_id[i]], cmap='viridis', square=True) 
        cbar = plt.gca().collections[-1].colorbar
        cbar.ax.tick_params(labelsize=16)  
        plt.xlabel('x', fontsize = 16)
        plt.ylabel('y', fontsize = 16)
        plt.xticks(ticks, [f"{label:.1f}" for label in x_labels], fontsize = 16)  
        plt.yticks(ticks, [f"{label:.1f}" for label in y_labels], fontsize = 16)
        plt.tight_layout() 
        plt.savefig(f"time_plot_{label}_{t_3_id[i]}_double_slit.pdf")
        plt.show()



def plot_prob_screen(x, t, T):

    """ Takes in the position and time of the screen """

    x_id = int( x * (M - 1) )
    t_id = int( (t / T) * (timesteps - 1) )
    
    p_screen = P[:, x_id, t_id]
    p_screen = p_screen / np.sum(p_screen)  #Needs to be normalized as we expect to find the particle here

    y = np.linspace(0, 1, M)

    plt.plot(y, p_screen)
    plt.ylabel('Probability', fontsize = 16)
    plt.xlabel('y', fontsize = 16)
    plt.xticks(fontsize = 16)
    plt.yticks(fontsize = 16)
    plt.grid()
    plt.tight_layout()
    plt.savefig("Screen_three_slit.pdf")
    plt.show()







def update(frame, data, heatmap, ax, label):

    heatmap.set_data(data[:, :, frame])
    ax.set_title(f"{label} - Time step {frame}", fontsize=16)  # Update the title for each frame
    return heatmap



def animate_heatmap(data, label, timesteps):

    fig, ax = plt.subplots()

    # Setup initial heatmap
    heatmap = ax.imshow(data[:, :, 0], cmap='viridis', interpolation='nearest')
    cbar = fig.colorbar(heatmap, ax=ax)
    cbar.ax.tick_params(labelsize=16)
    
    plt.xticks(ticks, [f"{label:.1f}" for label in x_labels], fontsize = 16)  
    plt.yticks(ticks, [f"{label:.1f}" for label in y_labels], fontsize = 16)
    
    # Plot x and y labels only once
    ax.set_xlabel('x', fontsize=16)
    ax.set_ylabel('y', fontsize=16)
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{label:.1f}" for label in x_labels], fontsize=16)
    ax.set_yticks(ticks)
    ax.set_yticklabels([f"{label:.1f}" for label in y_labels], fontsize=16)
    
    anim = FuncAnimation(fig, update, fargs=(data, heatmap, ax, label),
                         frames=timesteps, interval=200/3, blit=False, repeat=False)
    
    plt.tight_layout()
    plt.show()
    return anim




plot_probabilities(T=0.002)
plot_potential()
plot_three_times(P, 'P')
plot_three_times(U_real, 'real')
plot_three_times(U_imag, 'imag')
plot_prob_screen(x=0.8, t=0.002, T=0.002)

label = 'Animation of Double Slit'
anim = animate_heatmap(P, label, timesteps)
anim.save('animation.mp4')



