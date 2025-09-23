import numpy as np
import matplotlib.pyplot as plt


# ---------- Reading variables ----------
time_analytical, x_analytical, y_analytical, z_analytical = np.loadtxt("xyz_analytical.txt", unpack=True, skiprows=1)

time_RK4_1, x_RK4_1, y_RK4_1, z_RK4_1 = np.loadtxt("xyz_RK4_1.txt", unpack=True, skiprows=1)
time_RK4_1_velocity, x_RK4_1_velocity, y_RK4_1_velocity, z_RK4_1_velocity = np.loadtxt("xyz_RK4_1_velocity.txt", unpack=True, skiprows=1)

time_RK4_2, x_RK4_2, y_RK4_2, z_RK4_2 = np.loadtxt("xyz_RK4_2.txt", unpack=True, skiprows=1)
time_RK4_2_velocity, x_RK4_2_velocity, y_RK4_2_velocity, z_RK4_2_velocity = np.loadtxt("xyz_RK4_2_velocity.txt", unpack=True, skiprows=1)

time_RK4_interactions_1, x_RK4_interactions_1, y_RK4_interactions_1, z_RK4_interactions_1 = np.loadtxt("xyz_RK4_interactions_1.txt", unpack=True, skiprows=1)
time_RK4_interactions_1_velocity, x_RK4_interactions_1_velocity, y_RK4_interactions_1_velocity, z_RK4_interactions_1_velocity = np.loadtxt("xyz_RK4_interactions_1_velocity.txt", unpack=True, skiprows=1)

time_RK4_interactions_2, x_RK4_interactions_2, y_RK4_interactions_2, z_RK4_interactions_2 = np.loadtxt("xyz_RK4_interactions_2.txt", unpack=True, skiprows=1)
time_RK4_interactions_2_velocity, x_RK4_interactions_2_velocity, y_RK4_interactions_2_velocity, z_RK4_interactions_2_velocity = np.loadtxt("xyz_RK4_interactions_2_velocity.txt", unpack=True, skiprows=1)

time_FE_1, x_FE_1, y_FE_1, z_FE_1 = np.loadtxt("xyz_FE_1.txt", unpack=True, skiprows=1)
time_FE_1_interactions, x_FE_1_interactions, y_FE_1_interactions, z_FE_1_interactions = np.loadtxt("xyz_FE_1_interactions.txt", unpack=True, skiprows=1)

time_FE_2, x_FE_2, y_FE_2, z_FE_2 = np.loadtxt("xyz_FE_2.txt", unpack=True, skiprows=1)
time_FE_2_interactions, x_FE_2_interactions, y_FE_2_interactions, z_FE_2_interactions = np.loadtxt("xyz_FE_2_interactions.txt", unpack=True, skiprows=1)

time_4000, error_RK4_4000 = np.loadtxt('relative_error_RK4_4000.txt', unpack=True, skiprows=1)
time_8000, error_RK4_8000 = np.loadtxt('relative_error_RK4_8000.txt', unpack=True, skiprows=1)
time_16000, error_RK4_16000 = np.loadtxt('relative_error_RK4_16000.txt', unpack=True, skiprows=1)
time_32000, error_RK4_32000 = np.loadtxt('relative_error_RK4_32000.txt', unpack=True, skiprows=1)

time_4000, error_FE_4000 = np.loadtxt('relative_error_FE_4000.txt', unpack=True, skiprows=1)
time_8000, error_FE_8000 = np.loadtxt('relative_error_FE_8000.txt', unpack=True, skiprows=1)
time_16000, error_FE_16000 = np.loadtxt('relative_error_FE_16000.txt', unpack=True, skiprows=1)
time_32000, error_FE_32000 = np.loadtxt('relative_error_FE_32000.txt', unpack=True, skiprows=1)

particles_07_non, omega_v_07_non = np.loadtxt('txt_files/f_fine_non_interaction_0.700000.txt', unpack=True, skiprows=1)
particles_07_interaction, omega_v_07_interaction = np.loadtxt('txt_files/f_fine_interaction_0.700000.txt', unpack=True, skiprows=1)

particles_01_non, omega_v_01_non = np.loadtxt('txt_files/f_fine_non_interaction_0.100000.txt', unpack=True, skiprows=1)
particles_01_interaction, omega_v_01_interaction = np.loadtxt('txt_files/f_fine_interaction_0.100000.txt', unpack=True, skiprows=1)

particles_01, omega_v_01 = np.loadtxt('txt_files/f0.100000.txt', unpack=True, skiprows=1)
particles_04, omega_v_04 = np.loadtxt('txt_files/f0.400000.txt', unpack=True, skiprows=1)
particles_07, omega_v_07 = np.loadtxt('txt_files/f0.700000.txt', unpack=True, skiprows=1)


# ---------- Plotting - z(t) ----------
plt.figure(figsize=(10, 6))
plt.plot(time_RK4_1, z_RK4_1, label="RK4", color="blue", alpha=0.7)
plt.plot(time_FE_1, z_FE_1, label="FE", color="green", alpha=0.5)
plt.plot(time_analytical, z_analytical, label="Analytical", color="red", alpha=0.5)
plt.xlabel(r'Time [$\mu$s]', fontsize=16)
plt.ylabel(r'$z$ [$\mu$m]', fontsize=16)
plt.legend(fontsize=16)
plt.grid(True)
plt.tight_layout()
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.savefig("8.1_z_t.pdf", format="pdf")
plt.show()

# ---------- Plotting - xy-plane without interactions ----------
plt.figure(figsize=(10, 6))

plt.plot(x_RK4_1, y_RK4_1, label="Particle 1 (RK4)", color="blue")
plt.plot(x_RK4_2, y_RK4_2, label="Particle 2 (RK4)", color="red")

plt.scatter(x_FE_1, y_FE_1, label="Particle 1 (FE)", color="lightsteelblue")
plt.scatter(x_FE_2, y_FE_2, label="Particle 2 (FE)", color="plum")

plt.xlabel(r'$x$ [$\mu$m]', fontsize=16)
plt.ylabel(r'$y$ [$\mu$m]', fontsize=16)
plt.legend(fontsize=16, loc='lower right')
plt.grid(True)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.tight_layout()
plt.savefig("XY-Plane_without_interactions.pdf", format="pdf")
plt.show()

# ---------- Plotting - Single plot of analytical solution of xy-plane for particle 1 ----------

plt.figure(figsize=(10, 6))

plt.plot(x_analytical, y_analytical, label = "Particle 1", color="red")

plt.xlabel(r'$x$ [$\mu$m]', fontsize=16)
plt.ylabel(r'$y$ [$\mu$m]', fontsize=16)
plt.legend(fontsize=16, loc='lower right')
plt.grid(True)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.tight_layout()
plt.savefig("XY-Plane_analytical.pdf", format="pdf")
plt.show()

# ---------- Plotting - xy-plane with interactions ----------
plt.figure(figsize=(10, 6))

plt.plot(x_RK4_interactions_1, y_RK4_interactions_1, label="Particle 1 (RK4)", color="blue")
plt.plot(x_RK4_interactions_2, y_RK4_interactions_2, label="Particle 2 (RK4)", color="red")

plt.xlabel(r'$x$ [$\mu$m]', fontsize=16)
plt.ylabel(r'$y$ [$\mu$m]', fontsize=16)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.legend(fontsize=16)
plt.grid(True)
plt.tight_layout()
plt.savefig("XY-Plane_RK4_interactions.pdf", format="pdf")
plt.show()

plt.figure(figsize=(10, 6))
plt.plot(x_FE_1_interactions, y_FE_1_interactions, label="Particle 1 (FE)", color="blue")
plt.plot(x_FE_2_interactions, y_FE_2_interactions, label="Particle 2 (FE)", color="red")

plt.xlabel(r'$x$ [$\mu$m]', fontsize=16)
plt.ylabel(r'$y$ [$\mu$m]', fontsize=16)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.legend(fontsize=16)
plt.grid(True)
plt.tight_layout()
plt.savefig("XY-Plane_FE_interactions.pdf", format="pdf")
plt.show()


# ---------- Plotting - trajectory phase space plot without interactions (x, v_x) ----------
plt.figure(figsize=(10, 6))
plt.plot(x_RK4_1, x_RK4_1_velocity, label="Particle 1", color="red")
plt.plot(x_RK4_2, x_RK4_2_velocity, label="Particle 2", color="blue")
plt.xlabel(r'$x$ [$\mu$m]', fontsize=16)
plt.ylabel(r'$v_x$ [m/s]', fontsize=16)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.legend(fontsize=16)
plt.grid(True)
plt.tight_layout()
plt.savefig("x_vx.pdf", format="pdf")
plt.show()

#---------- Plotting - trajectory phase space plot without interactions (z, v_z)----------
plt.figure(figsize=(10, 6))
plt.plot(z_RK4_1, z_RK4_1_velocity, label="Particle 1", color="red")
plt.plot(z_RK4_2, z_RK4_2_velocity, label="Particle 2", color="blue")
plt.xlabel(r'$z$ [$\mu$m]', fontsize=16)
plt.ylabel(r'$v_z$ [m/s]', fontsize=16)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.legend(fontsize=16)
plt.grid(True)
plt.tight_layout()
plt.savefig("z_vz.pdf", format="pdf")
plt.show()

# ---------- Plotting - trajectory phase space plot with interactions (x, v_x) ----------
plt.figure(figsize=(10, 6))
plt.plot(x_RK4_interactions_1, x_RK4_interactions_1_velocity, label="Particle 1", color="red")
plt.plot(x_RK4_interactions_2, x_RK4_interactions_2_velocity, label="Particle 2", color="blue")
plt.xlabel(r'$x$ [$\mu$m]', fontsize=16)
plt.ylabel(r'$v_x$ [m/s]', fontsize=16)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.legend(fontsize=16)
plt.grid(True)
plt.tight_layout()
plt.savefig("x_vx_interactions.pdf", format="pdf")
plt.show()

#---------- Plotting - trajectory phase space plot with interactions (z, v_z)----------
plt.figure(figsize=(10, 6))
plt.plot(z_RK4_interactions_1, z_RK4_interactions_1_velocity, label="Particle 1", color="red")
plt.plot(z_RK4_interactions_2, z_RK4_interactions_2_velocity, label="Particle 2", color="blue")
plt.xlabel(r'$z$ [$\mu$m]', fontsize=16)
plt.ylabel(r'$v_z$ [m/s]', fontsize=16)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.legend(fontsize=16)
plt.grid(True)
plt.tight_layout()
plt.savefig("z_vz_interactions.pdf", format="pdf")
plt.show()


# ---------- 3D Plotting - Without Interactions ----------
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection='3d')
ax.plot(x_RK4_1, y_RK4_1, z_RK4_1, label="Particle 1", color="red")
ax.plot(x_RK4_2, y_RK4_2, z_RK4_2, label="Particle 2", color="blue")
ax.set_xlabel(r'$x$ [$\mu$m]', fontsize=16)
ax.set_ylabel(r'$y$ [$\mu$m]', fontsize=16)
ax.set_zlabel(r'$z$ [$\mu$m]', fontsize=16)
ax.tick_params(axis='x', labelsize=14)
ax.tick_params(axis='y', labelsize=14)
ax.tick_params(axis='z', labelsize=14)
ax.legend(fontsize=16)
plt.tight_layout()
plt.savefig("3D.pdf", format="pdf")
plt.show()

# ---------- 3D Plotting - With Interactions ----------
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection='3d')
ax.plot(x_RK4_interactions_1, y_RK4_interactions_1, z_RK4_interactions_1, label="Particle 1", color="red")
ax.plot(x_RK4_interactions_2, y_RK4_interactions_2, z_RK4_interactions_2, label="Particle 2", color="blue")
ax.set_xlabel(r'$x$ [$\mu$m]', fontsize=16)
ax.set_ylabel(r'$y$ [$\mu$m]', fontsize=16)
ax.set_zlabel(r'$z$ [$\mu$m]', fontsize=16)
ax.tick_params(axis='x', labelsize=14)
ax.tick_params(axis='y', labelsize=14)
ax.tick_params(axis='z', labelsize=14)
ax.legend(fontsize=16)
plt.tight_layout()
plt.savefig("3D_interactions.pdf", format="pdf")
plt.show()


fig, axs = plt.subplots(1, 2, figsize=(12, 6))

# RK4 Plot
axs[0].plot(time_4000, np.log10(error_RK4_4000), label="4000 steps")
axs[0].plot(time_8000, np.log10(error_RK4_8000), label="8000 steps")
axs[0].plot(time_16000, np.log10(error_RK4_16000), label="16000 steps")
axs[0].plot(time_32000, np.log10(error_RK4_32000), label="32000 steps")

axs[0].set_xlabel(r'Time [$\mu$s]', fontsize=16)
axs[0].set_ylabel(r'log Error', fontsize=16)
axs[0].tick_params(axis='x', labelsize=14)
axs[0].tick_params(axis='y', labelsize=14)
axs[0].legend(loc="lower right", fontsize=14)
axs[0].set_title("Runge Kutta 4", fontsize=16)
axs[0].grid()

# FE Plot
axs[1].plot(time_4000, np.log10(error_FE_4000), label="4000 steps")
axs[1].plot(time_8000, np.log10(error_FE_8000), label="8000 steps")
axs[1].plot(time_16000, np.log10(error_FE_16000), label="16000 steps")
axs[1].plot(time_32000, np.log10(error_FE_32000), label="32000 steps")
axs[1].set_xlabel(r'Time [$\mu$s]', fontsize=16)
axs[1].tick_params(axis='x', labelsize=14)
axs[1].tick_params(axis='y', labelsize=14)
axs[1].set_title("Forward Euler", fontsize=16)
axs[1].legend(loc="lower right", fontsize=14)
axs[1].grid()

plt.tight_layout()
fig.savefig("relative_error.pdf", format="pdf")
plt.show()

#---------- Error difference plot ----------
fig, ax = plt.subplots(2, 2, figsize=(12, 8))
diff_4000 = abs(error_RK4_4000 - error_FE_4000)
diff_8000 = abs(error_RK4_8000 - error_FE_8000)
diff_16000 = abs(error_RK4_16000 - error_FE_16000)
diff_32000 = abs(error_RK4_32000 - error_FE_32000)

ax[0][0].plot(time_4000, diff_4000, label="steps=4000", color='blue')
ax[0][0].set_xlabel(r'Time [$\mu$s]', fontsize=16)
ax[0][0].set_ylabel(r'Error', fontsize=16)
ax[0][0].tick_params(axis='x', labelsize=16)
ax[0][0].tick_params(axis='y', labelsize=16)
ax[0][0].legend(fontsize=16)
ax[0][0].grid()

ax[0][1].plot(time_8000, diff_8000, label="steps=8000", color='orange')
ax[0][1].set_xlabel(r'Time [$\mu$s]', fontsize=16)
ax[0][1].set_ylabel(r'Error', fontsize=16)
ax[0][1].tick_params(axis='x', labelsize=16)
ax[0][1].tick_params(axis='y', labelsize=16)
ax[0][1].legend(fontsize=16)
ax[0][1].grid()

ax[1][0].plot(time_16000, diff_16000, label="steps=16000", color='green')
ax[1][0].set_xlabel(r'Time [$\mu$s]', fontsize=16)
ax[1][0].set_ylabel(r'Error', fontsize=16)
ax[1][0].tick_params(axis='x', labelsize=16)
ax[1][0].tick_params(axis='y', labelsize=16)
ax[1][0].legend(fontsize=16)
ax[1][0].grid()

ax[1][1].plot(time_32000, diff_32000, label="steps=32000", color='red')
ax[1][1].set_xlabel(r'Time [$\mu$s]', fontsize=16)
ax[1][1].set_ylabel(r'Error', fontsize=16)
ax[1][1].tick_params(axis='x', labelsize=16)
ax[1][1].tick_params(axis='y', labelsize=16)
ax[1][1].legend(fontsize=16)
ax[1][1].grid()

plt.tight_layout()
fig.savefig("difference_error_plot.pdf", format="pdf")
plt.show()
#---------- plot of particles left in trap ----------
# First plot
plt.figure()
plt.plot(omega_v_01, particles_01, label=r'f = 0.1 [$\mu$m]')
plt.xlabel(r"$\omega_v$ [Mhz]", fontsize=16)
plt.ylabel(r'Particles left in trap', fontsize=16)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.legend(fontsize=16)
plt.grid()
plt.tight_layout()
plt.savefig("Particle_f_01.pdf", format="pdf")
plt.show()

# Second plot
plt.figure()
plt.plot(omega_v_04, particles_04, label=r'f = 0.4 [$\mu$m]')
plt.xlabel(r"$\omega_v$ [Mhz]", fontsize=16)
plt.ylabel(r'Particles left in trap', fontsize=16)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.legend(fontsize=16)
plt.grid()
plt.tight_layout()
plt.savefig("Particle_f_04.pdf", format="pdf")
plt.show()

# Third plot
plt.figure()
plt.plot(omega_v_07, particles_07, label=r'f = 0.7 [$\mu$m]')
plt.xlabel(r"$\omega_v$ [Mhz]", fontsize=16)
plt.ylabel(r'Particles left in trap', fontsize=16)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.legend(fontsize=16)
plt.grid()
plt.tight_layout()
plt.savefig("Particle_f_07.pdf", format="pdf")
plt.show()


# Fourth plot with both for f=0.7
plt.figure()
plt.plot(omega_v_07_interaction, particles_07_interaction, label=r'Interactions')
plt.plot(omega_v_07_non, particles_07_non, label=r'No interactions')
plt.xlabel(r"$\omega_v$ [Mhz]", fontsize=16)
plt.ylabel(r'Particles left in trap', fontsize=16)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.legend(fontsize=16)
plt.grid()
plt.tight_layout()
plt.savefig("Particle_f_07_both.pdf", format="pdf")
plt.show()

# Fourth plot with both for f=0.1
plt.figure()
plt.plot(omega_v_01_interaction, particles_01_interaction, label=r'Interactions')
plt.plot(omega_v_01_non, particles_01_non, label=r'No interactions')
plt.xlabel(r"$\omega_v$ [Mhz]", fontsize=16)
plt.ylabel(r'Particles left in trap', fontsize=16)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.legend(fontsize=16)
plt.grid()
plt.tight_layout()
plt.savefig("Particle_f_01_both.pdf", format="pdf")
plt.show()
