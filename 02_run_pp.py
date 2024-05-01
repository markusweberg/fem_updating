# -- coding utf-8 --

"""
02_run_pp.py

Postprocessing of results from the FEM updating.
"""


import numpy as np
import matplotlib.pyplot as plt
import pickle


# ---------------------------------------------------------------------------#
# DATA
# ---------------------------------------------------------------------------#

# Import data
with open('info.pkl', 'rb') as info_file:
    info = pickle.load(info_file)


# ---------------------------------------------------------------------------#
# PLOT
# ---------------------------------------------------------------------------#

# ------------------------------------ #
# SETTINGS
# ------------------------------------ #

fig_format1 = '.png'
fig_format2 = '.svg'
fig_folder = '04_Figures/'

# Matplotlib customization
plt.rcParams.update({'mathtext.fontset': 'stix'})
plt.rcParams.update({'font.size': 8})
plt.rcParams.update({'lines.linewidth': 0.75})


# ------------------------------------ #
# OBJECTIVE FUNCTION
# ------------------------------------ #

# Plot
iterations = np.arange(0, len(info['J_star_list']))

plt.figure(1, figsize=(10, 8))
plt.title('Objective function', fontweight='bold')
plt.plot(iterations, info['J_star_list'], color='black', linewidth=0.75)
plt.xticks(iterations)
plt.xlabel('Iterations')
plt.ylabel(r'$J(\Delta\theta)$')
plt.tight_layout()

# Save
plt.savefig(fig_folder + 'fig1' + fig_format1)
plt.savefig(fig_folder + 'fig1' + fig_format2)


# ------------------------------------ #
# SENSITIVITY MATRIX
# ------------------------------------ #

# Normalize
G = np.abs(info['G_matrices']['G50'])
# G = G / np.amax(G)
rows, cols = np.shape(G)
G[0:int(rows/2),:] = G[0:int(rows/2),:] / np.amax(G[0:int(rows/2),:])
G[int(rows/2):,:] = G[int(rows/2):,:] / np.amax(G[int(rows/2):,:])

# 3D plot
fig = plt.figure(2, figsize=(10, 8))
ax = plt.axes(projection='3d')
ax.set_title('Sensitivity matrix (Natural frequencies)', fontweight='bold')

for j in range(int(rows/2)):

    # Centering of the bars
    dx = 0.25
    dy = 0.25

    # Coordinates of the bars
    x = []
    y = []
    z = np.zeros(len(info["theta"][0]))

    # Size of the bars 
    dx_vec = []
    dy_vec = []
    dz_vec = G[j, :]

    for i in range(len(info["theta"][0])):
        x.append(i - dx)
        y.append(j - dy)

        dx_vec.append(0.5)
        dy_vec.append(0.5)
    
    # Colormap
    cmap_alt1 = plt.cm.viridis(dz_vec)
    cmap_alt2 = plt.cm.rainbow(dz_vec)

    # Plot
    bar1 = ax.bar3d(x, y, z, dx_vec, dy_vec, dz_vec, color=cmap_alt1,
                    zsort='max', alpha=1.0)

ax.set_xticks([i for i in range(len(info["theta"][0]))])
ax.set_xticklabels([r'$E_s$', r'$g_s$', r'$E_c$', r'$g_c$', r'$k_y$', r'$k_z$', 
                    r'$k_{\phi_x}$', r'$k_{\phi_y}$', r'$k_{\phi_z}$'])
ax.set_yticks([i for i in range(int(rows/2))])
ax.set_yticklabels([f"$f_\u007b{i+1}\u007d$" for i in range(int(rows/2))])
ax.set_zticks([0, 0.25, 0.5, 0.75, 1.0])

# Save
plt.savefig(fig_folder + 'fig2' + fig_format1)
plt.savefig(fig_folder + 'fig2' + fig_format2)

# 3D plot
fig = plt.figure(3, figsize=(10, 8))
ax = plt.axes(projection='3d')
ax.set_title('Sensitivity matrix (Mode shapes)', fontweight='bold')

for j in range(int(rows/2)):

    # Centering of the bars
    dx = 0.25
    dy = 0.25

    # Coordinates of the bars
    x = []
    y = []
    z = np.zeros(len(info["theta"][0]))

    # Size of the bars 
    dx_vec = []
    dy_vec = []
    dz_vec = G[j + int(rows/2), :]

    for i in range(len(info["theta"][0])):
        x.append(i - dx)
        y.append(j - dy)

        dx_vec.append(0.5)
        dy_vec.append(0.5)
    
    # Colormap
    cmap_alt1 = plt.cm.viridis(dz_vec)
    cmap_alt2 = plt.cm.rainbow(dz_vec)

    # Plot
    bar1 = ax.bar3d(x, y, z, dx_vec, dy_vec, dz_vec, color=cmap_alt1,
                    zsort='max', alpha=1.0)

ax.set_xticks([i for i in range(len(info["theta"][0]))])
ax.set_xticklabels([r'$E_s$', r'$g_s$', r'$E_c$', r'$g_c$', r'$k_y$', r'$k_z$', 
                    r'$k_{\phi}$', r'$k_{\phi}$', r'$k_{\phi}$'])
ax.set_yticks([i for i in range(int(rows/2))])
ax.set_yticklabels([f"$MAC_\u007b{i+1}\u007d$" for i in range(int(rows/2))])
ax.set_zticks([0, 0.25, 0.5, 0.75, 1.0])

# Save
plt.savefig(fig_folder + 'fig3' + fig_format1)
plt.savefig(fig_folder + 'fig3' + fig_format2)

# ------------------------------------ #
# MAC
# ------------------------------------ #

MAC_initial = info['MAC_initial']
MAC_updated = info['MAC_updated']['MAC_updated50']
MAC_updated = np.round(MAC_updated, 2)
rows = np.shape(MAC_initial)[1]
cols = np.shape(MAC_initial)[1]

# Plot parameters
x_tick_list = np.arange(0, cols, 1)
y_tick_list = np.arange(0, rows, 1)
naming_cols = []
naming_rows_initial = []
naming_rows_updated = []
ind_list_initial = []
ind_list_updated = []

for k in range(cols):
    name = str(k+1)
    naming_cols.append(name)

# Sorting initial MAC matrix
for i in range(len(MAC_initial[1])):
    ind = np.argsort(MAC_initial[:,i])[::-1]
    ind_list_initial.append(ind[0])
    naming_rows_initial.append(ind[0] + 1)

MAC_initial_sorted = MAC_initial[ind_list_initial,:]

# Sorting updated MAC matrix
for i in range(len(MAC_updated[1])):
    ind = np.argsort(MAC_updated[:,i])[::-1]
    ind_list_updated.append(ind[0])
    naming_rows_updated.append(ind[0] + 1)

MAC_updated_sorted = MAC_updated[ind_list_updated,:]

# Initial MAC plot
plt.figure(4, figsize=(10, 8))
plt.imshow(MAC_initial_sorted, cmap='viridis', interpolation=None, vmin=0, vmax=1,
           origin='lower')
plt.title('Initial MAC', fontweight='bold')
plt.colorbar()
plt.xlabel('Identified modes', fontweight='bold')
plt.ylabel('Numerical modes', fontweight='bold')
plt.xticks(x_tick_list, naming_cols)
plt.yticks(y_tick_list, naming_rows_initial,
           rotation=90, va='center')

for i in range(len(MAC_initial_sorted[0])):
    for j in range(len(MAC_initial_sorted[1])):
        plt.text(j, i, round(MAC_initial_sorted[i, j], 2), ha="center", va="center",
                 color="white", fontsize=6)

plt.tight_layout()

# Save
plt.savefig(fig_folder + 'fig4' + fig_format1)
plt.savefig(fig_folder + 'fig4' + fig_format2)

# Updated MAC plot
plt.figure(5, figsize=(10, 8))
plt.imshow(MAC_updated_sorted, cmap='viridis', interpolation=None, vmin=0, vmax=1,
           origin='lower')
plt.title('Updated MAC', fontweight='bold')
plt.colorbar()
plt.xlabel('Identified modes', fontweight='bold')
plt.ylabel('Numerical modes', fontweight='bold')
plt.xticks(x_tick_list, naming_cols)
plt.yticks(y_tick_list, naming_rows_updated,
           rotation=90, va='center')

for i in range(len(MAC_updated_sorted[0])):
    for j in range(len(MAC_updated_sorted[1])):
        plt.text(j, i, round(MAC_updated_sorted[i, j], 2), ha="center", va="center",
                 color="white", fontsize=6)

plt.tight_layout()

# Save
plt.savefig(fig_folder + 'fig5' + fig_format1)
plt.savefig(fig_folder + 'fig5' + fig_format2)

# ------------------------------------ #
# MMI
# ------------------------------------ #

MMI_initial = info['MMI_initial']
MMI_updated = info['MMI_updated']['MMI_updated50']
MMI_updated = np.round(MMI_updated, 2)
rows = np.shape(MMI_initial)[1]
cols = np.shape(MMI_initial)[1]

# Plot parameters
x_tick_list = np.arange(0, cols, 1)
y_tick_list = np.arange(0, rows, 1)
naming_cols = []
naming_rows_initial = []
naming_rows_updated = []
ind_list_initial = []
ind_list_updated = []

for k in range(cols):
    name = str(k+1)
    naming_cols.append(name)

# Sorting initial MMI matrix
for i in range(len(MMI_initial[1])):
    ind = np.argsort(MMI_initial[:,i])[::-1]
    ind_list_initial.append(ind[0])
    naming_rows_initial.append(ind[0] + 1)

MMI_initial_sorted = MMI_initial[ind_list_initial,:]

# Sorting updated MMI matrix
for i in range(len(MMI_updated[1])):
    ind = np.argsort(MMI_updated[:,i])[::-1]
    ind_list_updated.append(ind[0])
    naming_rows_updated.append(ind[0] + 1)

MMI_updated_sorted = MMI_updated[ind_list_updated,:]

# Initial MMI plot
plt.figure(6, figsize=(10, 8))
plt.imshow(MMI_initial_sorted, cmap='viridis', interpolation=None, vmin=0, vmax=0.5,
           origin='lower')
plt.title('Initial MMI', fontweight='bold')
plt.colorbar()
plt.xlabel('Identified modes', fontweight='bold')
plt.ylabel('Numerical modes', fontweight='bold')
plt.xticks(x_tick_list, naming_cols)
plt.yticks(y_tick_list, naming_rows_initial,
           rotation=90, va='center')

for i in range(len(MMI_initial_sorted[0])):
    for j in range(len(MMI_initial_sorted[1])):
        plt.text(j, i, round(MMI_initial_sorted[i, j], 2), ha="center", va="center",
                 color="white", fontsize=6)

plt.tight_layout()

# Save
plt.savefig(fig_folder + 'fig6' + fig_format1)
plt.savefig(fig_folder + 'fig6' + fig_format2)

# Updated MMI plot
plt.figure(7, figsize=(10, 8))
plt.imshow(MMI_updated_sorted, cmap='viridis', interpolation=None, vmin=0, vmax=0.5,
           origin='lower')
plt.title('Updated MMI', fontweight='bold')
plt.colorbar()
plt.xlabel('Identified modes', fontweight='bold')
plt.ylabel('Numerical modes', fontweight='bold')
plt.xticks(x_tick_list, naming_cols)
plt.yticks(y_tick_list, naming_rows_updated,
           rotation=90, va='center')

for i in range(len(MMI_updated_sorted[0])):
    for j in range(len(MMI_updated_sorted[1])):
        plt.text(j, i, round(MMI_updated_sorted[i, j], 2), ha="center", va="center",
                 color="white", fontsize=6)

plt.tight_layout()

# Save
plt.savefig(fig_folder + 'fig7' + fig_format1)
plt.savefig(fig_folder + 'fig7' + fig_format2)

plt.show()
