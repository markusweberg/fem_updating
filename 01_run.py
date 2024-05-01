# -- coding utf-8 --

"""
01_run.py

File for performing FEM updating. Results from the FE model updating are
stored in the file info.pkl.
"""

import os
import shutil
import pickle
import time
import numpy as np

from scipy.optimize import lsq_linear
from collections import OrderedDict
from functions_py3 import *


# Total time - start
total_time_start = time.time()

# -------------------------------------------------------------------------- #
# ANALYSIS SETTINGS
# -------------------------------------------------------------------------- #

# Number of iterations to perform
iterations = 50

# Measured frequencies to perform FE model updating on

mode_list = []

for i in range(0, 17):
    mode_list.append('Mode ' + str('{:02}'.format(i+1)))

# -------------------------------------------------------------------------- #
# INITIAL
# -------------------------------------------------------------------------- #

# Folder structure
folder_1 = '02_Analysis files'
folder_2 = '03_Results/01_Perturbed'
folder_3 = '04_Figures'

folders = [folder_1, folder_2, folder_3]

# Establish folder structure
for folder in folders:

    # Remove existing folders with content
    if os.path.isdir(folder):
        shutil.rmtree(folder)
        time.sleep(1)

    # Create folder
    os.makedirs(folder)

# Import measured frequencies and mode shapes
freqs_est = np.load('01_Initial/02_Measured/grenland_bridge_frequencies.npy')
modes_est = np.load('01_Initial/02_Measured/grenland_bridge_modes.npy')

# Copy relevant model files to directory
src = ('01_Initial/01_FE model/')

files_to_copy = ['grenland_bridge.inp', '1-Skrastag.inp', '2-Viadukt.inp']

# Copy
for file in files_to_copy:

    shutil.copy(src + file, file)

# Set work directory
cwd = os.getcwd()
os.chdir(cwd)

# ---------------------------------------------------------------------------#
# ABAQUS MODEL
# ---------------------------------------------------------------------------#

# ------------------------------------ #
# FEM UPDATING PARAMETERS
# ------------------------------------ #

# FEM initial updating parameters
steel_e_modulus = 2.1E+11
steel_density = 7850
concrete_e_modulus = 0.27E+11
concrete_density = 2500
joint_stiffness_y = 1E+10
joint_stiffness_z = 1E+10
joint_stiffness_phix = 1E+12
joint_stiffness_phiy = 1E+12
joint_stiffness_phiz = 1E+12

# Parameter list
theta_0 = np.array([steel_e_modulus, steel_density, concrete_e_modulus, concrete_density,
                    joint_stiffness_y, joint_stiffness_z, joint_stiffness_phix, joint_stiffness_phiy, 
                    joint_stiffness_phiz])

# ------------------------------------ #
# FEM ANALYSIS
# ------------------------------------ #

# ---------------- #
# EXPORT VARIABLES
# ---------------- #

# Job name
job_name = 'grenland_bridge'
with open('job_name_list' + '.txt', 'w') as fout:
    fout.write(job_name)

# Parameter list
np.save('parameter_list' + '.npy', theta_0, allow_pickle=True,
        fix_imports=True)

# ---------------- #
# RUN ANALYSIS
# ---------------- #

# Variables
script_name = 'abaqus_upd'

# Run script
print('--------------------------------------------')
print('Initial analysis (1)')
t0 = time.time()

os.system('abaqus cae noGUI=' + script_name)

t1 = time.time()
print('Initial analysis (1) done - ' + str(round(t1 - t0, 3)) + ' sec.')
print('--------------------------------------------')
print()

# Wait 3 sec for the Abaqus analysis to properly close
time.sleep(3)

# ------------------------------------ #
# POST-PROCESS
# ------------------------------------ #

# Import
import_folder = '03_Results/'

freqs_num = np.load(import_folder + job_name + '_frequencies_all.npy')
modes_num = np.load(import_folder + job_name + '_modes_all.npy')

# Sorted MAC and MMI
(MAC_initial_model, MMI_initial_model, results_MAC, results_MMI) = get_MAC_MMI(freqs_est, modes_est,
                                                     freqs_num, modes_num,
                                                     filtering=False)

# General parameters
rows_est, cols_est = np.shape(modes_est)
rows_num, cols_num = np.shape(modes_num)


# ---------------------------------------------------------------------------#
# FEM UPDATING
# ---------------------------------------------------------------------------#

# ------------------------------------ #
# INITIALIZATION
# ------------------------------------ #

# Variables
q = len(mode_list) * 2
p = len(theta_0)

# Declaration of variables
# Variables saved for each iteration
theta = np.zeros((iterations + 1, p), dtype=float)
lamda = np.zeros((iterations + 1, q), dtype=float)
lamda_all = np.zeros((iterations + 1, len(freqs_est)), dtype=float)
lamda_m = np.zeros(q, dtype=float)
lamda_start = np.zeros(q, dtype=float)

G_matrices = OrderedDict([])
MAC_updated_model = OrderedDict([])
MMI_updated_model = OrderedDict([])
results_MAC_updated_model = OrderedDict([])
results_MAC_pt_all_matrices = OrderedDict([])
results_MMI_updated_model = OrderedDict([])
results_MMI_pt_all_matrices = OrderedDict([])

J_star_list = []

# Weighting matrix
# W = np.linalg.inv(np.diag(lamda_m))**2
# W = np.identity(len(lamda_m))
W = np.zeros((len(lamda_m), len(lamda_m)))

# Start values for the iterations
for idx, mode in enumerate(mode_list):

    # Establish frequencies based on the mode_list string
    string_idx = mode[-2:]

    # Correct for 0 in the string
    if string_idx[0] == 0:

        f_idx = int(mode[-1]) - 1

    else:
        f_idx = int(string_idx) - 1

    # Measured frequencies (to perform FEM updating on)
    lamda_m[idx] = freqs_est[:, 1][f_idx]

    # Numerical frequencies (to perform FEM updating on)
    # Top match based on MMI
    lamda[0, idx] = results_MMI[mode][0, 2]

    # MAC value (to perform FEM updating on)
    lamda[0, idx + len(mode_list)] = results_MMI[mode][0, 1]

    # Target MAC value
    lamda_m[idx + len(mode_list)] = 1

    # Start values for lamda
    lamda_start[idx] = results_MMI[mode][0, 2]
    lamda_start[idx + len(mode_list)] = 1

    # Weighting matrix
    if idx < 4:
        W[idx, idx] = (2/3) / len(mode_list) + (1/40) * ((len(mode_list)-4)/len(mode_list))
        W[idx + len(mode_list), idx + len(mode_list)] = (1/3) / len(mode_list) + (1/40) * ((len(mode_list)-4)/len(mode_list))
    else:
        W[idx, idx] = (2/3) / len(mode_list) - (1/40) * (4/len(mode_list))
        W[idx + len(mode_list), idx + len(mode_list)] = (1/3) / len(mode_list) - (1/40) * (4/len(mode_list))

# All numerical frequencies
for idx, mode in enumerate(results_MMI.keys()):

    lamda_all[0, idx] = results_MMI[mode][0, 2]

# Updating parameters
theta[0, :] = np.copy(theta_0)
d_theta_start = 0.001*theta_0

# ------------------------------------ #
# OBJECTIVE FUNCTION
# ------------------------------------ #

# Evaluate the objective function for the initial analysis
J_star = np.sum(np.diag(W)*((lamda_m - lamda[0, :])/lamda_start)**2)
J_star_list.append(J_star)


# ------------------------------------ #
# ITERATIONS
# ------------------------------------ #

for itr in range(iterations):

    # First iteration step
    if itr == 0:

        # Initialize the parameter incremenet
        d_theta = d_theta_start

    # All other iteration steps
    else:

        # Update the parameter increment
        # d_theta = d_theta_upd
        d_theta = d_theta_start

    # Variable
    results_MAC_pt_matrices = OrderedDict([])
    results_MMI_pt_matrices = OrderedDict([])

    # -------------------------- #
    # SENSITIVITY MATRIX
    # -------------------------- #

    # Sensitivity matrix (scaled)
    G = np.zeros((q, p), dtype=float)

    # Loop through the parameters (columns)
    for ii in range(p):

        # -------------------------- #
        # PERTURBATION - THETA
        # -------------------------- #

        # Perturbed parameter vector
        theta_pt = np.copy(theta[itr, :])
        theta_pt[ii] = theta_pt[ii] + d_theta[ii]

        # -------------------------- #
        # PERTURBATION - LAMDA
        # -------------------------- #

        # -------------------------- #
        # FE ANALYSIS
        # -------------------------- #

        # ---------------- #
        # EXPORT VARIABLES
        # ---------------- #

        # Job name
        job_name = 'grenland_bridge' + str(itr+1) + '_pt' + str(ii+1)
        with open('job_name_list' + '.txt', 'a') as fout:
            fout.write('\n' + job_name)

        # Parameter list
        np.save('parameter_list_pt' + '.npy', theta_pt, allow_pickle=True,
                fix_imports=True)

        # ---------------- #
        # RUN ANALYSIS
        # ---------------- #

        # Variables
        script_name = 'abaqus_upd_pt'

        # Run script
        print('Perturbation analysis ' + str(ii+1))
        t0_pt = time.time()

        os.system('abaqus cae noGUI=' + script_name)

        t1_pt = time.time()
        print('Perturbation analysis ' + str(ii+1) + ' done - ' +
              str(round(t1_pt - t0_pt, 3)) + ' sec.')
        print('--------------------------------------------')
        print()

        # Wait 3 sec for the Abaqus analysis to properly close
        time.sleep(3)

        # -------------------------- #
        # POST-PROCESS
        # -------------------------- #

        # Import
        import_folder = '03_Results/01_Perturbed/'

        freqs_num_pt = np.load(import_folder + job_name +
                               '_frequencies_all.npy')

        modes_num_pt = np.load(import_folder + job_name + '_modes_all.npy')

        # Sorted MAC and MMI
        (_, _, results_MAC_pt, results_MMI_pt) = get_MAC_MMI(freqs_est, modes_est,
                                                freqs_num_pt, modes_num_pt,
                                                filtering=False)

        # Append perturbed results for validation
        results_MAC_pt_matrices.update({'results_MAC_pt' + str(ii+1) : 
                                        results_MAC_pt})
        
        results_MMI_pt_matrices.update({'results_MMI_pt' + str(ii+1) : 
                                        results_MMI_pt})        

        # Perturbed numerical frequencies
        lamda_pt = np.zeros(q, dtype=float)
        lamda_pt_all = np.zeros(cols_est, dtype=float)

        for idx, mode in enumerate(mode_list):

            # Establish frequencies based on the mode_list string
            string_idx = mode[-2:]

            # Correct for 0 in the string
            if string_idx[0] == 0:

                f_idx = int(mode[-1]) - 1

            else:
                f_idx = int(string_idx) - 1

            # Perturbed numerical frequencies
            # Top match based on MMI
            lamda_pt[idx] = results_MMI_pt[mode][0, 2]
            
            # MAC value (to perform FEM updating on)
            lamda_pt[idx + len(mode_list)] = results_MMI_pt[mode][0, 1]

        # All perturbed numerical frequencies
        for idx, mode in enumerate(results_MMI_pt.keys()):

            lamda_pt_all[idx] = results_MMI_pt[mode][0, 2]

        # Sensitivities
        # Loop through the frequencies (rows)
        for jj in range(q):

            # Numerator
            nmr = (lamda_pt[jj] - lamda[itr, jj]) / lamda_start[jj]

            # Denominator
            dnm = (theta_pt[ii] - theta[itr, ii]) / theta[0, ii]

            # Update sensitivity matrix
            G[jj, ii] = nmr/dnm

    # Append
    G_matrices.update({'G' + str(itr+1): G})

    # Append
    results_MAC_pt_all_matrices.update({'results_MAC_pt_matrices_' +
                                        str(itr+1): results_MAC_pt_matrices})
    
    results_MMI_pt_all_matrices.update({'results_MMI_pt_matrices_' +
                                        str(itr+1): results_MMI_pt_matrices})

    # -------------------------- #
    # OPTIMIZATION
    # -------------------------- #

    # Residual, r, scaled
    r_s = (lamda_m - lamda[itr, :]) / lamda_start

    # Global parameter bounds
    lower_allowable_theta = np.array([1.7E+11, 7000, 0.15E+11, 2000, 0, 0, 0, 0, 0])
    upper_allowable_theta = np.array([2.7E+11, 9000, 0.35E+11, 3000, 1E+14, 1E+14, 1E+16, 1E+16, 1E+16])

    # Local parameter bounds
    lb_local = np.array([0.90, 0.98, 0.90, 0.98, 0.80, 0.80, 0.80, 0.80, 0.80])
    ub_local = np.array([1.10, 1.02, 1.10, 1.02, 1.20, 1.20, 1.20, 1.20, 1.20])

    # Local parameter bounds, unscaled
    theta_min = lb_local*theta[itr, :]
    theta_max = ub_local*theta[itr, :]

    # Check for exceedance of global bounds
    for index in range(p):

        if theta_min[index] <= lower_allowable_theta[index]:
            theta_min[index] = np.copy(lower_allowable_theta[index])

        if theta_max[index] >= upper_allowable_theta[index]:
            theta_max[index] = np.copy(upper_allowable_theta[index])

    # Upper and lower parameter bounds
    d_theta_min = theta_min - theta[itr, :]
    d_theta_max = theta_max - theta[itr, :]

    # Upper and lower parameter bounds, scaled
    d_theta_min_s = d_theta_min / theta[0, :]
    d_theta_max_s = d_theta_max / theta[0, :]

    # Minimization variables
    A = -np.sqrt(W) @ G
    b = -np.sqrt(W) @ r_s
    lb = np.copy(d_theta_min_s)
    ub = np.copy(d_theta_max_s)

    # Minimization algorithm
    res = lsq_linear(A, b, bounds=(lb, ub), method='trf',
                     lsq_solver=None, lsmr_tol='auto', verbose=1)

    # -------------------------- #
    # UPDATE
    # -------------------------- #

    # Parameter incremental, d_theta, scaled
    d_theta_s_upd = np.copy(res.x)

    # Update d_theta_upd (unscaled)
    d_theta_upd = d_theta_s_upd*theta[0, :]

    # Update theta
    theta[itr+1, :] = theta[itr, :] + d_theta_upd

    # ------------------------------------ #
    # FEM ANALYSIS
    # ------------------------------------ #

    # ---------------- #
    # EXPORT VARIABLES
    # ---------------- #

    # Job name
    job_name = 'grenland_bridge' + str(itr+2)
    with open('job_name_list' + '.txt', 'a') as fout:
        fout.write('\n' + job_name)

    # Parameter list
    np.save('parameter_list_upd' + '.npy', theta[itr+1, :], allow_pickle=True,
            fix_imports=True)

    # ---------------- #
    # RUN ANALYSIS
    # ---------------- #

    # Variables
    script_name = 'abaqus_upd'

    # Run script
    print()
    print('--------------------------------------------')
    print('Updated analysis (' + str(itr+2) + ')')
    t0 = time.time()

    os.system('abaqus cae noGUI=' + script_name)

    t1 = time.time()
    print('Updated analysis (' + str(itr+2) + ') done - ' +
          str(round(t1 - t0, 3)) + ' sec.')
    print('--------------------------------------------')
    print()

    # Wait 3 sec for the Abaqus analysis to properly close
    time.sleep(3)

    # ------------------------------------ #
    # POST-PROCESS
    # ------------------------------------ #

    # Import
    import_folder = '03_Results/'

    freqs_num_upd = np.load(import_folder + job_name + '_frequencies_all.npy')
    modes_num_upd = np.load(import_folder + job_name + '_modes_all.npy')

    # Sorted MAC and MMI
    (MAC_updated, MMI_updated, results_MAC_upd, results_MMI_upd) = get_MAC_MMI(freqs_est,
                                                       modes_est,
                                                       freqs_num_upd,
                                                       modes_num_upd,
                                                       filtering=False)

    # Updated frequencies
    lamda_upd = np.zeros(q, dtype=float)
    lamda_upd_all = np.zeros(cols_est, dtype=float)

    for idx, mode in enumerate(mode_list):

        # Establish frequencies based on the mode_list string
        string_idx = mode[-2:]

        # Correct for 0 in the string
        if string_idx[0] == 0:

            f_idx = int(mode[-1]) - 1

        else:
            f_idx = int(string_idx) - 1

        # Updated numerical frequencies
        # Top match based on MMI
        lamda_upd[idx] = results_MMI_upd[mode][0, 2]

        # Updated MAC value
        lamda_upd[idx + len(mode_list)] = results_MMI_upd[mode][0, 1]

    # All updated numerical frequencies
    for idx, mode in enumerate(results_MMI_upd.keys()):

        lamda_upd_all[idx] = results_MMI_upd[mode][0, 2]

    # Update lamda
    lamda[itr+1, :] = np.copy(lamda_upd)
    lamda_all[itr+1, :] = np.copy(lamda_upd_all)

    # Append updated MAC and MMI results
    MAC_updated_model.update({'MAC_updated' + str(itr+1): MAC_updated})
    results_MAC_updated_model.update({'results_MAC_upd' + str(itr+1):
                                      results_MAC_upd})
    
    MMI_updated_model.update({'MMI_updated' + str(itr+1): MMI_updated})
    results_MMI_updated_model.update({'results_MMI_upd' + str(itr+1):
                                      results_MMI_upd})

    # -------------------------- #
    # OBJECTIVE FUNCTION
    # -------------------------- #

    # Objective function for evaluation
    J_star = np.sum(np.diag(W)*((lamda_m - lamda[itr+1, :])/lamda_start)**2)
    J_star_list.append(J_star)


# Wait 5 sec for the Abaqus analysis to properly close
print('Sleeping 5 sec for Abaqus to properly close.')
time.sleep(5)


# -------------------------------------------------------------------------- #
# MOVE FILES
# -------------------------------------------------------------------------- #

# Move analysis files to folder
dst = '02_Analysis files/'
src_file_list = os.listdir()

# Analysis files
for src_file in src_file_list:

    # Move Abaqus analysis files
    if src_file.startswith('grenland_bridge') or src_file.startswith('abaqus.rpy'):

        shutil.move(src_file, dst)

    elif src_file.startswith("1-Skrastag") or src_file.startswith("2-Viadukt"):

        shutil.move(src_file, dst)

    # Move Abaqus temporary files
    elif src_file.startswith('temp'):

        shutil.move(src_file, dst)

    # Move other files
    elif src_file.endswith('.npy'):

        shutil.move(src_file, dst)


# -------------------------------------------------------------------------- #
# EXPORT
# -------------------------------------------------------------------------- #

# Export relevant variables to postprocessing
info = {'freqs_est': freqs_est,
        'modes_est': modes_est,
        'freqs_num': freqs_num,
        'modes_num': modes_num,
        'theta': theta,
        'MAC_initial': MAC_initial_model,
        'MMI_initial': MMI_initial_model,
        'MAC_updated': MAC_updated_model,
        'MMI_updated': MMI_updated_model,
        'lamda_m': lamda_m,
        'lamda': lamda,
        'lamda_all': lamda_all,
        'G_matrices': G_matrices,
        'W': W,
        'J_star_list': J_star_list}

with open('info.pkl', 'wb') as file:
    pickle.dump(info, file)


# -------------------------------------------------------------------------- #
# END OF ANALYSIS
# -------------------------------------------------------------------------- #

# Total time - end
total_time_stop = time.time()
print()
print('Total analysis time is ' +
      str(round((total_time_stop - total_time_start), 1)) + ' sec.')
