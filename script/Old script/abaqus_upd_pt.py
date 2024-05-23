# -- coding utf-8 --
"""
abaqus_upd_pt.py

Perturbation analysis of simply supported beam model with updating parameters.
"""

from abaqus import *
from abaqusConstants import *

from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *

import os
import shutil
import numpy as np
from functions_py2 import *


# -------------------------------------------------------------------------- #
# INITIAL
# -------------------------------------------------------------------------- #

# Job names
job_name_list = np.loadtxt('job_name_list' + '.txt', dtype=str)
jobName = str(job_name_list[-1])
jobName_ending = jobName[15:]

skrastag_jobName = "1-Skrastag" + jobName_ending
viadukt_jobName = "2-Viadukt" + jobName_ending

# Creating the input files
shutil.copy('grenland_bridge.inp', jobName + '.inp')
shutil.copy("1-Skrastag.inp", skrastag_jobName + ".inp")
shutil.copy("2-Viadukt.inp", viadukt_jobName + ".inp")

# Updating parameters
theta_pt = np.load('parameter_list_pt' + '.npy')


# ---------------------------------------------------------------------------#
# ABAQUS MODEL
# ---------------------------------------------------------------------------#


# ------------------------------------ #
# FEM UPDATING PARAMETERS
# ------------------------------------ #

# FEM initial updating parameters
steel_e_modulus = theta_pt[0]
steel_density = theta_pt[1]
tower_concrete_e_modulus = theta_pt[2]
cable_steel_e_modulus = theta_pt[3]
joint_stiffness_y = theta_pt[4]
joint_stiffness_z = theta_pt[5]
joint_stiffness_phix = theta_pt[6]
joint_stiffness_phiy = theta_pt[7]
joint_stiffness_phiz = theta_pt[8]

# ------------------------------------ #
# FEM ANALYSIS
# ------------------------------------ #

skrastag_bsp_file = open(skrastag_jobName + ".inp", 'r').readlines()

skrastag_bsp_file[30405] = str(steel_e_modulus) + ", 0.3\n"
skrastag_bsp_file[30403] = str(steel_density) + "\n"

skrastag_bsp_file_upd = open(skrastag_jobName + ".inp", "w")
skrastag_bsp_file_upd.writelines(skrastag_bsp_file)
skrastag_bsp_file_upd.close()

os.system("abaqus job=" + skrastag_jobName + " interactive")

viadukt_bsp_file  = open(viadukt_jobName + ".inp", "r").readlines()

viadukt_bsp_file[26284] = str(steel_e_modulus) + ", 0.3\n"
viadukt_bsp_file[26282] = str(steel_density) + "\n"

viadukt_bsp_file_upd = open(viadukt_jobName + ".inp", "w")
viadukt_bsp_file_upd.writelines(viadukt_bsp_file)
viadukt_bsp_file_upd.close()

os.system("abaqus job=" + viadukt_jobName + " interactive")

inp_file = open(jobName + '.inp', 'r').readlines()

inp_file[9] = str(cable_steel_e_modulus) + ", 0.3\n"
inp_file[15] = str(tower_concrete_e_modulus) + ", 0.2\n"
inp_file[3681] = str(joint_stiffness_y) + "\n"
inp_file[3683] = str(joint_stiffness_z) + "\n"
inp_file[3685] = str(joint_stiffness_phix) + "\n"
inp_file[3687] = str(joint_stiffness_phiy) + "\n"
inp_file[3689] = str(joint_stiffness_phiz) + "\n"

inp_file[2378] = "*INCLUDE, INPUT=" + skrastag_jobName + ".bsp\n"
inp_file[3616] = "*INCLUDE, INPUT=" + viadukt_jobName + ".bsp\n"

inp_file_upd = open(jobName + '.inp', 'w')
inp_file_upd.writelines(inp_file)
inp_file_upd.close()

os.system("abaqus job=" + jobName + " interactive")

# ------------------------------------ #
# FEM POSTPROCESSING
# ------------------------------------ #

# Get data from ODB file
o3 = session.openOdb(name=jobName + '.odb')
odb = session.odbs[jobName + '.odb']

# Frequencies
freqs_num_pt = get_frequencies(odb, save=True, name=jobName,
                               folder='03_Results/01_Perturbed/')

# Mode shapes
modes_num_pt = get_modeshapes(odb, node_set_name_1='ACCELEROMETERS_1', node_set_name_2='ACCELEROMETERS_2', save=True,
                              name=jobName, folder='03_Results/01_Perturbed/')
