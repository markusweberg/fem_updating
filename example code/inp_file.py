# Opening the ".inp" file and inserting the parameter at the correct line
inp_file = open(jobName + '.inp', 'r').readlines()

inp_file[9] = str(cable_steel_e_modulus) + ", 0.3\n"
inp_file[15] = str(tower_concrete_e_modulus) + ", 0.2\n"
inp_file[3681] = str(joint_stiffness_y) + "\n"
inp_file[3683] = str(joint_stiffness_z) + "\n"
inp_file[3685] = str(joint_stiffness_phix) + "\n"
inp_file[3687] = str(joint_stiffness_phiy) + "\n"
inp_file[3689] = str(joint_stiffness_phiz) + "\n"

# Updating the ".bsp" file name
inp_file[2378] = "*INCLUDE, INPUT=" + skrastag_jobName + ".bsp\n"
inp_file[3616] = "*INCLUDE, INPUT=" + viadukt_jobName + ".bsp\n"

inp_file_upd = open(jobName + '.inp', 'w')
inp_file_upd.writelines(inp_file)
inp_file_upd.close()

# Starting the job with a terminal command
os.system("abaqus job=" + jobName + " interactive")