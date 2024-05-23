for idx, mode in enumerate(mode_list):
    # Measured frequencies 
    lamda_m[idx] = freqs_est[:, 1][f_idx]

    # Numerical frequencies 
    # Top match based on MMI
    lamda[0, idx] = results_MMI[mode][0, 2]

    # MAC value 
    lamda[0, idx + len(mode_list)] = results_MMI[mode][0, 1]

    # Target MAC value
    lamda_m[idx + len(mode_list)] = 1