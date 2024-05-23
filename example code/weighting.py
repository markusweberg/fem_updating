for idx, mode in enumerate(mode_list):
    # Weighting matrix
    if idx < 4:
        W[idx, idx] = (2/3) / len(mode_list) + (1/40) * ((len(mode_list)-4)/len(mode_list))
        W[idx + len(mode_list), idx + len(mode_list)] = (1/3) / len(mode_list) + (1/40) * ((len(mode_list)-4)/len(mode_list))
    else:
        W[idx, idx] = (2/3) / len(mode_list) - (1/40) * (4/len(mode_list))
        W[idx + len(mode_list), idx + len(mode_list)] = (1/3) / len(mode_list) - (1/40) * (4/len(mode_list))