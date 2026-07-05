#!/usr/bin/env python

import MDAnalysis as mda
import scipy.constants
import numpy as np
import warnings
warnings.filterwarnings(action="once")
import logging

logger = logging.getLogger("MDAnalysis.analysis.hbonds")


def test_mass_density_profile():
    #####----- Trj path -----
    w_path = "./tests/"
    u_if = mda.Universe(os.path.join(w_path,"run-pos.pdb"), os.path.join(w_path,"run-pos.dcd"))
    
    print("total nr. of frame: ", len(u_if.trajectory))
    tot_frames = len(u_if.trajectory)
    boxX = 48.57
    boxY = 15.667
    boxZ = 15.076
    box = [boxX, boxY, boxZ, 90, 90, 90]
    u_if.dimensions = box
    "Velesco angle H-O-O 35 degrees"
    HBs_criteria_input = "Sho"  # Luzar: a rectangule; Sho: Triangle
    start_stop_step = [0, 3, 1]  # q0.0-region2-new/ddec
    
    chemisorbed_cutoff_O = 10.5
    chemisorbed_cutoff_H = 9.85
    # path_results='./results/atomic_charge/'   # where the NAC results are stored
    print_results_path = w_path + "/results/"  # To save the results
    
    import mass_density
    dim = "x"
    bin_size = 0.02
    pbc = True
    if_q0_nac = mass_density.Mass_density(
        u_if,
        box,
        print_results_path,
        pbc,
        bin_size,
        dim,
        start=start_stop_step[0],
        stop=start_stop_step[1],
        step=start_stop_step[2],
    )
    
    his_edges, number_density_O = if_q0_nac._get_densityProfile("O")
    his_edges, number_density_H = if_q0_nac._get_densityProfile("H")
    his_edges, number_density_Au = if_q0_nac._get_densityProfile("Au")
    his_edges, number_density_Ne = if_q0_nac._get_densityProfile("Ne")

    assert len(his_edges) > 0, "Histogram edges should not be empty"
    assert len(number_density_O) == len(his_edges) - 1, "Density profile length must match bins length"
