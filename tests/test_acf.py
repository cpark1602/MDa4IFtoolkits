import MDAnalysis as mda
import os
import MDa4IFtoolkits.acf as acf

def test_acf():
    w_path = "./tests/"
    u_if = mda.Universe(os.path.join(w_path,"run-pos.pdb"), os.path.join(w_path,"run-pos.dcd"))
    boxX = 48.57
    boxY = 15.667
    boxZ = 15.076
    box = [boxX, boxY, boxZ, 90, 90, 90]
    u_if.dimensions = box
    start_stop_step = [0, 3, 1]
    print_results_path = "/results/"
    HBs_criteria_input = 'Sho'
    
    _u_if = acf.ACF(u_if, box, HBs_criteria_input, 'name O', 'name O', print_results_path, cutoff_dist_O_H =3.5, cutoff_dist_donor_acceptor = 3.5, cutoff_IF = [0, 12], cutoff_BULK = [19, 28], angle=35.0, pbc=True, start=start_stop_step[0], stop=start_stop_step[1], step=start_stop_step[2]) 
    
    _u_if.run()
