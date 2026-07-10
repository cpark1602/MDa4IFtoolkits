import MDAnalysis as mda
import os
import MDa4IFtoolkits.kirkwood_gk_interface as kirkwood_gk_interface

def test_dipole_angles():
    w_path = "./tests/"
    u_if = mda.Universe(os.path.join(w_path,"run-pos.pdb"), os.path.join(w_path,"run-pos.dcd"))
    boxX = 48.57
    boxY = 15.667
    boxZ = 15.076
    box = [boxX, boxY, boxZ, 90, 90, 90]
    u_if.dimensions = box
    print_results_path='results'
    pbc=True
    dim='x'
    bin_size = 0.02
    selection1 = 'name O'; selection2 = 'name O'
    start_stop_step = [0, 3, 1]
    _u_if = kirkwood_gk_interface.Kirkwood_Gk(u_if, box, print_results_path, pbc, bin_size, dim, selection1, selection2, cutoff_IF = [0, 12], cutoff_BULK = [19, 28], start=start_stop_step[0], stop=start_stop_step[1], step=start_stop_step[2])
    Kirkwood_gr = _u_if.run()
