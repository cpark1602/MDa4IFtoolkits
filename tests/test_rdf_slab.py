import MDAnalysis as mda
import os
import MDa4IFtoolkits.rdf_slab as rdf_slab

def test_dipole_angles():
    w_path = "./tests/"
    u_if = mda.Universe(os.path.join(w_path,"run-pos.pdb"), os.path.join(w_path,"run-pos.dcd"))
    boxX = 48.57
    boxY = 15.667
    boxZ = 15.076
    box = [boxX, boxY, boxZ, 90, 90, 90]
    u_if.dimensions = box
    bin_size = 0.1
    IdentityA = "O"
    IdentityB = "O"
    ag1 = u_if.select_atoms(f'name {IdentityA}')
    ag2 = u_if.select_atoms(f'name {IdentityB}')
    
    slab_analyzer = rdf_slab.InterSlabRDF(u_if, cutoff_slab=[8, 12], binsize=bin_size, exclusion_block=[1, 1])
    
    print("Starting trajectory analysis loop...")
    for ts in u_if.trajectory:
        slab_analyzer._single_frame(ag1, ag2)
    
    slab_analyzer.conclude()
    
    bins = slab_analyzer.bins
    rdf_data = slab_analyzer.rdf_slab_global
