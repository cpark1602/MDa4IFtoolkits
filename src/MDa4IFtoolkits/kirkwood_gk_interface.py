# Copyright (c) 2026 Chanbum Park chanbum.park@theochem.ruhr-uni-bochum.de 
# Distributed under the terms of the GNU General Public License.

import time
import logging
import warnings
import numpy as np
from MDAnalysis.lib.NeighborSearch import AtomNeighborSearch
from MDAnalysis.lib import distances

import MDa4IFtoolkits.kirkwood_gk as kw_gk

# Setup warnings and logging filter configurations to avoid stdout pollution
warnings.filterwarnings(action='once')
logger = logging.getLogger('MDAnalysis.analysis.hbonds')

class Kirkwood_Gk:
    """
    Trajectory analyzer managing data workflow chunk scheduling, spatial filtering, 
    and mapping of raw hydrogen atom groups to explicit parent water oxygen atoms 
    for multi-layer local Kirkwood Gk factor processing.
    """
    def __init__(self, universe, box,  
                 print_results_path, pbc, bin_size, dim,
                 selection1, selection2, cutoff_IF, cutoff_BULK,  
                 start=None, stop=None, step=None): 
        
        self.u = universe
        self.seconds1 = time.time()  # Benchmark start timestamp

        # --- Trajectory Slicing Parameters ---
        self.total_frame = len(self.u.trajectory)
        print("Total frame:", self.total_frame)
        self.box = box
        self.selection1 = selection1  # Target water group selection string (e.g., 'resname SOL and name OW')
        self.selection2 = selection2  # Reference neighbor search base selection string
       
        # --- Spatial & Coordination Cutoffs (Angstroms) ---
        self.OH_dist_cutoff = 1.2     # Max intra-molecular O-H covalent bond distance
        self.cutoff_IF = cutoff_IF     # Spatial bounds [min, max] identifying the Interface region
        self.cutoff_BULK = cutoff_BULK # Spatial bounds [min, max] identifying the Bulk solvent core
        self.cutoff_NE = 31            # Spatial threshold bound mapping Far-End layer zones
        
        # --- Simulation Setup ---
        self.pbc = pbc and all(self.u.dimensions[:3])
        self.start = start if start is not None else 0
        self.stop = stop if stop is not None else self.total_frame
        self.step = step if step is not None else 1
        
        # --- Output Target Validation ---
        self.print_results_path = print_results_path
        #if not os.path.exists(self.print_results_path):
        #    os.makedirs(self.print_results_path)

    def _get_bonded_hydrogens_dist(self, atom):
        """
        Extracts structural target water hydrogens within a close distance cutoff 
        of a specific, parent oxygen atom.
        """
        try:
            # Query residue scope via geographic selection search query
            sel_h = atom.residue.atoms.select_atoms(
                "(name HW1 or name HW2) and around {0:f} index {1!s}".format(self.OH_dist_cutoff, atom.index))
            return sel_h
        except Exception:
            return []
        
    def _update_selection_1(self):
        """
        Loops through active selection targets to dynamically re-build local covalent dictionary 
        mappings for the currently selected oxygen atom indices.
        """
        self._s1 = self.u.select_atoms(self.selection1)
        self._s1_donors_h = {}

        for i, d in enumerate(self._s1):
            tmp = self._get_bonded_hydrogens_dist(d)
            if tmp:
                self._s1_donors_h[i] = tmp

    def _run_window_i(self, start, stop, step):
        """
        Core sub-trajectory processing iteration loop. Sets up neighborhood spatial lookup grid trees 
        and extracts geometric assignments for each target simulation frame.
        """
        self._s1 = self.u.select_atoms(self.selection1)
        self.dict_Oatom_index_to_i = {}

        n_frames_i = 0
        self.u.dimensions = self.box  # Enforce runtime system unit dimensions bounds

        logger.info("HBond analysis: starting window step loop")

        # Instantiate underlying multi-dimensional core calculation calculator
        kw_gk_u = kw_gk.Kirkwood_Gk(self.u)

        # Loop through a chunk window slice of the trajectory
        for ts in self.u.trajectory[start:stop:step]:
            frame = ts.frame
            print(f"Frame: {frame:5d}\n-----------------------------")

            self.dict_oxygen_with_h_tagged = {}
            n_frames_i += 1
            
            # Re-evaluate dynamic atom group assignments per frame
            self._s1 = self.u.select_atoms(self.selection1)
            self._s2 = self.u.select_atoms(self.selection2)

            # Map the absolute MDAnalysis system index to sequential loop indices
            for i_s1 in range(len(self._s1)):
                self.dict_Oatom_index_to_i[self._s1[i_s1].index] = i_s1 

            self._update_selection_1()
            
            # Generate high-performance grid data lookup tree for neighbor queries
            self.ns_acceptors = AtomNeighborSearch(self._s2, self.box)

            # --- Spatial Density and Count Filtering ---
            locpres = self.u.coord.positions
            JJ = self._s1.indices
            rxo = locpres[JJ, 0]  # Isolate coordinate coordinates along the X-axis
            
            # Determine density threshold count in interface layer bounds
            nr_water_in_if = len(np.where(rxo < self.cutoff_IF[1])[0])

            # Process dipole profiles only if density parameters satisfy solvent limits
            if nr_water_in_if > 2:
                self.Hatoms = self.u.select_atoms('name H')
                
                # --- Map Hydrogens to Closest Parent Oxygens ---
                for h_tag in self.Hatoms:
                    # Query neighbor lookup tree for all oxygens within 5.0 Angstroms
                    oxygens_for_Htagging = self.ns_acceptors.search(h_tag, 5.0)   
                    dict_dist_o_h_tagged = {} 
                    
                    # Compute true minimum distance taking PBC box wrap metrics into account
                    for o_near_h in oxygens_for_Htagging:
                        dist_o_h_tagged = distances.calc_bonds(h_tag.position, o_near_h.position, box=self.box)
                        dict_dist_o_h_tagged[o_near_h.index] = dist_o_h_tagged
                    
                    # Sort by distance and pop out the absolute closest oxygen index assignment
                    dict_dist_o_h_tagged_sorted = {k: v for k, v in sorted(dict_dist_o_h_tagged.items(), key=lambda item: item[1])}
                    closest_o_near_h_tagged_index = list(dict_dist_o_h_tagged_sorted.keys())[0]

                    # Append the checked hydrogen to its corresponding unique oxygen molecular map entry
                    if closest_o_near_h_tagged_index in self.dict_oxygen_with_h_tagged:
                        self.dict_oxygen_with_h_tagged[closest_o_near_h_tagged_index].append(h_tag)
                    else:
                        self.dict_oxygen_with_h_tagged[closest_o_near_h_tagged_index] = [h_tag]

                # --- Execute Local Kirkwood Gk Processing ---
                # Iterate through targeted spherical structural lookup coordination shells
                for kw_gk_cutoff in [0.5, 1, 2]: 
                    kw_gk_u.run(ts, n_frames_i, self._s1, self.dict_Oatom_index_to_i, 
                                self.dict_oxygen_with_h_tagged, kw_gk_cutoff, self.box)

        print("=== End of ts window loop ===") 
        
        # Extrapolate normalized system average properties evaluated across this chunk
        kw_gk_mu_aver_global = kw_gk_u.mu_aver / kw_gk_u.count  
        return kw_gk_mu_aver_global

    def _slice_trj(self):
        """
        Splits global trajectory parameters up into standardized block tracking arrays 
        containing frame index boundary indicators spaced by 100-frame chunks.
        """
        sampling_numbers = 3 
        tmp = int((self.stop - self.start) / self.step)
        nruns = np.ceil(tmp / sampling_numbers)

        window_list = []
        windows = 0
        for i in range(int(nruns + 1)):
            window_list.append(windows)
            windows += sampling_numbers
        return window_list

    def run(self, **kwargs):
        """
        Main execution manager endpoint handling trajectory slicing logic, 
        running the window processing blocks, and tracking calculation execution speeds.
        """
        window_list = self._slice_trj()
        nruns = len(window_list)
        print('Total window segments to process: ', nruns - 1)

        kw_gk_mu_aver_global = 0.0

        # Run through data segments individually to keep memory overhead clean
        for i in range(nruns - 1):
            print('Processing window block: ', i)
            kw_gk_mu_aver = self._run_window_i(window_list[i], window_list[i+1], self.step) 
            kw_gk_mu_aver_global += kw_gk_mu_aver

        # Final normalization calculation step over processed execution sets
        kw_gk_mu_aver_global /= (nruns - 1)

        # --- Benchmark Diagnostics ---
        self.seconds2 = time.time()
        run_time = self.seconds2 - self.seconds1
        print(f"Analysis complete. Total runtime: {run_time / 60.0:.2f} Minutes")
        return kw_gk_mu_aver_global
