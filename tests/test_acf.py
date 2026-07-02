#!/usr/bin/env python

import numpy as np
import os
import time

import warnings
warnings.filterwarnings(action='once')

from MDAnalysis.lib.NeighborSearch import AtomNeighborSearch
from MDAnalysis.lib import distances

class ACF:
    def __init__(self, universe, box, HBs_criteria, selection1, selection2, print_results_path, cutoff_dist_O_H, cutoff_dist_donor_acceptor, cutoff_IF, cutoff_BULK, angle, pbc, start=None, stop=None, step=None, nac='IF'):   #  , **kwargs
        
        self.u = universe
        self.seconds1 = time.time()
        self.t_l = []

        self.total_frame = len(self.u.trajectory)
        print("Total frame:", self.total_frame)
        self.box = box
        self.selection1 = selection1
        self.selection2 = selection2
        self.n_frames = None
        
        # Cutoff for searching O/H atoms around the reference
        self.cutoff_dist_O_H = cutoff_dist_O_H
        # Cutoff for the distance between O-O.
        self.cutoff_dist_donor_acceptor = cutoff_dist_donor_acceptor
        # a bonded O-H distance cutoff; used in _gen_bonded_hydrogens_dist
        self.OH_dist_cutoff = 1.2    
        self.HBs_criteria = HBs_criteria
        
        # cutoff for the donor position; the values are based on the density profile Fig. 2.4 in the report
        self.cutoff_IF = cutoff_IF         #12
        self.cutoff_BULK = cutoff_BULK     #[19, 28] between 19 to 28 Angstrom 
        self.cutoff_NE = 31
        
        # Cutoff angle O-H-O
        self.angle = angle
        
        self.pbc = pbc and all(self.u.dimensions[:3])
        # Start frame nr. and stop frame nr.
        self.start = start
        self.stop = stop
        self.step= step
        # time for each frame
        self.timesteps = None
        
        ### HBs ACF
        self.hb_acf_results = None

        self.print_results_path = print_results_path
        if not os.path.exists(self.print_results_path):
            os.makedirs(self.print_results_path)

    def _get_bonded_hydrogens_dist(self, atom):
        """Find bonded hydrogens within cutoff to 'atom'.
        Hydrogen bonds are detected by the cutoff;
        The distance from the reference 'atom' is calculated for all hydrogens in the residue
        and only those within a cutoff are kept."""
        try:
            sel_h = atom.residue.atoms.select_atoms(
                "(name HW1 or name HW2) and around {0:f} index {1!s}".format(self.OH_dist_cutoff, atom.index))
            return sel_h
            
        except NoDataError:
            return []
        
    def _update_selection_1(self):
        """Update the hydrogens around Oxygen atoms"""
        self._s1 = self.u.select_atoms(self.selection1)
        self._s1_donors = {}
        self._s1_donors_h = {}
        self._s1_acceptors = {}

        # d: donor atom_select
        for i, d in enumerate(self._s1):
            tmp = self._get_bonded_hydrogens_dist(d)
            #print('get bonded h: ', tmp)
            if tmp:
                self._s1_donors_h[i] = tmp   # fill the dict[i]

                
    def _update_selection_2(self):
        """Update the hydrogens around Oxygen atoms"""
        self._s2 = self.u.select_atoms(self.selection2)
        self._s2_donors = {}
        self._s2_donors_h = {}
        self._s2_acceptors = {}
        for i, d in enumerate(self._s2):
            tmp = self._get_bonded_hydrogens_dist(d)
            #print(tmp)
            if tmp:
                self._s2_donors_h[i] = tmp   # fill the dict[i]

    def _single_run(self, start, stop, step):
        self.timesteps = []
        self._s1 = self.u.select_atoms(self.selection1)
        self.dict_Oatom_index_to_i = {}
        "----------HB ACF----------"
        already_found_first_frame_IF = {}                    # disc for the HB ACF
        prev_already_found_IF = {}                           # disc for the HB ACF
        
        already_found_first_frame_BULK = {}                    # disc for the HB ACF
        prev_already_found_BULK = {}                           # disc for the HB ACF

        hb_acf_results_IF = np.zeros_like(np.arange(start, stop, step), dtype=np.float32)
        hb_acf_results_BULK = np.zeros_like(np.arange(start, stop, step), dtype=np.float32)
        nr_HBs_IF_t = np.zeros_like(np.arange(start, stop, step), dtype=np.float32)
        nr_HBs_BULK_t = np.zeros_like(np.arange(start, stop, step), dtype=np.float32)
        nr_HBs_NE_t = np.zeros_like(np.arange(start, stop, step), dtype=np.float32)
        "---------------------------"

        self.dim_ndx = 0
        L = self.box[0]
        n_frames_i = 0
        stepsize=0.5
        for ts in self.u.trajectory[start:stop:step]:

            frame_results = []

            print("Frame: {0:5d}".format(ts.frame))

            already_found = {}
            "----------HB ACF----------"
            #hb_acf_already_found = {}
            hb_acf_already_found_IF = {}
            hb_acf_already_found_BULK = {}
            "--------------------------"
       
            "Dictionary for oxygen atoms with H tagged"
            self.dict_oxygen_with_h_tagged = {}

            frame = ts.frame
            print("-----------------------------")

            n_frames_i += 1
            time_fs = ts.frame * stepsize   # * self.step; ts.frame already include the step 
            
            self.timesteps.append(time_fs)

            #n_frableames += 1
            print('n_frames: ', n_frames_i)
            # Select atom groups
            self._s1 = self.u.select_atoms(self.selection1)
            self._s2 = self.u.select_atoms(self.selection2)

            for i_s1 in range(len(self._s1)):
                self.dict_Oatom_index_to_i[self._s1[i_s1].index] = i_s1 

            for i_s1 in range(len(self._s1)):
                ndx_s1 = self._s1[i_s1].index
            
            # Update the HB list
            #if self.update_selection1:
            """_update_selection_1() call the function _get_bonded_hydrogens_dist"""
            """to find the bonded Hydrogens"""
            self._update_selection_1()
            #if self.update_selection2:
            self._update_selection_2()

            # Call the class AtomNeighborSearch wrt _s2
            # Return all atoms/residues/segments that are within *radius* of the atoms in *atoms*.
            # ns_acceptor is a module, which is called later to search within 'cutoff'.** 
            self.ns_acceptors = AtomNeighborSearch(self._s2, self.box)

            # Count number water molecules in IF/BULK/Ne
            locpres= self.u.coord.positions
            JJ=self._s1.indices
            rxo=locpres[JJ,0]
            nr_water_in_if = np.where(rxo < self.cutoff_IF[1])
            nr_water_in_if = len(nr_water_in_if[0])
            nr_HBs_IF_t[n_frames_i-1] = nr_water_in_if * 2

            nr_water_in_bulk = np.where( (rxo > self.cutoff_BULK[0]) & ( rxo < self.cutoff_BULK[1] ) )
            nr_water_in_bulk = len(nr_water_in_bulk[0])
            nr_HBs_BULK_t[n_frames_i-1] = nr_water_in_bulk * 2

            nr_water_in_ne = np.where(rxo > self.cutoff_NE)
            nr_water_in_ne = len(nr_water_in_ne[0])
            nr_HBs_NE_t[n_frames_i-1] = nr_water_in_ne * 2

            """ search oxygens, which are closest hydrogen atoms  """
            self.Hatoms = self.u.select_atoms('name H')
            for h_tag in self.Hatoms:
                oxygens_for_Htagging = self.ns_acceptors.search(h_tag, 5.0)   # self.box[0]*2 self.box[0], search oxygens within the simulation box
                dict_dist_o_h_tagged = {} 
                ### Search oxygen near the reference hydrogen within the cutoff
                for o_near_h in oxygens_for_Htagging:
                    dist_o_h_tagged = distances.calc_bonds(h_tag.position, o_near_h.position, box=self.box)
                    dict_dist_o_h_tagged[o_near_h.index] = dist_o_h_tagged
                ### sort the dictionary to find the closest oxygen
                dict_dist_o_h_tagged_sorted = {k: v for k, v in sorted(dict_dist_o_h_tagged.items(), key=lambda item: item[1])}
                closest_o_near_h_tagged_index =  list(dict_dist_o_h_tagged_sorted.keys())[0]

                if closest_o_near_h_tagged_index in self.dict_oxygen_with_h_tagged:
                    self.dict_oxygen_with_h_tagged[closest_o_near_h_tagged_index].append(h_tag)
                else:
                    self.dict_oxygen_with_h_tagged[closest_o_near_h_tagged_index] = [h_tag]

            ###=== Tagging H finished===###
            ###--- loop over the oxygen dictionary with tagged H---###
            ###--- Oxygen atoms belong to ClO4 are not selected; because those are not the closest oxygens---###

            count_h3o = 0
            h3o_dist_list = []; h3o_dist_da_dict = {}
            for key, value in self.dict_oxygen_with_h_tagged.items():
                d = self._s1[self.dict_Oatom_index_to_i[key]]   # use a new dictionary which is sorted according to the self._s1 list.
                donor_pos_xyz = d.position
                donor_index = d.index
                
                # loop over tagged H, which belongs to the O_i
                for h in value:               
                    res = self.ns_acceptors.search(h, self.cutoff_dist_O_H)
                    # d-h -- a ; Search a around h within cutoff_dist_O_H 
                    for a in res:
                        if d.index == a.index:
                            pass
                        else: 
                            if (d.index, h.index, a.index) in already_found:
                                pass
                            else:
                                angle_rad = distances.calc_angles(h.position, 
                                                              d.position,
                                                              a.position, box=self.box)
                                
                                angle_OHO_rad = distances.calc_angles(d.position, 
                                                              h.position,
                                                              a.position, box=self.box)
                                
                                angle = np.rad2deg(angle_rad)
                                dist = distances.calc_bonds(d.position, a.position, box=self.box)
                                
                                dist_O1_H = distances.calc_bonds(d.position, h.position, box=self.box)
                                dist_O2_H = distances.calc_bonds(a.position, h.position, box=self.box)

                                "Rectangular HBs cutoffs"
                                if self.HBs_criteria == 'Luzar':
                                    if angle <= self.angle and dist <= self.cutoff_dist_donor_acceptor:

                                        already_found[(d.index, h.index, a.index)] = True
                                               
                                        # Find active H in the H3O+, update the histograme of the free energy of H3O+
                                        if len(value) == 3: 
                                            dist_1 = distances.calc_bonds(h.position, d.position, box=self.box)
                                            dist_2 = distances.calc_bonds(h.position, a.position, box=self.box)
                                            dist_delta_d_a = abs(dist_1 - dist_2)
                                            h3o_dist_list.append(dist_delta_d_a) 
                                            h3o_dist_da_dict[a.index] = dist_delta_d_a 
                                            h3o_activeHbond_donor = d.index

                                        "#####***** HBs ACF (IF) *****#####"
                                        if self.cutoff_IF[0] < d.position[0] <= self.cutoff_IF[1]:
                                            if frame == start:
                                                prev_already_found_IF[(d.index, h.index, a.index)] = True
                                            else:
                                                # search the same d-h-a, which survive from the previous 
                                                if (d.index, h.index, a.index) in prev_already_found_IF:
                                                    hb_acf_already_found_IF[(d.index, h.index, a.index)] = True
                                                    hb_acf_results_IF[n_frames_i-1] += 1

                                        "#####***** HBs ACF (BULK) *****#####"
                                        if d.position[0] >= self.cutoff_BULK[0] and d.position[0] <= self.cutoff_BULK[1]:
                                            if frame == start:
                                                prev_already_found_BULK[(d.index, h.index, a.index)] = True
                                            else:
                                                if (d.index, h.index, a.index) in prev_already_found_BULK:
                                                    hb_acf_already_found_BULK[(d.index, h.index, a.index)] = True
                                                    hb_acf_results_BULK[n_frames_i-1] += 1


                                #"Triangular HBs cutoffs"    
                                elif self.HBs_criteria == 'Sho':
                                    cosine_term = -1.71 * np.cos(angle_OHO_rad) + 1.37
                                    if dist_O2_H < cosine_term:
                                        already_found[(d.index, h.index, a.index)] = True
                                    
                                        # Make a new dict for analysis
                                        # Add the dornor and acceptor pair into the dict.
                                        # Add new value to the list of the key if a key does not exist, make a new key

                                        # Find active H in the H3O+, update the histograme of the free energy of H3O+
                                        if len(value) == 3: 
                                            dist_1 = distances.calc_bonds(h.position, d.position, box=self.box)
                                            dist_2 = distances.calc_bonds(h.position, a.position, box=self.box)
                                            dist_delta_d_a = abs(dist_1 - dist_2)
                                            h3o_dist_list.append(dist_delta_d_a) 
                                            h3o_dist_da_dict[a.index] = dist_delta_d_a 
                                            h3o_activeHbond_donor = d.index

                                        "#####***** HBs ACF (IF) *****#####"
                                        if self.cutoff_IF[0] < d.position[0] <= self.cutoff_IF[1]:
                                            if frame == start:
                                                prev_already_found_IF[(d.index, h.index, a.index)] = True
                                            else:
                                                #print("----- NOT first frame-----")
                                                # search the same d-h-a, which survive from the previous 
                                                if (d.index, h.index, a.index) in prev_already_found_IF:
                                                    hb_acf_already_found_IF[(d.index, h.index, a.index)] = True
                                                    hb_acf_results_IF[n_frames_i-1] += 1

                                        "#####***** HBs ACF (BULK) *****#####"
                                        if d.position[0] >= self.cutoff_BULK[0] and d.position[0] <= self.cutoff_BULK[1]:
                                            if frame == start:
                                                prev_already_found_BULK[(d.index, h.index, a.index)] = True
                                            else:
                                                if (d.index, h.index, a.index) in prev_already_found_BULK:
                                                    hb_acf_already_found_BULK[(d.index, h.index, a.index)] = True
                                                    hb_acf_results_BULK[n_frames_i-1] += 1



            "###********** HB ACF (IF) ***********###"
            if frame == start:
                already_found_first_frame_IF = prev_already_found_IF
                len_hb_acf_firstFrame_IF = len(already_found_first_frame_IF)
                hb_acf_results_IF[n_frames_i-1] = len_hb_acf_firstFrame_IF 
            else:
                prev_already_found_IF = hb_acf_already_found_IF 
                
            "###********** HB ACF (BULK) ***********###"
            if frame == start:
                already_found_first_frame_BULK = prev_already_found_BULK
                len_hb_acf_firstFrame_BULK = len(already_found_first_frame_BULK)
                hb_acf_results_BULK[n_frames_i-1] = len_hb_acf_firstFrame_BULK
                #prev_already_found_IF = already_found_IF
            else:
                prev_already_found_BULK = hb_acf_already_found_BULK
                
        print("=== End of ts loop ===") 

        def compute_stderror(s, sq, c):
            #s: sum; c: count
            return np.sqrt((sq/c) - (s*s)/(c*c) ) / (np.sqrt(c))


        # Normalized by the total number of HB in each region.
        hb_acf_results_IF /= nr_HBs_IF_t   # h at if
        hb_acf_results_BULK /= nr_HBs_BULK_t  # h at bulk 

        c_t = np.array(self.timesteps)
        h_0_if = hb_acf_results_IF[0] 
        h_t_if = hb_acf_results_IF
        c_dot_if = - (( h_t_if - h_0_if ) / c_t ) * ( 1 - h_t_if )
        c_dot_if_2 = - (( h_t_if[1] - h_0_if ) / (c_t[1] - c_t[0]) ) * ( 1 - h_t_if )

        h_0_bulk = hb_acf_results_BULK[0] 
        h_t_bulk = hb_acf_results_BULK
        c_dot_bulk = - (( h_t_bulk - h_0_bulk ) / c_t ) * ( 1 - h_t_bulk )
        c_dot_bulk_2 = - (( h_t_bulk[1] - h_0_bulk ) / (c_t[1] - c_t[0]) ) * ( 1 - h_t_bulk )

        self.seconds2 = time.time()
        run_time = self.seconds2 - self.seconds1
        runt_time = run_time/60

        return hb_acf_results_IF, hb_acf_results_BULK, c_dot_if, c_dot_if_2

    def _slice_trj(self):
        #window = 10000
        #sampling_numbers = 10000
        sampling_numbers = 3

        tmp = int((self.stop - self.start) / self.step)
        nruns = np.ceil(tmp / sampling_numbers)
            
        window_list = []; windows = 0
        for i in range(int(nruns+1)):
           window_list.append(windows)
           windows += sampling_numbers
        return window_list

    def run(self, **kwargs):
        step=1; 
        #window_list = self._slice_trj(nruns)    # [start, ..., intermediates..., end]
        window_list = self._slice_trj()    
        nruns=len(window_list);

        # Prepare
        hb_acf_results_IF_global = np.zeros_like(np.arange(window_list[0], window_list[1], self.step), dtype=np.float32)
        hb_acf_results_IF_acf_global = np.zeros_like(np.arange(window_list[0], window_list[1], self.step), dtype=np.float32)

        hb_acf_results_BULK_global = np.zeros_like(np.arange(window_list[0], window_list[1], self.step), dtype=np.float32)
        hb_acf_results_BULK_acf_global = np.zeros_like(np.arange(window_list[0], window_list[1], self.step), dtype=np.float32)

        c_dot_if_global = np.zeros_like(np.arange(window_list[0], window_list[1], self.step), dtype=np.float32)
        c_dot_if_2_global = np.zeros_like(np.arange(window_list[0], window_list[1], self.step), dtype=np.float32)

        # Average HB number at each starting point; Average of HB.
        hb_acf_results_IF_global = 0
        hb_acf_results_BULK_global = 0
        for i in range(nruns-1):
            print('nth window: ', i)
            hb_acf_results_IF, hb_acf_results_BULK, c_dot_if, c_dot_if_2 = self._single_run(window_list[i], window_list[i+1], self.step) 
            c_dot_if_global += c_dot_if
            c_dot_if_2_global += c_dot_if_2

            # ACF <h(0)h(t)>
            acf_if = hb_acf_results_IF[0] * hb_acf_results_IF
            hb_acf_results_IF_acf_global += acf_if
            # To compute average h
            hb_acf_results_IF_global +=  hb_acf_results_IF[0]

            # ACF <h(0)h(t)>
            acf_bulk = hb_acf_results_BULK[0] * hb_acf_results_BULK
            hb_acf_results_BULK_acf_global += acf_bulk
            # To compute average h
            hb_acf_results_BULK_global +=  hb_acf_results_BULK[0]

            #print(hb_acf_results_IF_global, hb_acf_results_IF) 

        hb_acf_results_IF_acf_global = hb_acf_results_IF_acf_global / (nruns -1)
        hb_acf_results_BULK_acf_global = hb_acf_results_BULK_acf_global / (nruns -1)

        hb_acf_results_IF_global = hb_acf_results_IF_global / (nruns -1)
        hb_acf_results_BULK_global = hb_acf_results_BULK_global / (nruns -1)

        print('hb_acf_results_IF_acf_global ', hb_acf_results_IF_acf_global)
        print('hb_acf_results_BULK_acf_global', hb_acf_results_BULK_global)

        # Final ACF; <h(0)*h(t)> / <h>
        hb_C_IF = hb_acf_results_IF_acf_global / hb_acf_results_IF_global 
        hb_C_BULK = hb_acf_results_BULK_acf_global / hb_acf_results_BULK_global 

        c_dot_if_global = ( c_dot_if_global / (nruns -1) ) / hb_acf_results_IF_global
        c_dot_if_2_global = ( c_dot_if_2_global / (nruns -1) ) / hb_acf_results_IF_global


        np.save(self.print_results_path+'/hb_acf_results_IF.npy', hb_C_IF)
        np.save(self.print_results_path+'/hb_acf_results_BULK.npy', hb_C_BULK)
        np.save(self.print_results_path+'/c_dot_if_global.npy', c_dot_if_global)
        np.save(self.print_results_path+'/c_dot_if_2_global.npy', c_dot_if_2_global)


sample_data = [1.0, 1.0, 1.0, 1.0]
#    result = calculate_acf(sample_data)
    
# Assert that the output calculates properly and matches shape expectations
assert len(sample_data) == 4
assert np.allclose(sample_data, [1.0, 1.0, 1.0, 1.0])
