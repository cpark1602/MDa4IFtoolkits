#!/usr/bin/env python
import os
import time
import warnings
import numpy as np
from MDAnalysis.lib.NeighborSearch import AtomNeighborSearch
from MDAnalysis.lib import distances
from MDAnalysis import NoDataError

# Suppress repetitive warnings during continuous trajectory parsing loops
warnings.filterwarnings(action="once")


class ACF:
    """Computes the Hydrogen Bond Autocorrelation Function (ACF) for interfacial

    and bulk regions using specified criteria (Luzar or Sho).
    """

    def __init__(
        self,
        universe,
        box,
        HBs_criteria,
        selection1,
        selection2,
        print_results_path,
        cutoff_dist_O_H,
        cutoff_dist_donor_acceptor,
        cutoff_IF,
        cutoff_BULK,
        angle,
        pbc,
        start=None,
        stop=None,
        step=None,
        nac="IF",
    ):
        """Initializes the ACF analyzer with simulation universes, geometries,

        and hydrogen bond criteria cutoffs.
        """
        self.u = universe
        self.seconds1 = time.time()  # Track the analytical start time
        self.t_l = []

        self.total_frame = len(self.u.trajectory)
        print("Total frame:", self.total_frame)
        self.box = box
        self.selection1 = selection1
        self.selection2 = selection2
        self.n_frames = None

        # Cutoff parameters for identifying eligible H-bond components
        self.cutoff_dist_O_H = cutoff_dist_O_H  # Max distance for O...H pairs
        self.cutoff_dist_donor_acceptor = (
            cutoff_dist_donor_acceptor  # Max distance for O...O pairs
        )
        self.OH_dist_cutoff = 1.2  # Hard cutoff defining structural O-H covalent binding
        self.HBs_criteria = HBs_criteria  # Geometry switch: 'Luzar' or 'Sho'

        # Interfacial spatial slicing regions based on system density profile boundaries
        self.cutoff_IF = cutoff_IF  # Interfacial zone boundaries [min, max]
        self.cutoff_BULK = cutoff_BULK  # Bulk liquid region boundaries [min, max]
        self.cutoff_NE = 31  # Coordinate limit identifying non-electrolyte/vacuum zone

        self.angle = angle  # Max angle cutoff for H-bond tracking validation
        self.pbc = pbc and all(
            self.u.dimensions[:3]
        )  # Confirm true Periodic Boundary Conditions status

        # Trajectory slice boundaries
        self.start = start
        self.stop = stop
        self.step = step
        self.timesteps = None

        self.hb_acf_results = None
        self.print_results_path = print_results_path
        if not os.path.exists(self.print_results_path):
            os.makedirs(self.print_results_path)

    def _get_bonded_hydrogens_dist(self, atom):
        """Finds covalently bonded hydrogens within a rigid proximity cutoff of

        the reference Oxygen atom.

        Args:
            atom (Atom): MDAnalysis Atom object representing the reference oxygen.

        Returns:
            AtomGroup or list: Selected hydrogen atoms inside the covalent radius,
                              or an empty list if none are resolved.
        """
        try:
            # Query hydrogens belonging strictly to the same local water residue
            sel_h = atom.residue.atoms.select_atoms(
                "(name HW1 or name HW2) and around {0:f} index {1!s}".format(
                    self.OH_dist_cutoff, atom.index
                )
            )
            return sel_h
        except NoDataError:
            return []

    def _update_selection_1(self):
        """Updates and queries the AtomGroups and structural hydrogen lookups

        for Selection 1 (Donors).
        """
        self._s1 = self.u.select_atoms(self.selection1)
        self._s1_donors = {}
        self._s1_donors_h = {}
        self._s1_acceptors = {}

        # Loop over all selected oxygen atoms to identify intramolecular bonded hydrogens
        for i, d in enumerate(self._s1):
            tmp = self._get_bonded_hydrogens_dist(d)
            if tmp:
                self._s1_donors_h[i] = tmp  # Map local donor index to hydrogen group

    def _update_selection_2(self):
        """Updates and queries the AtomGroups and structural hydrogen lookups

        for Selection 2 (Acceptors).
        """
        self._s2 = self.u.select_atoms(self.selection2)
        self._s2_donors = {}
        self._s2_donors_h = {}
        self._s2_acceptors = {}

        for i, d in enumerate(self._s2):
            tmp = self._get_bonded_hydrogens_dist(d)
            if tmp:
                self._s2_donors_h[i] = tmp  # Map local acceptor index to hydrogen group

    def _single_run(self, start, stop, step):
        """Processes a single trajectory window block to trace time-correlated

        H-bond survival probabilities.

        Args:
            start (int): Starting frame index of the trajectory block.
            stop (int): Ending frame index of the trajectory block.
            step (int): Trajectory sampling stride step size.

        Returns:
            tuple: Arrays containing raw H-bond populations and continuous decay indicators.
        """
        self.timesteps = []
        self._s1 = self.u.select_atoms(self.selection1)
        self.dict_Oatom_index_to_i = {}

        # Track history dictionaries to check H-bond survival relative to the time origin
        already_found_first_frame_IF = {}
        prev_already_found_IF = {}
        already_found_first_frame_BULK = {}
        prev_already_found_BULK = {}

        # Initialize result storage arrays matching the window frame count size
        frame_range_len = len(np.arange(start, stop, step))
        hb_acf_results_IF = np.zeros(frame_range_len, dtype=np.float32)
        hb_acf_results_BULK = np.zeros(frame_range_len, dtype=np.float32)

        # Track populations to normalize configurations based on regional capacities
        nr_HBs_IF_t = np.zeros(frame_range_len, dtype=np.float32)
        nr_HBs_BULK_t = np.zeros(frame_range_len, dtype=np.float32)
        nr_HBs_NE_t = np.zeros(frame_range_len, dtype=np.float32)

        self.dim_ndx = 0
        n_frames_i = 0
        stepsize = 0.5  # Simulation time calibration step size parameter (in fs/ps)

        # Iterate step-by-step through the designated slice window of the simulation trajectory
        for ts in self.u.trajectory[start:stop:step]:
            print("Frame: {0:5d}".format(ts.frame))

            already_found = {}
            hb_acf_already_found_IF = {}
            hb_acf_already_found_BULK = {}

            # Dictionary to dynamically pair specific hydrogen objects to nearest parent oxygens
            self.dict_oxygen_with_h_tagged = {}
            frame = ts.frame
            print("-----------------------------")

            n_frames_i += 1
            time_fs = ts.frame * stepsize
            self.timesteps.append(time_fs)

            print("n_frames: ", n_frames_i)
            self._s1 = self.u.select_atoms(self.selection1)
            self._s2 = self.u.select_atoms(self.selection2)

            # Generate quick reverse lookup index mapping for fast selection matching
            for i_s1 in range(len(self._s1)):
                self.dict_Oatom_index_to_i[self._s1[i_s1].index] = i_s1

            # Re-evaluate covalent topology selections for the current frame step
            self._update_selection_1()
            self._update_selection_2()

            # Initialize fast tree search structure for target H-bond acceptors
            self.ns_acceptors = AtomNeighborSearch(self._s2, self.box)

            # Count water distributions across spatial slabs using primary coordinate arrays
            locpres = self.u.coord.positions
            JJ = self._s1.indices
            rxo = locpres[JJ, 0]  # Read coordinates along the target surface normal (X-axis)

            # Count population inside the Interface (IF) boundaries
            nr_water_in_if = len(np.where(rxo < self.cutoff_IF[1])[0])
            nr_HBs_IF_t[n_frames_i - 1] = nr_water_in_if * 2

            # Count population inside the Bulk zone boundaries
            nr_water_in_bulk = len(
                np.where(
                    (rxo > self.cutoff_BULK[0]) & (rxo < self.cutoff_BULK[1])
                )[0]
            )
            nr_HBs_BULK_t[n_frames_i - 1] = nr_water_in_bulk * 2

            # Count population inside the Non-Electrolyte (NE) zone boundaries
            nr_water_in_ne = len(np.where(rxo > self.cutoff_NE)[0])
            nr_HBs_NE_t[n_frames_i - 1] = nr_water_in_ne * 2

            # Dynamic assignment step: Tag every hydrogen to its closest parent oxygen
            self.Hatoms = self.u.select_atoms("name H")
            for h_tag in self.Hatoms:
                # Search all nearby oxygens within an expanded 5.0 Å bounding check area
                oxygens_for_Htagging = self.ns_acceptors.search(h_tag, 5.0)
                dict_dist_o_h_tagged = {}

                for o_near_h in oxygens_for_Htagging:
                    dist_o_h_tagged = distances.calc_bonds(
                        h_tag.position, o_near_h.position, box=self.box
                    )
                    dict_dist_o_h_tagged[o_near_h.index] = dist_o_h_tagged

                # Sort the distance map to resolve the authentic parent oxygen index
                dict_dist_o_h_tagged_sorted = {
                    k: v
                    for k, v in sorted(
                        dict_dist_o_h_tagged.items(), key=lambda item: item[1]
                    )
                }
                closest_o_near_h_tagged_index = list(
                    dict_dist_o_h_tagged_sorted.keys()
                )[0]

                if closest_o_near_h_tagged_index in self.dict_oxygen_with_h_tagged:
                    self.dict_oxygen_with_h_tagged[
                        closest_o_near_h_tagged_index
                    ].append(h_tag)
                else:
                    self.dict_oxygen_with_h_tagged[
                        closest_o_near_h_tagged_index
                    ] = [h_tag]

            # Hydrogen bonding verification step
            h3o_dist_list = []
            h3o_dist_da_dict = {}
            for key, value in self.dict_oxygen_with_h_tagged.items():
                d = self._s1[self.dict_Oatom_index_to_i[key]]

                # Evaluate all target hydrogens associated with parent donor molecule d
                for h in value:
                    # Look for nearby H-bond acceptors (a) surrounding the active hydrogen
                    res = self.ns_acceptors.search(h, self.cutoff_dist_O_H)
                    for a in res:
                        if d.index == a.index:
                            pass  # Skip intramolecular assignments
                        else:
                            if (d.index, h.index, a.index) in already_found:
                                pass  # Skip duplicate tracking entries
                            else:
                                # Compute internal angle arrays using minimum image convention
                                angle_rad = distances.calc_angles(
                                    h.position,
                                    d.position,
                                    a.position,
                                    box=self.box,
                                )
                                angle_OHO_rad = distances.calc_angles(
                                    d.position,
                                    h.position,
                                    a.position,
                                    box=self.box,
                                )

                                angle_deg = np.rad2deg(angle_rad)
                                dist = distances.calc_bonds(
                                    d.position, a.position, box=self.box
                                )
                                dist_O2_H = distances.calc_bonds(
                                    a.position, h.position, box=self.box
                                )

                                # --- Condition A: Geometric Luzar H-Bond Criteria ---
                                if self.HBs_criteria == "Luzar":
                                    if (
                                        angle_deg <= self.angle
                                        and dist
                                        <= self.cutoff_dist_donor_acceptor
                                    ):
                                        already_found[
                                            (d.index, h.index, a.index)
                                        ] = True

                                        # Optional: Analyze extra bonding properties if hydronium (H3O+) states exist
                                        if len(value) == 3:
                                            dist_1 = distances.calc_bonds(
                                                h.position,
                                                d.position,
                                                box=self.box,
                                            )
                                            dist_2 = distances.calc_bonds(
                                                h.position,
                                                a.position,
                                                box=self.box,
                                            )
                                            dist_delta_d_a = abs(
                                                dist_1 - dist_2
                                            )
                                            h3o_dist_list.append(dist_delta_d_a)
                                            h3o_dist_da_dict[a.index] = (
                                                dist_delta_d_a
                                            )

                                        # Interfacial (IF) track validation & persistence checking
                                        if (
                                            self.cutoff_IF[0]
                                            < d.position[0]
                                            <= self.cutoff_IF[1]
                                        ):
                                            if frame == start:
                                                prev_already_found_IF[
                                                    (d.index, h.index, a.index)
                                                ] = True
                                            else:
                                                if (
                                                    d.index,
                                                    h.index,
                                                    a.index,
                                                ) in prev_already_found_IF:
                                                    hb_acf_already_found_IF[
                                                        (
                                                            d.index,
                                                            h.index,
                                                            a.index,
                                                        )
                                                    ] = True
                                                    hb_acf_results_IF[
                                                        n_frames_i - 1
                                                    ] += 1

                                        # Bulk zone track validation & persistence checking
                                        if (
                                            d.position[0] >= self.cutoff_BULK[0]
                                            and d.position[0]
                                            <= self.cutoff_BULK[1]
                                        ):
                                            if frame == start:
                                                prev_already_found_BULK[
                                                    (d.index, h.index, a.index)
                                                ] = True
                                            else:
                                                if (
                                                    d.index,
                                                    h.index,
                                                    a.index,
                                                ) in prev_already_found_BULK:
                                                    hb_acf_already_found_BULK[
                                                        (
                                                            d.index,
                                                            h.index,
                                                            a.index,
                                                        )
                                                    ] = True
                                                    hb_acf_results_BULK[
                                                        n_frames_i - 1
                                                    ] += 1

                                # --- Condition B: Non-linear/Triangular Sho Criteria ---
                                elif self.HBs_criteria == "Sho":
                                    cosine_term = (
                                        -1.71 * np.cos(angle_OHO_rad) + 1.37
                                    )
                                    if dist_O2_H < cosine_term:
                                        already_found[
                                            (d.index, h.index, a.index)
                                        ] = True

                                        if len(value) == 3:
                                            dist_1 = distances.calc_bonds(
                                                h.position,
                                                d.position,
                                                box=self.box,
                                            )
                                            dist_2 = distances.calc_bonds(
                                                h.position,
                                                a.position,
                                                box=self.box,
                                            )
                                            dist_delta_d_a = abs(
                                                dist_1 - dist_2
                                            )
                                            h3o_dist_list.append(dist_delta_d_a)
                                            h3o_dist_da_dict[a.index] = (
                                                dist_delta_d_a
                                            )

                                        # Interfacial (IF) persistence tracking (Sho model)
                                        if (
                                            self.cutoff_IF[0]
                                            < d.position[0]
                                            <= self.cutoff_IF[1]
                                        ):
                                            if frame == start:
                                                prev_already_found_IF[
                                                    (d.index, h.index, a.index)
                                                ] = True
                                            else:
                                                if (
                                                    d.index,
                                                    h.index,
                                                    a.index,
                                                ) in prev_already_found_IF:
                                                    hb_acf_already_found_IF[
                                                        (
                                                            d.index,
                                                            h.index,
                                                            a.index,
                                                        )
                                                    ] = True
                                                    hb_acf_results_IF[
                                                        n_frames_i - 1
                                                    ] += 1

                                        # Bulk zone persistence tracking (Sho model)
                                        if (
                                            d.position[0] >= self.cutoff_BULK[0]
                                            and d.position[0]
                                            <= self.cutoff_BULK[1]
                                        ):
                                            if frame == start:
                                                prev_already_found_BULK[
                                                    (d.index, h.index, a.index)
                                                ] = True
                                            else:
                                                if (
                                                    d.index,
                                                    h.index,
                                                    a.index,
                                                ) in prev_already_found_BULK:
                                                    hb_acf_already_found_BULK[
                                                        (
                                                            d.index,
                                                            h.index,
                                                            a.index,
                                                        )
                                                    ] = True
                                                    hb_acf_results_BULK[
                                                        n_frames_i - 1
                                                    ] += 1

            # Finalize block state storage transitions at the time origin boundary
            if frame == start:
                already_found_first_frame_IF = prev_already_found_IF
                len_hb_acf_firstFrame_IF = len(already_found_first_frame_IF)
                hb_acf_results_IF[n_frames_i - 1] = len_hb_acf_firstFrame_IF
            else:
                prev_already_found_IF = hb_acf_already_found_IF

            if frame == start:
                already_found_first_frame_BULK = prev_already_found_BULK
                len_hb_acf_firstFrame_BULK = len(already_found_first_frame_BULK)
                hb_acf_results_BULK[n_frames_i - 1] = len_hb_acf_firstFrame_BULK
            else:
                prev_already_found_BULK = hb_acf_already_found_BULK

        print("=== End of ts loop ===")

        def compute_stderror(s, sq, c):
            """Calculates the statistical standard error from cumulative sum

            variances.
            """
            return np.sqrt((sq / c) - (s * s) / (c * c)) / (np.sqrt(c))

        # Normalize remaining survival configurations by regional densities
        hb_acf_results_IF /= nr_HBs_IF_t
        hb_acf_results_BULK /= nr_HBs_BULK_t

        # Evaluate differential decay kinetics curves
        c_t = np.array(self.timesteps)
        h_0_if = hb_acf_results_IF[0]
        h_t_if = hb_acf_results_IF
        c_dot_if = -((h_t_if - h_0_if) / c_t) * (1 - h_t_if)
        c_dot_if_2 = -((h_t_if[1] - h_0_if) / (c_t[1] - c_t[0])) * (1 - h_t_if)

        h_0_bulk = hb_acf_results_BULK[0]
        h_t_bulk = hb_acf_results_BULK
        c_dot_bulk = -((h_t_bulk - h_0_bulk) / c_t) * (1 - h_t_bulk)
        c_dot_bulk_2 = -((h_t_bulk[1] - h_0_bulk) / (c_t[1] - c_t[0])) * (
            1 - h_t_bulk
        )

        self.seconds2 = time.time()
        run_time = self.seconds2 - self.seconds1
        runt_time = run_time / 60

        return hb_acf_results_IF, hb_acf_results_BULK, c_dot_if, c_dot_if_2

    def _slice_trj(self):
        """Slices the trajectory range into block windows to apply ensemble

        averages across distinct time origins.

        Returns:
            list: Index markers defining sequential block boundaries.
        """
        sampling_numbers = 3  # Frames per averaging block window
        tmp = int((self.stop - self.start) / self.step)
        nruns = np.ceil(tmp / sampling_numbers)

        window_list = []
        windows = 0
        for i in range(int(nruns + 1)):
            window_list.append(windows)
            windows += sampling_numbers
        return window_list

    def run(self, **kwargs):
        """Executes the complete window slicing ensemble routine to average the

        interfacial and bulk H-bond ACFs.
        """
        step = 1
        window_list = self._slice_trj()
        nruns = len(window_list)

        # Allocate global tracking arrays corresponding to singular block window sizes
        window_size = len(
            np.arange(window_list[0], window_list[1], self.step)
        )
        hb_acf_results_IF_acf_global = np.zeros(window_size, dtype=np.float32)
        hb_acf_results_BULK_acf_global = np.zeros(window_size, dtype=np.float32)
        c_dot_if_global = np.zeros(window_size, dtype=np.float32)
        c_dot_if_2_global = np.zeros(window_size, dtype=np.float32)

        hb_acf_results_IF_global = 0
        hb_acf_results_BULK_global = 0

        # Run sequential calculation passes across the trajectory block segments
        for i in range(nruns - 1):
            print("nth window: ", i)
            hb_acf_results_IF, hb_acf_results_BULK, c_dot_if, c_dot_if_2 = (
                self._single_run(window_list[i], window_list[i + 1], self.step)
            )
            c_dot_if_global += c_dot_if
            c_dot_if_2_global += c_dot_if_2

            # Compute actual autocorrelation statistics: <h(0) . h(t)>
            acf_if = hb_acf_results_IF[0] * hb_acf_results_IF
            hb_acf_results_IF_acf_global += acf_if
            hb_acf_results_IF_global += hb_acf_results_IF[0]

            acf_bulk = hb_acf_results_BULK[0] * hb_acf_results_BULK
            hb_acf_results_BULK_acf_global += acf_bulk
            hb_acf_results_BULK_global += hb_acf_results_BULK[0]

        # Calculate averages across all processed slice blocks
        hb_acf_results_IF_acf_global /= nruns - 1
        hb_acf_results_BULK_acf_global /= nruns - 1
        hb_acf_results_IF_global /= nruns - 1
        hb_acf_results_BULK_global /= nruns - 1

        # Final step normalization: Compute final correlation values C(t) = <h(0)*h(t)> / <h>
        hb_C_IF = hb_acf_results_IF_acf_global / hb_acf_results_IF_global
        hb_C_BULK = hb_acf_results_BULK_acf_global / hb_acf_results_BULK_global

        c_dot_if_global = (c_dot_if_global / (nruns - 1)) / hb_acf_results_IF_global
        c_dot_if_2_global = (
            c_dot_if_2_global / (nruns - 1)
        ) / hb_acf_results_IF_global

        # Un-comment the lines below whenever writing data straight to disk files is needed:
        # np.save(self.print_results_path + "/hb_acf_results_IF.npy", hb_C_IF)
        # np.save(self.print_results_path + "/hb_acf_results_BULK.npy", hb_C_BULK)
