import os
import numpy as np
from MDAnalysis import NoDataError
from MDAnalysis.lib import distances
import math as m


class Dangling_bonds:
    """
    A class to analyze the structural orientations of solvent molecules (e.g., water)
    near interfaces using MD trajectories. It calculates dipole vectors, OH bond vectors, 
    and classifies molecules into distinct regions (Interface, Bulk, etc.) along the X-axis.
    """
    def __init__(
        self,
        universe,
        box,
        print_results_path,
        pbc,
        bin_size,
        dim,
        selection1,
        cutoff_IF,
        cutoff_BULK,
        start=None,
        stop=None,
        step=None,
    ):
        # ---------------------------------------------------------------------
        # Simulation Universe & Geometric Parameters
        # ---------------------------------------------------------------------
        self.bin_size = bin_size
        self.u = universe
        self.pbc = pbc
        self.dim = dim
        self.box = box  # System dimensions array [Lx, Ly, Lz, alpha, beta, gamma]

        # ---------------------------------------------------------------------
        # Output directory configuration
        # ---------------------------------------------------------------------
        self.print_results_path = print_results_path
        #if not os.path.exists(self.print_results_path):
        #    os.makedirs(self.print_results_path)

        # ---------------------------------------------------------------------
        # Analysis Cutoffs & Trajectory Slabbing Boundaries
        # ---------------------------------------------------------------------
        self.OH_dist_cutoff = 1.2    # Intramolecular O-H distance cutoff (Angstroms)
        self.cutoff_IF = cutoff_IF      # [X_min, X_max] boundary for the Interface region
        self.cutoff_BULK = cutoff_BULK  # [X_min, X_max] boundary for the Bulk region

        # Check PBC conditions across the three dimensions
        self.pbc = pbc and all(self.u.dimensions[:3])
        self.start = start
        self.stop = stop
        self.step = step
        self.start_end_skip = [start, stop, step]

        self.selection1 = selection1    # Atom selection string for the central atoms (typically Oxygen)

        # =====================================================================
        # Storage Lists for Orientational Results (Angles & Direction Cosines)
        # =====================================================================
        
        # --- Dipole Vectors Orientation Arrays (Degrees) ---
        self.valdip_IF_xaxis = []
        self.valdip_IM_xaxis = []
        self.valdip_BULK_xaxis = []
        self.valdip_NE_xaxis = []

        self.valdip_IF_yaxis = []
        self.valdip_IM_yaxis = []
        self.valdip_BULK_yaxis = []
        self.valdip_NE_yaxis = []

        self.valdip_IF_zaxis = []
        self.valdip_IM_zaxis = []
        self.valdip_BULK_zaxis = []
        self.valdip_NE_zaxis = []

        # --- Dipole Orientation Arrays (Cosine Theta / Direction Cosines) ---
        self.cosTheta_valdip_IF_xaxis = []
        self.cosTheta_valdip_IM_xaxis = []
        self.cosTheta_valdip_BULK_xaxis = []
        self.cosTheta_valdip_NE_xaxis = []

        self.cosTheta_valdip_IF_yaxis = []
        self.cosTheta_valdip_IM_yaxis = []
        self.cosTheta_valdip_BULK_yaxis = []
        self.cosTheta_valdip_NE_yaxis = []

        self.cosTheta_valdip_IF_zaxis = []
        self.cosTheta_valdip_IM_zaxis = []
        self.cosTheta_valdip_BULK_zaxis = []
        self.cosTheta_valdip_NE_zaxis = []

        # --- OH Bond Vectors Orientation Arrays (Degrees) ---
        self.vecOH_IF_xaxis = []
        self.vecOH_IM_xaxis = []
        self.vecOH_BULK_xaxis = []
        self.vecOH_NE_xaxis = []

        self.vecOH_IF_yaxis = []
        self.vecOH_IM_yaxis = []
        self.vecOH_BULK_yaxis = []
        self.vecOH_NE_yaxis = []

        self.vecOH_IF_zaxis = []
        self.vecOH_IM_zaxis = []
        self.vecOH_BULK_zaxis = []
        self.vecOH_NE_zaxis = []

        # --- OH Bond Orientation Arrays (Cosine Theta / Direction Cosines) ---
        self.cosTheta_vecOH_IF_xaxis = []
        self.cosTheta_vecOH_IM_xaxis = []
        self.cosTheta_vecOH_BULK_xaxis = []
        self.cosTheta_vecOH_NE_xaxis = []

        self.cosTheta_vecOH_IF_yaxis = []
        self.cosTheta_vecOH_IM_yaxis = []
        self.cosTheta_vecOH_BULK_yaxis = []
        self.cosTheta_vecOH_NE_yaxis = []

        self.cosTheta_vecOH_IF_zaxis = []
        self.cosTheta_vecOH_IM_zaxis = []
        self.cosTheta_vecOH_BULK_zaxis = []
        self.cosTheta_vecOH_NE_zaxis = []

    def _get_bonded_hydrogens_dist(self, atom):
        """
        Finds covalently bonded hydrogen atoms belonging to the same residue 
        that lie within a strict intramolecular cutoff of the reference atom.
        """
        try:
            # Selects H atoms in the current residue within the intramolecular cutoff distance
            sel_h = atom.residue.atoms.select_atoms(
                "(name H) and around {0:f} index {1!s}".format(
                    self.OH_dist_cutoff, atom.index
                )
            )
            return sel_h
        except NoDataError:
            return []

    def _update_selection_1(self):
        """
        Updates the primary atom groups and builds a dictionary mapping 
        each central atom index to its respective bonded hydrogen atoms.
        """
        self._s1 = self.u.select_atoms(self.selection1)
        self._s1_donors = {}
        self._s1_donors_h = {}
        self._s1_acceptors = {}
        
        # Loop through all selected central atoms (e.g., water Oxygens)
        for i, d in enumerate(self._s1):
            tmp = self._get_bonded_hydrogens_dist(d)
            if tmp:
                self._s1_donors_h[i] = tmp  # Populate dictionary with matched Hydrogens

    def unit_vector(self, vector):
        """Returns the normalized unit vector of a given 3D array vector."""
        return vector / np.linalg.norm(vector)

    def angle_between(self, v1, v2):
        """
        Computes and returns the cosine of the angle (cosTheta) between 
        two vectors using their dot product.
        """
        v1_u = self.unit_vector(v1)
        v2_u = self.unit_vector(v2)
        cosTheta = v1.dot(v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
        return cosTheta

    def vector2degree(self, vec, axis):
        """
        Calculates both the spatial angle in degrees and the cosine value (cosTheta) 
        of a targeted vector relative to one of the primary laboratory axes ('x', 'y', or 'z').
        """
        if axis == "x":
            axis_ndx = 0
            unit_vec = np.array([1, 0, 0])
        elif axis == "y":
            axis_ndx = 1
            unit_vec = np.array([0, 1, 0])
        else:
            axis_ndx = 2
            unit_vec = np.array([0, 0, 1])

        cosTheta_inRadian = self.angle_between(vec, unit_vec)
        # Bound clip to prevent potential floating point errors outside [-1, 1] during arccos
        degree = np.rad2deg(np.arccos(np.clip(cosTheta_inRadian, -1.0, 1.0)))

        return degree, cosTheta_inRadian

    def _getCosTheta(self, ndx):
        """
        Core structural engine:
        1. Manually unwraps coordinates across Periodic Boundary Conditions (PBC) 
           along Y and Z axes to keep water molecules geometrically intact.
        2. Computes the individual OH bond vectors.
        3. Approximates the molecular dipole vector from the bisector.
        4. Profiles orientations into geometric arrays based on X-axis spatial bins.
        """
        d = self._s1[ndx]  # Current central atom (Oxygen)

        dcoordO = []
        xcoordH = []
        ycoordH = []
        zcoordH = []
        tmp_H_position = []
        h_index_tmp_hb = []
        ohb = 0
        
        # Collect raw coordinates of the bonded Hydrogens
        for i in range(len(self._s1_donors_h[ndx])):
            new_H_position = self._s1_donors_h[ndx][i].position[:3]
            h_index_tmp_hb.append(self._s1_donors_h[ndx][i].index)
            xcoordH.append(self._s1_donors_h[ndx][i].position[0])
            ycoordH.append(self._s1_donors_h[ndx][i].position[1])
            zcoordH.append(self._s1_donors_h[ndx][i].position[2])

            # Calculate raw bond distance without simulation box wraps
            ohbond_dist = distances.calc_bonds(
                d.position, self._s1_donors_h[ndx][i].position, box=None
            )  
            if ohbond_dist > 2:  # If distance exceeds a standard bond length, PBC wrap occurred
                ohb = 1

        coordO = d.position
        dcoordO.append(0.0)  # No correction applied to the X-direction
        bcell = self.box[1]         # Box Y length
        ccell = self.box[2]         # Box Z length
        
        # --- Manual Minimum Image Convention / Unwrapping Handling for Oxygen ---
        # Check and correct Y coordinate boundary crossing
        if coordO[1] >= bcell:
            dcoordO.append(m.floor(coordO[1] / bcell) * bcell)
        elif coordO[1] < float(0):
            dcoordO.append(-m.ceil(abs(coordO[1]) / bcell) * bcell)
        else:
            dcoordO.append(0.0)
            
        # Check and correct Z coordinate boundary crossing
        if coordO[2] > ccell:
            dcoordO.append(m.floor(coordO[2] / ccell) * ccell)
        elif coordO[2] <= float(0):
            dcoordO.append(-m.ceil(abs(coordO[2]) / ccell) * ccell)
        else:
            dcoordO.append(0.0)

        # Shift Hydrogen 1 positions relative to the corrected Oxygen coordinate grid
        xcoordHnew = xcoordH[0] - dcoordO[0]
        ycoordHnew = ycoordH[0] - dcoordO[1]
        zcoordHnew = zcoordH[0] - dcoordO[2]
        tmp_H_position.append(np.array([xcoordHnew, ycoordHnew, zcoordHnew]))

        # Shift Hydrogen 2 positions relative to the corrected Oxygen coordinate grid
        xcoordHnew = xcoordH[1] - dcoordO[0]
        ycoordHnew = ycoordH[1] - dcoordO[1]
        zcoordHnew = zcoordH[1] - dcoordO[2]
        tmp_H_position.append(np.array([xcoordHnew, ycoordHnew, zcoordHnew]))

        # Apply spatial corrections back to the actual Oxygen position tracking vector
        coordO[:] = [float(p - q) for (p, q) in zip(coordO, dcoordO)]

        # --- Advanced Hydrogen Unwrapping (If molecule was found broken across PBC) ---
        if ohb == 1:
            for hb_iii in range(len(tmp_H_position)):
                dcoordH_a = []
                ohbond_dist = distances.calc_bonds(
                    coordO, tmp_H_position[hb_iii], box=None
                )  
                if ohbond_dist > 2.0:
                    ohbond_dist_pbc = distances.calc_bonds(
                        coordO, tmp_H_position[hb_iii], box=self.box
                    )  
                    dcoordH_a.append(0.0)
                    bcell = self.box[1]
                    ccell = self.box[2]
                    
                    # Check Hydrogen Y-dimension crossing relative to its shifted coordinates
                    if tmp_H_position[hb_iii][1] >= bcell:
                        dcoordH_a.append(m.floor(tmp_H_position[hb_iii][1] / bcell) * bcell)
                    elif tmp_H_position[hb_iii][1] < float(0):
                        dcoordH_a.append(-m.ceil(abs(tmp_H_position[hb_iii][1]) / bcell) * bcell)
                    else:
                        dcoordH_a.append(0.0)
                        
                    # Check Hydrogen Z-dimension crossing relative to its shifted coordinates
                    if tmp_H_position[hb_iii][2] > ccell:
                        dcoordH_a.append(m.floor(tmp_H_position[hb_iii][2] / ccell) * ccell)
                    elif tmp_H_position[hb_iii][2] <= float(0):
                        dcoordH_a.append(-m.ceil(abs(tmp_H_position[hb_iii][2]) / ccell) * ccell)
                    else:
                        dcoordH_a.append(0.0)

                    # Calculate finalized unwrapped physical coordinates for Hydrogen
                    xcoordHnew = tmp_H_position[hb_iii][0] - dcoordH_a[0]
                    ycoordHnew = tmp_H_position[hb_iii][1] - dcoordH_a[1]
                    zcoordHnew = tmp_H_position[hb_iii][2] - dcoordH_a[2]
                    newH_pos = np.array([xcoordHnew, ycoordHnew, zcoordHnew])
                    tmp_H_position[hb_iii] = newH_pos
                    
                    ohbond_dist_new = distances.calc_bonds(coordO, newH_pos, box=None)
                    delta_ohbond_dist = abs(ohbond_dist_new - ohbond_dist_pbc)
                    
                    # Fine-grain boundary checks for extreme edge cases near half-box boundaries
                    if delta_ohbond_dist > 0.5:
                        ohbond_dist_tmp = distances.calc_bonds(
                            np.array([0, coordO[1], 0]), np.array([0, newH_pos[1], 0]), box=None
                        )
                        arbitary_cutoff = 2.0
                        if (
                            coordO[1] <= self.box[1] / 2
                            and newH_pos[1] > self.box[1] / 2
                            and ohbond_dist_tmp > arbitary_cutoff
                        ):
                            newH_pos[1] = self.box[1] - newH_pos[1]
                        if (
                            coordO[1] >= self.box[1] / 2
                            and newH_pos[1] < self.box[1] / 2
                            and ohbond_dist_tmp > arbitary_cutoff
                        ):
                            newH_pos[1] = self.box[1] + newH_pos[1]
                            
                        # Perform identical fine-grain check along Z-dimension
                        ohbond_dist_tmp = distances.calc_bonds(
                            np.array([0, 0, coordO[2]]), np.array([0, 0, newH_pos[2]]), box=None
                        )
                        if (
                            coordO[2] <= self.box[2] / 2
                            and newH_pos[2] > self.box[2] / 2
                            and ohbond_dist_tmp > arbitary_cutoff
                        ):
                            newH_pos[2] = self.box[2] - newH_pos[2]
                        if (
                            coordO[2] >= self.box[2] / 2
                            and newH_pos[2] < self.box[2] / 2
                            and ohbond_dist_tmp > arbitary_cutoff
                        ):
                            newH_pos[2] = self.box[2] + newH_pos[2]

        # --- Vector Extractions ---
        # Calculate OH bond directional vectors (pointing from Oxygen to Hydrogen)
        OH_Vector_1 = tmp_H_position[0] - coordO  
        OH_Vector_2 = tmp_H_position[1] - coordO  

        # Compute Molecular Dipole Vector approximation: (H1 + H2)/2 - O
        dipVector = (tmp_H_position[0] + tmp_H_position[1]) * 0.5 - coordO

        # --- Compute Angles & Directional Cosines relative to XYZ lab frame ---
        angle_dip_xaxis, cosTheta_dip_xaxis = self.vector2degree(dipVector, "x")
        angle_dip_yaxis, cosTheta_dip_yaxis = self.vector2degree(dipVector, "y")
        angle_dip_zaxis, cosTheta_dip_zaxis = self.vector2degree(dipVector, "z")

        angle_OH_1_xaxis, cosTheta_OH_1_xaxis = self.vector2degree(OH_Vector_1, "x")
        angle_OH_1_yaxis, cosTheta_OH_1_yaxis = self.vector2degree(OH_Vector_1, "y")
        angle_OH_1_zaxis, cosTheta_OH_1_zaxis = self.vector2degree(OH_Vector_1, "z")

        angle_OH_2_xaxis, cosTheta_OH_2_xaxis = self.vector2degree(OH_Vector_2, "x")
        angle_OH_2_yaxis, cosTheta_OH_2_yaxis = self.vector2degree(OH_Vector_2, "y")
        angle_OH_2_zaxis, cosTheta_OH_2_zaxis = self.vector2degree(OH_Vector_2, "z")

        # --- Spatial Slabbing / Sorting Logic based on Oxygen X-Coordinate ---
        
        # Region A: Interface Region (IF)
        if self.cutoff_IF[0] < d.position[0] <= self.cutoff_IF[1]:
            self.valdip_IF_xaxis.append(angle_dip_xaxis)
            self.valdip_IF_yaxis.append(angle_dip_yaxis)
            self.valdip_IF_zaxis.append(angle_dip_zaxis)

            self.cosTheta_valdip_IF_xaxis.append(cosTheta_dip_xaxis)
            self.cosTheta_valdip_IF_yaxis.append(cosTheta_dip_yaxis)
            self.cosTheta_valdip_IF_zaxis.append(cosTheta_dip_zaxis)

            self.vecOH_IF_xaxis.append(angle_OH_1_xaxis)
            self.vecOH_IF_yaxis.append(angle_OH_1_yaxis)
            self.vecOH_IF_zaxis.append(angle_OH_1_zaxis)
            self.vecOH_IF_xaxis.append(angle_OH_2_xaxis)
            self.vecOH_IF_yaxis.append(angle_OH_2_yaxis)
            self.vecOH_IF_zaxis.append(angle_OH_2_zaxis)

            self.cosTheta_vecOH_IF_xaxis.append(cosTheta_OH_1_xaxis)
            self.cosTheta_vecOH_IF_yaxis.append(cosTheta_OH_1_yaxis)
            self.cosTheta_vecOH_IF_zaxis.append(cosTheta_OH_1_zaxis)
            self.cosTheta_vecOH_IF_xaxis.append(cosTheta_OH_2_xaxis)
            self.cosTheta_vecOH_IF_yaxis.append(cosTheta_OH_2_yaxis)
            self.cosTheta_vecOH_IF_zaxis.append(cosTheta_OH_2_zaxis)

        # Region B: Bulk Solvent Region
        elif self.cutoff_BULK[0] < d.position[0] <= self.cutoff_BULK[1]:
            self.cosTheta_valdip_BULK_xaxis.append(cosTheta_dip_xaxis)
            self.cosTheta_vecOH_BULK_xaxis.append(cosTheta_OH_1_xaxis)
            self.cosTheta_vecOH_BULK_xaxis.append(cosTheta_OH_2_xaxis)

        # Region C: Beyond Bulk / Far-end Region (NE)
        elif d.position[0] > self.cutoff_BULK[1]:
            self.cosTheta_valdip_NE_xaxis.append(cosTheta_dip_xaxis)
            self.cosTheta_vecOH_NE_xaxis.append(cosTheta_OH_1_xaxis)
            self.cosTheta_vecOH_NE_xaxis.append(cosTheta_OH_2_xaxis)

    def run(self, start, stop, step):
        """
        Main execution loop over the sliced trajectory frame ranges.
        """
        for ts in self.u.trajectory[start:stop:step]:
            # Step 1: Update target atom selections for current frame
            self._s1 = self.u.select_atoms(self.selection1)

            # Step 2: Map covalent topologies and identify local Hydrogens
            self._update_selection_1()

            # Step 3: Loop over valid mapped solvent molecules to perform orientational checks
            for i, donor_h_set in self._s1_donors_h.items():
                self._getCosTheta(i)

        # Return interface dipole angle distributions for X, Y, and Z axes
        return [self.valdip_IF_xaxis, self.valdip_IF_yaxis, self.valdip_IF_zaxis]
