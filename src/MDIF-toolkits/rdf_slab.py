# Copyright (c) 2026 Chanbum Park chanbum.park@theochem.ruhr-uni-bochum.de 
# Distributed under the terms of the GNU General Public License.

"""
Compute Radial Distribution Function (RDF) restricted to a slab geometry.
Optimized for inhomogeneous systems along a specific coordinate axis.
"""

import os
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import distances
import matplotlib.pyplot as pl

class InterSlabRDF:
    def __init__(self, u, box, cutoff_slab, binsize=0.1, exclusion_block=None):
        """
        Parameters:
        -----------
        u : MDAnalysis.Universe
        cutoff_slab : list/tuple
            Boundaries [min, max] along the slicing axis (X-axis here).
        binsize : float
            Width of each RDF bin in Angstroms.
        exclusion_block : tuple, optional
            Intra-molecular exclusion masking pairs (e.g., [3, 3] for water).
        """
        self.u = u
        # Cache initial frame dimensions for settings setup
        self.init_box = self.u.dimensions
        
        self.cutoff_slab = np.asarray(cutoff_slab)
        self.slab_width = self.cutoff_slab[1] - self.cutoff_slab[0]
        
        # Max search radius is half the shortest transverse box length
        max_range = min(self.init_box[1], self.init_box[2]) / 2.0
        nbins = int(max_range / binsize)
        self.rdf_settings = {'bins': nbins, 'range': (0.0, max_range)}
        
        self._exclusion_block = exclusion_block
        
        # Setup bin configurations using a dummy histogram
        _, edges = np.histogram([-1], **self.rdf_settings)
        self.edges = edges
        self.bins = 0.5 * (edges[:-1] + edges[1:])
        
        # Initialize storage arrays
        self._prepare()

    def _prepare(self):
        """Initializes accumulator arrays across the trajectory."""
        self.rdf_slab = np.zeros(self.rdf_settings['bins'], dtype=np.float64)
        self.n_frames = 0 

    def _single_frame(self, ts, g1, g2):
        """Processes a single trajectory frame."""
        # Dynamically fetch current frame box dimension to handle NPT or NVT fluctuations
        box = ts.dimensions
        slab_volume = self.slab_width * self.init_box[1] * self.init_box[2]
        
        # Slicing atoms only located within the specified slab boundaries (along X-axis)
        pos_1_ndx = np.where((g1.positions[:, 0] > self.cutoff_slab[0]) & 
                             (g1.positions[:, 0] <= self.cutoff_slab[1]))[0]
        pos_2_ndx = np.where((g2.positions[:, 0] > self.cutoff_slab[0]) & 
                             (g2.positions[:, 0] <= self.cutoff_slab[1]))[0]
        
        pos_1 = g1.positions[pos_1_ndx]
        pos_2 = g2.positions[pos_2_ndx]
        
        pos_1_N = len(pos_1)
        pos_2_N = len(pos_2)
        
        # If either selection is empty in the slab for this frame, skip distance calculation
        if pos_1_N == 0 or pos_2_N == 0:
            return

        # Distance matrix engine utilizing accelerated cell-linked grids
        pairs, dist = distances.capped_distance(pos_1, pos_2,
                                                self.rdf_settings['range'][1],
                                                box=self.init_box)
        
        # Apply intra-molecular exclusion loops if requested
        if self._exclusion_block is not None:
            idxA = pairs[:, 0] // self._exclusion_block[0]
            idxB = pairs[:, 1] // self._exclusion_block[1]
            mask = np.where(idxA != idxB)[0]
            dist = dist[mask]
            
        # Get raw counts for this frame
        counts, _ = np.histogram(dist, **self.rdf_settings)

        # Ideal number of paired interactions in the subset
        N_pairs = pos_1_N * pos_2_N

        # Cylinder intersection shell volume calculation
        dr = self.edges[1] - self.edges[0]
        vol_slab = 2.0 * np.pi * dr * self.bins * self.slab_width 
        
        # Compute normalized frame density ratios
        density_slab = counts / vol_slab
        density_slab_ideal = N_pairs / slab_volume 

        # Protect against division by zero in empty layers
        if density_slab_ideal > 0:
            rdf_frame = (density_slab / density_slab_ideal)
            self.rdf_slab += rdf_frame
            self.n_frames += 1

    def conclude(self):
        """Averages data over total analyzed frames."""
        if self.n_frames > 0:
            self.rdf_slab_global = self.rdf_slab / self.n_frames
        else:
            self.rdf_slab_global = self.rdf_slab


#if __name__ == '__main__':
#    runName = "run-pos-spce-nvt-2"
#    top_ext, trj_ext = 'pdb', 'dcd'
#    
#    # Check paths locally
#    tpr_file = f"../{runName}.{top_ext}"
#    trj_file = f"../{runName}.{trj_ext}"
#    
#    if not (os.path.exists(tpr_file) and os.path.exists(trj_file)):
#        raise FileNotFoundError("Trajectory or Topology path invalid.")
#
#    u = mda.Universe(tpr_file, trj_file)
#    
#    IdentityA = "OW"
#    IdentityB = "OW"
#    
#    ag1 = u.select_atoms(f'name {IdentityA}')
#    ag2 = u.select_atoms(f'name {IdentityB}')
#    
#    # Initialize our slab analyzer (Slab bounds set from X=0 to X=4)
#    slab_analyzer = InterSlabRDF(u, cutoff_slab=[0, 4], binsize=0.1, exclusion_block=[1, 1])
#    
#    print("Starting trajectory analysis loop...")
#    # Loop across frames without manually clearing internal arrays
#    for ts in u.trajectory:
#        slab_analyzer._single_frame(ts, ag1, ag2)
#        
#    slab_analyzer.conclude()
#    
#    # Resolve your query about the first bin: 
#    # At r -> 0, Lennard-Jones repulsion means atoms cannot overlap. 
#    # Any signal in the first bin is an artifact of bin limits or self-interaction!
#    bins = slab_analyzer.bins
#    rdf_data = slab_analyzer.rdf_slab_global
#    
#    # Save arrays cleanly
#    np.save('rdf_mda_slab_results1.npy', np.array([bins, rdf_data], dtype=object))
#    
#    # Output to clean GitHub-friendly dat file
#    out_name = f"rdf_mda_slab_{IdentityA}-{IdentityB}.dat"
#    with open(out_name, 'w') as fout:
#        for b, r in zip(bins, rdf_data):
#            fout.write(f"{b:.4f} {r:.6f}\n")
#            
#    print(f"Analysis finished. Data saved to {out_name}")
