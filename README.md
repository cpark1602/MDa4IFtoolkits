# MDIF-toolkits
MDIF-toolkits is a analysis package for atomistic simulations with interfacial geometry.

# MDIF-toolkits

[![CI](https://github.com/cpark1602/MDIF-toolkits/actions/workflows/ci.yml/badge.svg)](https://github.com/cpark1602/MDIF-toolkits/actions)
[![Code Style: Ruff](https://img.shields.io/badge/code%20style-ruff-000000.svg)](https://github.com/astral-sh/ruff)

A specialized scientific Python package designed to post-process, analyze, and extract structural and transport properties from Molecular Dynamics (MD) trajectories of complex electrified interfaces and slab geometries.

---

## Installation

### 1. Requirements
* **Python** $\ge$ 3.9
* Core dependencies (automatically resolved during installation): `numpy`, `scipy`, `pandas`

### 2. Local Installation (Development Mode)
To install the toolkit locally on your laptop or office machine so that changes you make to the source files are instantly active, clone the repository and run an editable installation:

```bash
git clone [https://github.com/cpark1602/MDIF-toolkits.git](https://github.com/cpark1602/MDIF-toolkits.git)
cd MDIF-toolkits
pip install -e .
```

---


## Available Analysis Modules

The package contains targeted workflows optimized for spatial sorting and time-correlation analysis of interfacial systems:

### Mass Density Profiling

Calculates the spatial mass distribution of different atomic species (e.g., water, ions) along a specific axis perpendicular to the slab geometry surface.


#### Usage
```bash
import MDAnalysis as mda
import scipy.constants

u_if = mda.Universe("run-pos.pdb", "run-pos.dcd")
boxX = 48.57
boxY = 15.667
boxZ = 15.076
box = [boxX, boxY, boxZ, 90, 90, 90]
u_if.dimensions = box
start_stop_step = [0, -1, 1]
print_results_path = "/results/"

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
```


### Water Dipole Orientation Angles
Analyzes the structural ordering and polarization of water molecules in proximity to electrified interfaces by measuring the distribution of dipole vector angles relative to the surface normal.

#### Usage

```bash
u_if = mda.Universe("run-pos.pdb", "run-pos.dcd")
boxX = 48.57
boxY = 15.667
boxZ = 15.076
box = [boxX, boxY, boxZ, 90, 90, 90]
u_if.dimensions = box
start_stop_step = [0, -1, 1]
print_results_path = "/results/"

# ----- Load dipole angle analysis
import dipole_angles
pbc=True
dim='x'
bin_size = 0.02
selection1 = 'name O'
_u_if = dipole_angles.Dangling_bonds(u_if, box, print_results_path, pbc, bin_size, dim, selection1, cutoff_IF = [0, 12], cutoff_BULK = [19, 28], start=start_stop_step[0], stop=start_stop_step[1], step=start_stop_step[2])

valdip_IF =  _u_if.run( start_stop_step[0], start_stop_step[1], start_stop_step[2] )
```

### Hydrogen Bond Autocorrelation Function (ACF)
Tracks the dynamic lifetime, breaking, and forming kinetics of hydrogen bond networks near a surface using structural time-correlation frameworks.


### Ionic Conductivity
Evaluates charge transport dynamics across the electrolyte layer by parsing collective ionic current fluctuations under bias potentials.

### Mean Squared Displacement (MSD)

Calculates the translational diffusion coefficients of ionic species and solvent clusters in bulk or confined slab environments.

### Radial Distribution Function (RDF) in Slab Geometry
Computes a modified, spatially restricted g(r) function to evaluate pair correlation properties inside local thin-film slabs without bulk volume projection artifacts.

```bash
import rdf_slab
IdentityA = "O"
IdentityB = "O"

ag1 = u_if.select_atoms(f'name {IdentityA}')
ag2 = u_if.select_atoms(f'name {IdentityB}')

# Initialize our slab analyzer (Slab bounds set from X=0 to X=4)
slab_analyzer = rdf_slab.InterSlabRDF(u_if, box, cutoff_slab=[8, 12], binsize=0.1, exclusion_block=[1, 1])

print("Starting trajectory analysis loop...")
# Loop across frames without manually clearing internal arrays
for ts in u_if.trajectory:
    slab_analyzer._single_frame(ts, ag1, ag2)

slab_analyzer.conclude()

bins = slab_analyzer.bins
rdf_data = slab_analyzer.rdf_slab_global
```


## Kirkwood $G_K$ Factor

The local orientational order of molecular dipole moments is intrinsically linked to the macroscopic dielectric constant of polar liquids. This relationship can be expressed via the Kirkwood-Fröhlich equation:

$$\frac{4\pi \beta N \mu^2 G_K}{\Omega} = \frac{(\epsilon - 1)(2\epsilon + 1)}{\epsilon}$$

Where:
* $N$ is the number of polar molecules in a system of total volume $\Omega$.
* $\mu$ is the magnitude of the liquid-phase molecular dipole moment.
* $\beta$ represents the thermodynamic beta ($\frac{1}{k_B T}$).
* $\epsilon$ is the static dielectric constant of the system.

The distance-dependent Kirkwood $G_K$ factor approaches a constant value at long ranges, making it an exceptionally useful property to estimate the dielectric constant from short-range molecular interactions.

#### Mathematical Formulation
The distance-dependent Kirkwood $G_K$ factor evaluated for water configurations is computed using the following ensemble average:

$$G_K(r) = \frac{\langle \sum_{i=1}^{N} \mathbf{\mu}_i \cdot \mathbf{M}_i(r) \rangle}{\langle N \rangle \mu^2}$$

Where $\mathbf{M}_i(r)$ is the net total sum of all molecular dipoles $\mathbf{\mu}$ located within a cutoff sphere of radius $r$ centered around the reference dipole $\mathbf{\mu}_i$ (explicitly including $\mathbf{\mu}_i$ itself).

#### Usage

```bash
u_if = mda.Universe("run-pos.pdb", "run-pos.dcd")
boxX = 48.57
boxY = 15.667
boxZ = 15.076
box = [boxX, boxY, boxZ, 90, 90, 90]
u_if.dimensions = box
start_stop_step = [0, -1, 1]
print_results_path = "/results/"

import kirkwood_gk_interface

pbc=True
dim='x'
bin_size = 0.02
selection1 = 'name O'; selection2 = 'name O'
if_q0_nac = kirkwood_gk_interface.Kirkwood_Gk(u_if, box, print_results_path, pbc, bin_size, dim, selection1, selection2, cutoff_IF = [0, 12], cutoff_BULK = [19, 28], start=start_stop_step[0], stop=start_stop_step[1], step=start_stop_step[2])
```

### Ionic conductivity


#### Usage
```bash
import conductivity

pbc=True
Ex = 0.01 # External Field
atom1 = 'Li' # Lithium ion
atom2 = 'N3' # TFSI ion

u_if_q0_nac = mda.Universe(os.path.join(w_path,trj_file_pdb), os.path.join(w_path,trj_file_trj_1))

start_stop_step = [0, -1, 1] 
if_q0_nac = conductivity.Conductivity(u_if_q0_nac, 'name '+atom1, 'name '+atom2, print_results_path, Ex, start=start_stop_step[0], stop=start_stop_step[1], step=start_stop_step[2])
if_q0_nac.run()
kw_gk_mu_aver_global = if_q0_nac.run()
```

### Mean Squared Displacement

#### Usage
```bash
import msd

trj_file1 = "run-pos.pdb"
trj_file2 = "run-pos.dcd"

u = mda.Universe(os.path.join(trj_file1), os.path.join(trj_file2))
print("total nr. of frame: ", len(u.trajectory))
tot_frames = len(u.trajectory)
skip = int(1)
start_stop_step = [0, -1, skip]  # if xtc gro are loaded.

u_msd = msd.MSD(
    u,
    select="index 2",
    msd_type="xyz",
    fft=True,
    start=start_stop_step[0],
    stop=start_stop_step[1],
    step=start_stop_step[2],
)  # tot_frames

u_msd.run()
msd = u_msd.timeseries
nframes = u_msd.n_frames
timestep = 1  # 0.5
lagtimes = np.arange(nframes) * timestep * skip
```



---

## Automated Code Quality & Testing

This project uses a modern CI/CD pipeline powered by GitHub Actions:

    Linting / Formatting: Regulated via ruff inside pre-commit hooks to guarantee clean code style.

    Unit Testing: Powered by pytest across Python 3.9, 3.10, and 3.11 matrices.

To run tests locally on your machine, simply execute:



---

## Contributors & Contact

    Developer: Chanbum Park (Sorbonne Université, CNRS)
               chanbum.park@sorbonne-universite.fr
               chanbum.park@ruhr-uni-bochum.de
               chanbum.park@theochem.ruhr-uni-bochum.de

    GitHub Repository: cpark1602/MDIF-toolkits


