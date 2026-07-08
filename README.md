# MDa4IFtoolkits
**MDa4IFtoolkits** (Molecular Dynamics Analysis for Interfacial Systems Toolkits) is an open-source analysis package designed for atomistic simulations featuring interfacial geometries.

<!--- [![CI](https://github.com/cpark1602/MDa4IFtoolkits/actions/workflows/ci.yml/badge.svg)](https://github.com/cpark1602/MDa4IFtoolkits/actions) --->
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
git clone [https://github.com/cpark1602/MDa4IFtoolkits.git](https://github.com/cpark1602/MDa4IFtoolkits.git)
cd MDa4IFtoolkits
pip install -e .
```

---


## Available Analysis Modules

The package contains targeted workflows optimized for spatial sorting and time-correlation analysis of interfacial systems:

### Mass Density Profiling

Calculates the spatial mass distribution of different atomic species (e.g., water, ions) along a specific axis perpendicular to the slab geometry surface.

![The mass density profiles across the Au/water/Ne simulation cell (source: doi.org/10.1103/7ht1-xv82).](./figs/density.png)

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
_u_if = mass_density.Mass_density(
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

his_edges, number_density_O = _u_if._get_densityProfile("O")
his_edges, number_density_H = _u_if._get_densityProfile("H")
his_edges, number_density_Au = _u_if._get_densityProfile("Au")
his_edges, number_density_Ne = _u_if._get_densityProfile("Ne")
```


### Water Dipole Orientation Angles
Analyzes the structural ordering and polarization of water molecules in proximity to electrified interfaces by measuring the distribution of dipole vector angles relative to the surface normal.

![Snapshot of the water dipole angle orientation (source: doi.org/10.1103/7ht1-xv82).](./figs/water-dipole-angle.png)

#### Usage

```bash
import MDAnalysis as mda
u_if = mda.Universe("run-pos.pdb", "run-pos.dcd")
boxX = 48.57
boxY = 15.667
boxZ = 15.076
box = [boxX, boxY, boxZ, 90, 90, 90]
u_if.dimensions = box
start_stop_step = [0, -1, 1]
print_results_path = "/results/"

import dipole_angles
pbc=True
dim='x'
bin_size = 0.02
selection1 = 'name O'
_u_if = dipole_angles.Dangling_bonds(u_if, box, print_results_path, pbc, bin_size, dim, selection1, cutoff_IF = [0, 12], cutoff_BULK = [19, 28], start=start_stop_step[0], stop=start_stop_step[1], step=start_stop_step[2])

valdip_IF =  _u_if.run( start_stop_step[0], start_stop_step[1], start_stop_step[2] )
```

### Hydrogen Bond Autocorrelation Function (ACF)
Liquid water molecules form a hydrogen-bonded network characterized by continuous breaking and reformation due to thermal fluctuations. These dynamics are intrinsically linked to both the translational and rotational motion of the molecules.

To quantify these processes in our molecular dynamics (MD) simulations, we analyze the second-rank orientational time correlation function (OTCF), $C_2(t)$, of the $\text{O–H}$ bond:

$$C_2(t) = \frac{\langle P_2(\mathbf{u}(0) \cdot \mathbf{u}(t)) \rangle}{\langle P_2(\mathbf{u}(0) \cdot \mathbf{u}(0)) \rangle}$$,

where $P_2(x)$ is the second-rank Legendre polynomial defined as:

$$P_2(x) = \frac{3x^2 - 1}{2}$$,

where $\mathbf{u}(t)$ represents the unit vector along the $\text{O–H}$ bond axis at time $t$. The dot product $\mathbf{u}(0) \cdot \mathbf{u}(t) = \cos\theta(t)$ defines the cosine of the angle $\theta$ swept by the bond vector between the initial time $0$ and time $t$. The angled brackets $\langle \dots \rangle$ denote an ensemble average over multiple time origins and all active water molecules in the system.

The orientational time correlation function $C_2(t)$ is particularly useful because it is the function directly related to NMR relaxation rates. From this, the reorientational correlation time ($\tau_2$) can be obtained by integrating $C_2(t)$ over time:

$$\tau_2 = \int_{0}^{\infty} C_2(t) \, dt$$.


#### Usage
```bash
import MDAnalysis as mda
u_if = mda.Universe("run-pos.pdb", "run-pos.dcd")
boxX = 48.57
boxY = 15.667
boxZ = 15.076
box = [boxX, boxY, boxZ, 90, 90, 90]
u_if.dimensions = box
start_stop_step = [0, -1, 1]
print_results_path = "/results/"

import acf
_u_if = acf.ACF(u_if, box, HBs_criteria_input, 'name O', 'name O', print_results_path, cutoff_dist_O_H =3.5, cutoff_dist_donor_acceptor = 3.5, cutoff_IF = [0, 12], cutoff_BULK = [19, 28], angle=35.0, pbc=True, start=start_stop_step[0], stop=start_stop_step[1], step=start_stop_step[2], nac='IF') 

_u_if.run()
```


### Ionic Conductivity
Evaluates charge transport dynamics across the electrolyte layer by parsing collective ionic current fluctuations under bias potentials.

The electric conductivities of the solutions in our simulations are computed using Ohm’s law:

$$J = \kappa E$$

where $E$ is the applied electric field and $J = \sum_i J_i = \sum_i \kappa_i E$ represents the total current density. This relationship naturally defines the partial ionic conductivity, $\kappa_i$, for each individual ionic species. 

The individual conductivities ($\kappa_i$) are evaluated from the linear slope of $J$ versus $E$ strictly within the linear-response region.

### Radial Distribution Function (RDF) in Slab Geometry
When computing the Radial Distribution Function (RDF) near an interface, the system exhibits strong anisotropy along the surface normal (typically chosen as the $z$-axis or $x$-axis depending on your simulation setup). Because of this broken symmetry, the conventional isotropic 3D spherical RDF fails to properly capture the local structural changes.

To account for this non-spherical symmetry, the RDF within a slab geometry is evaluated by constraining the reference and target particles inside an inhomogeneous lateral slice (an $xy$-slab) of a fixed total width, $2h$ (or $\Delta x$).

The corresponding volume element $V(r)$ representing the intersection of a spherical shell of radius $r \rightarrow r + dr$ with the slab boundary is given by a cylindrical shell approximation:



$$V(r) = 4\pi h r dr$$

The local number density $\rho(r)$ within this localized shell volume is calculated as:

$$\rho(r) = \frac{N(r)}{V(r)}$$

Where $N(r)$ is the average number of particles found within the shell.

Finally, the normalized RDF within the slab geometry, $g_{\text{slab}}(r)$, is computed by scaling the local shell density by the average bulk density of the target slab region ($\rho_{\text{slab}}$):

$$g_{\text{slab}}(r) = \frac{\rho(r)}{\rho_{\text{slab}}}$$.


```bash
import MDAnalysis as mda
import rdf_slab
IdentityA = "O"
IdentityB = "O"

ag1 = u_if.select_atoms(f'name {IdentityA}')
ag2 = u_if.select_atoms(f'name {IdentityB}')

slab_analyzer = rdf_slab.InterSlabRDF(u_if, box, cutoff_slab=[8, 12], binsize=0.1, exclusion_block=[1, 1])

print("Starting trajectory analysis loop...")

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
import MDAnalysis as mda
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
_u_if = kirkwood_gk_interface.Kirkwood_Gk(u_if, box, print_results_path, pbc, bin_size, dim, selection1, selection2, cutoff_IF = [0, 12], cutoff_BULK = [19, 28], start=start_stop_step[0], stop=start_stop_step[1], step=start_stop_step[2])
```

### Ionic conductivity


#### Usage
```bash
import MDAnalysis as mda
import conductivity

pbc=True
Ex = 0.01 # External Field
atom1 = 'Li' # Lithium ion
atom2 = 'N3' # TFSI ion

u_if = mda.Universe(os.path.join(w_path,trj_file_pdb), os.path.join(w_path,trj_file_trj_1))

start_stop_step = [0, -1, 1] 
_u_if = conductivity.Conductivity(u_if, 'name '+atom1, 'name '+atom2, print_results_path, Ex, start=start_stop_step[0], stop=start_stop_step[1], step=start_stop_step[2])
_u_if.run()
kw_gk_mu_aver_global = _u_if.run()
```

### Mean Squared Displacement
This repository provides Python tools to compute the Mean Squared Displacement (MSD) of molecules from Molecular Dynamics (MD) trajectory data using both direct tracking and accelerated Fast Fourier Transform (FFT) methods.

#### Theoretical Framework

Mean Squared Displacement (MSD) is a standard statistical measure in physics and chemistry that quantifies how far a particle (or molecule) moves from its starting position over a given time interval ($t$). It tracks the spatial exploration of an object over time and serves as a fundamental metric to compute the self-diffusion coefficient ($D$) of fluids.

Mathematically, the MSD is defined as an ensemble average over all particles and time origins:

$$\text{MSD}(t) = \langle |\mathbf{r}_i(t) - \mathbf{r}_i(0)|^2 \rangle$$

Where $\mathbf{r}_i(t)$ represents the position vector of particle $i$ at time $t$, and the angled brackets $\langle \dots \rangle$ denote an average over all active particles and multiple time origins to maximize statistical sampling.

#### Calculation of Diffusion Coefficients
In the long-time limit (the linear diffusive regime), the diffusion coefficient $D$ can be extracted from the slope of the MSD curve via the **Einstein relation**:

$$D = \lim_{t \to \infty} \frac{1}{2d \cdot t} \langle |\mathbf{r}_i(t) - \mathbf{r}_i(0)|^2 \rangle$$

Where $d$ represents the dimensionality of the system (typically $d = 3$ for standard bulk 3D diffusion).

#### Usage
```bash
import MDAnalysis as mda
import msd

trj_file1 = "run-pos.pdb"
trj_file2 = "run-pos.dcd"

u_if = mda.Universe(os.path.join(trj_file1), os.path.join(trj_file2))
print("total nr. of frame: ", len(u.trajectory))
tot_frames = len(u_if.trajectory)
skip = int(1)
start_stop_step = [0, -1, skip]  # if xtc gro are loaded.

_u_if = msd.MSD(
    u_if,
    select="index 2",
    msd_type="xyz",
    fft=True,
    start=start_stop_step[0],
    stop=start_stop_step[1],
    step=start_stop_step[2],
)  

_u_if.run()
msd = _u_if.timeseries
nframes = _u_if.n_frames
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

    GitHub Repository: cpark1602/MDa4IFtoolkits


