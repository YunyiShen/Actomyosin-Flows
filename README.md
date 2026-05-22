## Code for the numerical experiments used in Fig. 6 of "Long-Range Actomyosin Flows Partition Bulk Cytoplasm After Mitosis in Large Embryos" by Bai et al. 2026+

### Structure and files of the repo:
- `package`
    - `cytoFD`: a FiPy-based numerical Navier–Stokes solver for 2D and 3D spherical cells with active stress. However, 3D simulations are not recommended, as the finite-volume method scales poorly with dimension.
- `script`: code used to conduct the experiments reported in Fig. 6
    - `2Dcell`: simulations in two dimensions
        - `hydrostatic.py`: numerically solves the Stokes equation with various setups for the mechanism proposed in the main paper by Bai et al. 2026+
        - `hydrostatic.sh`: SLURM file for submitting jobs
        - `Fig_6.py`: uses simulation results to generate several subpanels in Fig. 6
        - `Fig_6_video.py`: generates videos of the simulated time series
    - `3Dcell`: simulations in three dimensions. Mathematically, in our setup, 3D simulations can be obtained by rotating the 2D results, but nevertheless we performed simulations directly in 3D.
        - `hydrostatic3D.py`: numerically solves the Stokes equation with various setups for the mechanism proposed in the main paper by Bai et al. 2026+, in 3D
        - `hydrostatic3D.sh`: SLURM file for submitting jobs
        - `interpolating.py`: interpolates lower-resolution results to a slightly higher resolution
        - `interpolating.sh`: SLURM file for submitting jobs
        - `visual.py`: visualizes the 3D simulations


### To repeat the numerical experiments
First install package `cytoFD`, by e.g., `pip install -e ./package`. And get into `script/2Dcell`, run `python hydrostatic.py --run_id 0 --tmax 600 --dt 1 --N 91 --Stokes True` that will simulate Stokes equation up to 600s, with stepsize 1s, and a grid size of 91 by 91, `run_id` 0 corresponding to the setup with stress 1Pa, viscosity 8-40Pa.s and no drag. 

### Run your own experiment
The code is written to be used for quite general active stress schemes beyond the Gaussian bump we used in the paper. The key classes are `cytoFD.forward.planegeneral.ActinModel` and `cytoFD.forward.planegeneral.CellDivFlow2D`. The former can be used to define your own stress model that takes input of space and time that returns biology, namely viscosity, stress and drag. The latter is a general solver class that used to define solver states like grid size and similar things. 

For example, to simulate our Gaussian bump model one can do
```
from cytoFD.forward.planegeneral import growinggaussianbump2Dconc
cell_radius = 0.5
tmax = 600
stress_max = 1e3

actin = growinggaussianbump2Dconc(theta = (90/90) *  np.pi/2,
                                precision = np.array([[25, 0],[0, 30**2]]) / ((cell_radius * 2) ** 2),
                                timescale = tmax/2.)

```

Or we can have a surface tension version that is implemented in the package

```
from cytoFD.forward.planegeneral import growinggaussianbump2Dconc, cortex2Dconc
actin = cortex2Dconc(R = 0.95*cell_radius, width = 0.02,  # depth of the surface stress
                             aspect_ratio = 1, 
                             phase_rot = 0., timescale = tmax/2)

```
Now we make the biology class and numerically solve the equation using the solver
```
### define biology
biology = ActinModel(actin = actin,
                        stress_range = [1e-5, stress_max],
                        drag_range = [0, 0],
                        visc_range = [8000, 40000],
                        domain_size = 2 * cell_radius, # 1mm 
                        cell_radius = cell_radius
                        )

### define solver
myflow = CellDivFlow2D(N=int(2 * cell_radius/dx),
                        Stokes = True      
                        )

### solve biology 
res = myflow.solve(biology, dt = dt, steps = int(tmax/dt), 
                    save_every = max(int(tmax/dt/100), 1))

```