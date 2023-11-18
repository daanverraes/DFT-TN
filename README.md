# Bridging the gap between density functional theory and quantum tensor networks to accurately model strongly correlated nanostructured materials

## Prelude
**PyWannier90**: (HPC: 6.4.0 GCCore 10.3.0 IPython 7.25.0) \
module load pyWannier90/2021-12-07-foss-2021a\
module load matplotlib/3.4.2-foss-2021a \

**ASE**: (HPC: 6.4.0 GCCore 11.3.0 IPython 8.5.0) \
module load PySCF/2.1.1-foss-2022a \
module load ASE/3.22.1-foss-2022a \

**Local**: \
PySCF 2.4.0


## To Do
- Resolve issue with isosurfaces
- Disentanglement: hopping, pywannier90, cRPA, ... (benchmark with La2CuO4)
- cRPA
- cGW (self-energy corrections)
- speed-up code/calc vs. accuracy: symmetriseren, parallelliseren, smearing, cutoff, ...
- Alpha en beta (vb. KUKS -> W90: Welke golffuncties?) + spin: KRKS, KUKS + fix errors/warnings?
- heatmap: lm-decomp op bands
- init_guess: fixen?
- DMRG, DMET, DMFT, ... (paper Gagliardi)
