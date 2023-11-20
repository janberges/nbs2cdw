[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10144912.svg)](https://doi.org/10.5281/zenodo.10144912)

# Charge-density waves in NbS₂

This repository contains the data and source code associated with the paper:

- Timo Knispel, Jan Berges, Arne Schobert, Erik G. C. P. van Loon, Wouter Jolie,
  Tim Wehling, Thomas Michely, and Jeison Fischer,
  *Unconventional charge-density-wave gap in monolayer NbS₂*,
  [arXiv:2307.13791](https://arxiv.org/abs/2307.13791)

## Experimental setup

STM and STS were carried out at a base operating temperature of T₀ = 0.35 K
after in-situ transfer from the preparation chamber. STS was performed with the
lock-in technique. STM images were taken in constant current mode.

## Computational setup

All DFT and DFPT calculations were performed using Quantum ESPRESSO 7.1.
Similarly recent versions should work equally well. (Note that the ab initio
step can be skipped since we provide all relevant output files.) To interface
the EPW code to our *elphmod* package, the file `EPW/src/ephwann_shuffle.f90`
has to be modified. After the determination of the Wigner-Seitz points, just
before the Bloch-to-Wannier transform, add the following lines:

    IF (ionode) THEN
      OPEN (13, file='wigner.dat', action='write', status='replace', &
        access='stream')
      WRITE (13) dims, dims2
      WRITE (13) nrr_k, irvec_k, ndegen_k
      WRITE (13) nrr_g, irvec_g, ndegen_g
      CLOSE (13)
    ENDIF

Our *elphmod* and *StoryLines* packages and all other Python requirements can be
installed in a virtual environment:

    python3 -m venv venv
    source venv/bin/activate
    python3 -m pip install -r requirements.txt

A LaTeX installation, preferably TeX Live, is required to typeset the figure.

## List of scripts

Workflow of the DFT, DFPT, and model calculations:

- `run.sh`

Initial structural relaxation:

- `cdw_relax.py` (create contents of directory `xyz`)

Variation of parameters:

- `change_smearing.py` (rerelax CDW structures for different smearings)
- `change_doping.py` (rerelax CDW structures for different dopings)
- `change_scaling.py` (rescale CDW displacements)

Characterization of selected CDW phases:

- `compute_dos.py` (calculate electron density of states)
- `compute_a2f.py` (calculate phonon density of states and Eliashberg function)
- `compute_modes.py` (calculate and identify amplitude and phase modes)

Plotting:

- `plot.py` (create Figures 5, S13, S14, and S15)
- `plot_phases.py` (create Figure S11)
- `plot_doping.py` (create Figure S12)

## Contents of directory `exp`

Experimental data shown in the main text:

- `Figure 1a.txt`
- `Figure 1b.txt`
- `Figure 2a.txt`
- `Figure 2a inset.txt`
- `Figure 2b.txt`
- `Figure 2b inset.txt`
- `Figure 2c.txt`
- `Figure 2c inset.txt`
- `Figure 2d.txt`
- `Figure 2d left inset.txt`
- `Figure 2d right inset.txt`
- `Figure 2e.txt`
- `Figure 2e inset.txt`
- `Figure 2f.txt`
- `Figure 2f inset.txt`
- `Figure 2g.txt`
- `Figure 2g inset.txt`
- `Figure 2h.txt`
- `Figure 3a.txt`
- `Figure 3b.txt`
- `Figure 3c.txt`

In Figure 1a, the size was reduced from 250 nm × 250 nm to 250 nm × 166 nm. This
reduction is not included in the source data.

In Figure 1b, the minimum height on the graphene surface is corrected to 0 nm.

In the insets of Figure 2a–c and Figure 2e–g, the images were clipped with a
400% zoom.

In the insets of Figure 2d, the images were clipped as a circle in the same
location.

Figure 2h shows the FFT intensity of first-order 3 × 3 CDW spots
(I<sub>CDW</sub>) vs. first-order Bragg spots (I<sub>atom</sub>) at different
temperatures. The spot intensity was extracted by line profiles along the six
high symmetry directions in the FFT. For each temperature, the average ratio
I<sub>CDW</sub>/I<sub>atom</sub> of all six directions was taken. FFTs were
extracted from constant current dI/dV maps as shown in Figure 2e–g.

In Figure 3a, the size was reduced from 3.5 nm × 3.5 nm to 2.7 nm × 2.7 nm. This
reduction is not included in the source data.

In Figure 3b, the data was plotted with an arbitrary offset between the spectra
at different positions, for clarity. This offset is not included in the source
data.

In Figure 3c, the data was plotted with an offset of 10 nS between the spectra,
for clarity. This offset is not included in the source data.

## Contents of directory `dft`

Quantum ESPRESSO input files:

- `scf.in` (self-consistent DFT calculation with `pw.x`)
- `nscf.in` (non-self-consistent DFT calculation with `pw.x`)
- `ph.in` (DFPT calculation with `ph.x`)
- `epw.in` (Wannierization with `epw.x`)

Relevant Quantum ESPRESSO output files:

- `NbS2_hr.dat` (hopping parameters)
- `NbS2_wsvec.dat` (superlattice vectors to correct Wigner-Seitz cell)
- `dyn0` ... `dyn19` (force constants)
- `NbS2.epmatwp` (electron-phonon coupling)
- `wigner.dat` (lattice vectors of electron-phonon coupling)

## Contents of directory `xyz`

Relaxed CDW structures:

- `t1.xyz`
- `t2.xyz` (at doping of -0.1 e/f.u.)
- `t1p.xyz`
- `t2p.xyz`
- `hexagons.xyz`
- `stars.xyz` (at doping of -0.1 e/f.u.)

## Contents of directory `model`

Data shown in Figure 4 (T1 phase):

- `smearing+0.0050_t1_hr.dat`
- `scaling+0.0_dos.dat` (also shown in Figures S13, S14, and S15)
- `scaling+0.1_t1_dos.dat`
- `scaling+0.2_t1_dos.dat`
- `scaling+0.3_t1_dos.dat`
- `scaling+0.4_t1_dos.dat`
- `scaling+0.5_t1_dos.dat`
- `scaling+0.6_t1_dos.dat`
- `scaling+0.7_t1_dos.dat`
- `scaling+0.8_t1_dos.dat`
- `scaling+0.9_t1_dos.dat`
- `smearing+0.0050_t1_dos.dat`
- `scaling+0.0_dos_zoom.dat` (also shown in Figures S13, S14, and S15)
- `scaling+0.1_t1_dos_zoom.dat`
- `scaling+0.2_t1_dos_zoom.dat`
- `scaling+0.3_t1_dos_zoom.dat`
- `scaling+0.4_t1_dos_zoom.dat`
- `scaling+0.5_t1_dos_zoom.dat`
- `scaling+0.6_t1_dos_zoom.dat`
- `scaling+0.7_t1_dos_zoom.dat`
- `scaling+0.8_t1_dos_zoom.dat`
- `scaling+0.9_t1_dos_zoom.dat`
- `smearing+0.0050_t1_dos_zoom.dat` (also shown in Figure S12)
- `smearing+0.0150_t1_dyn.dat` (to identify amplitude and phase modes)
- `smearing+0.0130_t1_ifc.dat`
- `modes_t1.dat`
- `smearing+0.0050_t1_a2f.dat` (also shown in Figure S12)

Data shown in Figure S11 (CDW phases):

- `doping.dat`
- `smearing.dat`

Data shown in Figure S12 (doping dependence):

- `doping-0.02_t1_dyn.dat` (only atomic positions)
- `doping-0.01_t1_dyn.dat` (only atomic positions)
- `smearing+0.0050_t1_dyn.dat` (only atomic positions)
- `doping+0.01_t1_dyn.dat` (only atomic positions)
- `doping+0.02_t1_dyn.dat` (only atomic positions)
- `doping-0.02_t1_dos_zoom.dat`
- `doping-0.01_t1_dos_zoom.dat`
- `smearing+0.0050_t1_dos_zoom.dat` (also shown in Figure 4)
- `doping+0.01_t1_dos_zoom.dat`
- `doping+0.02_t1_dos_zoom.dat`
- `doping-0.02_t1_a2f.dat`
- `doping-0.01_t1_a2f.dat`
- `smearing+0.0050_t1_a2f.dat` (also shown in Figure 4)
- `doping+0.01_t1_a2f.dat`
- `doping+0.02_t1_a2f.dat`

Data shown in Figure S13 (hexagons phase):

- `smearing+0.0050_hexagons_hr.dat`
- `scaling+0.0_dos.dat` (also shown in Figures 4, S14, and S15)
- `scaling+0.1_hexagons_dos.dat`
- `scaling+0.2_hexagons_dos.dat`
- `scaling+0.3_hexagons_dos.dat`
- `scaling+0.4_hexagons_dos.dat`
- `scaling+0.5_hexagons_dos.dat`
- `scaling+0.6_hexagons_dos.dat`
- `scaling+0.7_hexagons_dos.dat`
- `scaling+0.8_hexagons_dos.dat`
- `scaling+0.9_hexagons_dos.dat`
- `smearing+0.0050_hexagons_dos.dat`
- `scaling+0.0_dos_zoom.dat` (also shown in Figures 4, S14, and S15)
- `scaling+0.1_hexagons_dos_zoom.dat`
- `scaling+0.2_hexagons_dos_zoom.dat`
- `scaling+0.3_hexagons_dos_zoom.dat`
- `scaling+0.4_hexagons_dos_zoom.dat`
- `scaling+0.5_hexagons_dos_zoom.dat`
- `scaling+0.6_hexagons_dos_zoom.dat`
- `scaling+0.7_hexagons_dos_zoom.dat`
- `scaling+0.8_hexagons_dos_zoom.dat`
- `scaling+0.9_hexagons_dos_zoom.dat`
- `smearing+0.0050_hexagons_dos_zoom.dat`
- `smearing+0.0150_hexagons_dyn.dat` (to identify amplitude and phase modes)
- `smearing+0.0130_hexagons_ifc.dat`
- `modes_hexagons.dat`
- `smearing+0.0050_hexagons_a2f.dat`

Data shown in Figure S14 (T1' phase):

- `smearing+0.0050_t1p_hr.dat`
- `scaling+0.0_dos.dat` (also shown in Figures 4, S13, and S15)
- `scaling+0.1_t1p_dos.dat`
- `scaling+0.2_t1p_dos.dat`
- `scaling+0.3_t1p_dos.dat`
- `scaling+0.4_t1p_dos.dat`
- `scaling+0.5_t1p_dos.dat`
- `scaling+0.6_t1p_dos.dat`
- `scaling+0.7_t1p_dos.dat`
- `scaling+0.8_t1p_dos.dat`
- `scaling+0.9_t1p_dos.dat`
- `smearing+0.0050_t1p_dos.dat`
- `scaling+0.0_dos_zoom.dat` (also shown in Figures 4, S13, and S15)
- `scaling+0.1_t1p_dos_zoom.dat`
- `scaling+0.2_t1p_dos_zoom.dat`
- `scaling+0.3_t1p_dos_zoom.dat`
- `scaling+0.4_t1p_dos_zoom.dat`
- `scaling+0.5_t1p_dos_zoom.dat`
- `scaling+0.6_t1p_dos_zoom.dat`
- `scaling+0.7_t1p_dos_zoom.dat`
- `scaling+0.8_t1p_dos_zoom.dat`
- `scaling+0.9_t1p_dos_zoom.dat`
- `smearing+0.0050_t1p_dos_zoom.dat`
- `smearing+0.0150_t1p_dyn.dat` (to identify amplitude and phase modes)
- `smearing+0.0115_t1p_ifc.dat`
- `modes_t1p.dat`
- `smearing+0.0050_t1p_a2f.dat`

Data shown in Figure S15 (T2' phase):

- `smearing+0.0050_t2p_hr.dat`
- `scaling+0.0_dos.dat` (also shown in Figures 4, S13, and S14)
- `scaling+0.1_t2p_dos.dat`
- `scaling+0.2_t2p_dos.dat`
- `scaling+0.3_t2p_dos.dat`
- `scaling+0.4_t2p_dos.dat`
- `scaling+0.5_t2p_dos.dat`
- `scaling+0.6_t2p_dos.dat`
- `scaling+0.7_t2p_dos.dat`
- `scaling+0.8_t2p_dos.dat`
- `scaling+0.9_t2p_dos.dat`
- `smearing+0.0050_t2p_dos.dat`
- `scaling+0.0_dos_zoom.dat` (also shown in Figures 4, S13, and S14)
- `scaling+0.1_t2p_dos_zoom.dat`
- `scaling+0.2_t2p_dos_zoom.dat`
- `scaling+0.3_t2p_dos_zoom.dat`
- `scaling+0.4_t2p_dos_zoom.dat`
- `scaling+0.5_t2p_dos_zoom.dat`
- `scaling+0.6_t2p_dos_zoom.dat`
- `scaling+0.7_t2p_dos_zoom.dat`
- `scaling+0.8_t2p_dos_zoom.dat`
- `scaling+0.9_t2p_dos_zoom.dat`
- `smearing+0.0050_t2p_dos_zoom.dat`
- `smearing+0.0150_t2p_dyn.dat` (to identify amplitude and phase modes)
- `smearing+0.0070_t2p_ifc.dat`
- `modes_t2p.dat`
- `smearing+0.0050_t2p_a2f.dat`

Relevant force-constants data (not needed for plotting):

- `smearing+0.0050_t1_ifc.dat`
- `smearing+0.0050_hexagons_ifc.dat`
- `smearing+0.0050_t1p_ifc.dat`
- `smearing+0.0050_t2p_ifc.dat`
