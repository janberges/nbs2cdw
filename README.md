# Charge-density waves in NbS₂

This repository contains the source code and data associated with the paper:

- Timo Knispel, Jan Berges, Arne Schobert, Erik G. C. P. van Loon, Wouter Jolie,
  Tim Wehling, Thomas Michely, and Jeison Fischer,
  *Unconventional charge-density-wave gap in monolayer NbS₂*,
  [arXiv:2307.13791](https://arxiv.org/abs/2307.13791)

We provide all relevant output files of the first-principles calculations and
all data files needed to create the figures. Some intermediate files have been
omitted but can be generated using the provided scripts.

## Software requirements

We have performed all DFT and DFPT calculations using version 7.1 of Quantum
ESPRESSO. Any similarly recent version should work equally well. To interface
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

- `plot.py` (create Figs. 5, S13, S14, and S15)
- `plot_phases.py` (create Fig. S11)
- `plot_doping.py` (create Fig. S12)

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

Data shown in Fig. 4 (T1 phase):

- `smearing+0.0050_t1_hr.dat`
- `scaling+0.0_dos.dat` (also shown in Figs. S13, S14, and S15)
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
- `scaling+0.0_dos_zoom.dat` (also shown in Figs. S13, S14, and S15)
- `scaling+0.1_t1_dos_zoom.dat`
- `scaling+0.2_t1_dos_zoom.dat`
- `scaling+0.3_t1_dos_zoom.dat`
- `scaling+0.4_t1_dos_zoom.dat`
- `scaling+0.5_t1_dos_zoom.dat`
- `scaling+0.6_t1_dos_zoom.dat`
- `scaling+0.7_t1_dos_zoom.dat`
- `scaling+0.8_t1_dos_zoom.dat`
- `scaling+0.9_t1_dos_zoom.dat`
- `smearing+0.0050_t1_dos_zoom.dat` (also shown in Fig. S12)
- `smearing+0.0150_t1_dyn.dat` (to identify amplitude and phase modes)
- `smearing+0.0130_t1_ifc.dat`
- `modes_t1.dat`
- `smearing+0.0050_t1_a2f.dat` (also shown in Fig. S12)

Data shown in Fig. S11 (CDW phases):

- `doping.dat`
- `smearing.dat`

Data shown in Fig. S12 (doping dependence):

- `doping-0.02_t1_dyn.dat` (only atomic positions)
- `doping-0.01_t1_dyn.dat` (only atomic positions)
- `smearing+0.0050_t1_dyn.dat` (only atomic positions)
- `doping+0.01_t1_dyn.dat` (only atomic positions)
- `doping+0.02_t1_dyn.dat` (only atomic positions)
- `doping-0.02_t1_dos_zoom.dat`
- `doping-0.01_t1_dos_zoom.dat`
- `smearing+0.0050_t1_dos_zoom.dat` (also shown in Fig. 4)
- `doping+0.01_t1_dos_zoom.dat`
- `doping+0.02_t1_dos_zoom.dat`
- `doping-0.02_t1_a2f.dat`
- `doping-0.01_t1_a2f.dat`
- `smearing+0.0050_t1_a2f.dat` (also shown in Fig. 4)
- `doping+0.01_t1_a2f.dat`
- `doping+0.02_t1_a2f.dat`

Data shown in Fig. S13 (hexagons phase):

- `smearing+0.0050_hexagons_hr.dat`
- `scaling+0.0_dos.dat` (also shown in Figs. 4, S14, and S15)
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
- `scaling+0.0_dos_zoom.dat` (also shown in Figs. 4, S14, and S15)
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

Data shown in Fig. S14 (T1' phase):

- `smearing+0.0050_t1p_hr.dat`
- `scaling+0.0_dos.dat` (also shown in Figs. 4, S13, and S15)
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
- `scaling+0.0_dos_zoom.dat` (also shown in Figs. 4, S13, and S15)
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

Data shown in Fig. S15 (T2' phase):

- `smearing+0.0050_t2p_hr.dat`
- `scaling+0.0_dos.dat` (also shown in Figs. 4, S13, and S14)
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
- `scaling+0.0_dos_zoom.dat` (also shown in Figs. 4, S13, and S14)
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
