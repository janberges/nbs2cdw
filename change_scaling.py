#!/usr/bin/env python3

import elphmod
import numpy as np

scalings = np.linspace(0.1, 0.9, 9)
phases = 't1', 't1p', 't2p', 'hexagons'

pw = elphmod.bravais.read_pwi('dft/scf.in')

el = elphmod.el.Model('dft/NbS2', rydberg=True)
ph = elphmod.ph.Model('dft/dyn', divide_mass=False, apply_asr_simple=True)
elph = elphmod.elph.Model('dft/NbS2.epmatwp', 'dft/wigner.dat', el, ph,
    divide_mass=False, shared_memory=True)

driver = elphmod.md.Driver(elph, nk=(48, 48), n=1.0,
    kT=pw['degauss'], f=elphmod.occupations.smearing(pw['smearing']))

driver.kT = 0.005

E0 = 9 * driver.free_energy()

el = driver.electrons(seedname='model/scaling%+4.1f' % 0.0, dk1=4, dk2=4)

elph = elph.supercell(3, 3, shared_memory=True)

driver = elphmod.md.Driver(elph, nk=(16, 16), n=9.0,
    kT=pw['degauss'], f=elphmod.occupations.smearing(pw['smearing']))

driver.kT = 0.005

for scaling in scalings:
    for phase in phases:
        driver.from_xyz('xyz/%s.xyz' % phase)
        driver.u *= scaling

        driver.diagonalize()

        driver.electrons(seedname='model/scaling%+4.1f_%s' % (scaling, phase),
            dk1=4, dk2=4)
