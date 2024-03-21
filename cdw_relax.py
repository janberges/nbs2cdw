#!/usr/bin/env python3

import elphmod
import numpy as np
import scipy.optimize
import sys

phase = sys.argv[1] if len(sys.argv) > 1 else 'random'

pw = elphmod.bravais.read_pwi('dft/scf.in')

el = elphmod.el.Model('dft/NbS2', rydberg=True)
ph = elphmod.ph.Model('dft/dyn', divide_mass=False, apply_asr_simple=True)
elph = elphmod.elph.Model('dft/NbS2.epmatwp', 'dft/wigner.dat', el, ph,
    divide_mass=False, shared_memory=True)

elph = elph.supercell(3, 3, shared_memory=True)

driver = elphmod.md.Driver(elph, nk=(16, 16), n=9.0,
    kT=pw['degauss'], f=elphmod.occupations.smearing(pw['smearing']))

driver.kT = 0.005

guess = dict(
    t1=[(2, [0, 18, 24]), (1, [9, 12, 21]), (1, [3, 6, 15])],
    t2=[(1, [0, 18, 24]), (-1, [9, 12, 21]), (2, [3, 6, 15])],
    t1p=[(2, [0, 6, 24]), (1, [3, 12, 15]), (1, [9, 18, 21])],
    t2p=[(1, [0, 6, 24]), (-1, [3, 12, 15]), (2, [9, 18, 21])],
    hexagons=[(2, [0, 3, 9, 15, 21, 24])],
    stars=[(2, [0, 3]), (2, [9, 21]), (2, [15, 24])],
    )

if phase in guess:
    if phase in {'t2', 'stars'}:
        driver.n = 9.9

    for scale, atoms in guess[phase]:
        center = np.average(driver.elph.ph.r[atoms], axis=0)

        for atom in atoms:
            u = center - driver.elph.ph.r[atom]
            u *= 0.05 * scale / np.linalg.norm(u)

            driver.u[3 * atom:3 * atom + 3] = u
else:
    driver.random_displacements()

driver.plot(label=False, interactive=True)

scipy.optimize.minimize(driver.free_energy, driver.u, jac=driver.jacobian,
    method='BFGS', options=dict(gtol=1e-6, norm=np.inf))

driver.to_xyz('xyz/%s.xyz' % phase)
