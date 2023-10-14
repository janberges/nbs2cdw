#!/usr/bin/env python3

import elphmod
import numpy as np
import scipy.optimize

comm = elphmod.MPI.comm

dopings = np.linspace(-0.2, 0.2, 41)
phases = 't1', 't2', 't1p', 't2p', 'hexagons', 'stars'

pw = elphmod.bravais.read_pwi('dft/scf.in')

el = elphmod.el.Model('dft/NbS2', rydberg=True)
ph = elphmod.ph.Model('dft/dyn', divide_mass=False, apply_asr_simple=True)
elph = elphmod.elph.Model('dft/NbS2.epmatwp', 'dft/wigner.dat', el, ph,
    divide_mass=False, shared_memory=True)

elph = elph.supercell(3, 3, shared_memory=True)

driver = elphmod.md.Driver(elph, nk=(16, 16), nq=(4, 4), n=9.0,
    kT=pw['degauss'], f=elphmod.occupations.smearing(pw['smearing']))

driver.kT = 0.005

if comm.rank == 0:
    data = open('model/doping.dat', 'w')

for doping in dopings:
    driver.n = 9 * (1 - doping)
    driver.u[:] = 0.0

    E0 = driver.free_energy()

    for phase in phases:
        driver.from_xyz('xyz/%s.xyz' % phase)

        u0 = driver.u / np.linalg.norm(driver.u)

        scipy.optimize.minimize(driver.free_energy, driver.u,
            jac=driver.jacobian, method='BFGS',
            options=dict(gtol=1e-6, norm=np.inf))

        u = driver.u / np.linalg.norm(driver.u)

        E = driver.free_energy()

        C = driver.hessian(fildyn='model/doping%+5.2f_%s_dyn.dat'
            % (doping, phase), apply_asr_simple=True)

        if comm.rank == 0:
            data.write('%5.2f %8s %12.9f %12.9f %12.9f\n' % (doping, phase,
                E - E0, np.dot(u0, u), np.linalg.eigvalsh(C).min()))
            data.flush()

        if phase == 't1' and any(np.isclose(doping, x)
                for x in [-0.02, -0.01, 0.01, 0.02]):

            driver.electrons(seedname='model/doping%+5.2f_%s'
                % (doping, phase), dk1=4, dk2=4)

            driver.phonons(flfrc='model/doping%+5.2f_%s_ifc.dat'
                % (doping, phase), apply_asr_simple=True)

if comm.rank == 0:
    data.close()
