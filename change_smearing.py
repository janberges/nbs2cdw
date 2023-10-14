#!/usr/bin/env python3

import elphmod
import numpy as np
import scipy.optimize

comm = elphmod.MPI.comm

smearings = np.linspace(0.002, 0.020, 37)
phases = 't1', 't2', 't1p', 't2p', 'hexagons', 'stars'
selection = dict(t1=0.013, t1p=0.0115, t2p=0.007, hexagons=0.013)

def get_nk(smearing):
    return 4 * int(np.ceil(0.02 / smearing)) # 4 originally

pw = elphmod.bravais.read_pwi('dft/scf.in')

el = elphmod.el.Model('dft/NbS2', rydberg=True)
ph = elphmod.ph.Model('dft/dyn', divide_mass=False, apply_asr_simple=True)
elph = elphmod.elph.Model('dft/NbS2.epmatwp', 'dft/wigner.dat', el, ph,
    divide_mass=False, shared_memory=True)

elph = elph.supercell(3, 3, shared_memory=True)

if comm.rank == 0:
    data = open('model/smearing.dat', 'w')

nk = 0

for smearing in smearings:
    if nk != get_nk(smearing):
        nk = get_nk(smearing)

        driver = elphmod.md.Driver(elph, nk=(nk, nk),
            nq=(1, 1) if smearing < 0.005 else (4, 4), n=9.0,
            kT=pw['degauss'], f=elphmod.occupations.smearing(pw['smearing']))

    driver.kT = smearing
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

        C = driver.hessian(fildyn='model/smearing%+7.4f_%s_dyn.dat'
            % (smearing, phase), apply_asr_simple=True)

        if comm.rank == 0:
            data.write('%6.4f %8s %12.9f %12.9f %12.9f\n' % (smearing, phase,
                E - E0, np.dot(u0, u), np.linalg.eigvalsh(C).min()))
            data.flush()

        if phase in selection:
            selected = np.isclose(smearing, 0.005)

            if selected:
                driver.electrons(seedname='model/smearing%+7.4f_%s'
                    % (smearing, phase), dk1=nk // 4, dk2=nk // 4)

            if selected or np.isclose(smearing, selection[phase]):
                driver.phonons(flfrc='model/smearing%+7.4f_%s_ifc.dat'
                    % (smearing, phase), apply_asr_simple=True)

if comm.rank == 0:
    data.close()
