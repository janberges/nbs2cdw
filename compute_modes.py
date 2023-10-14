#!/usr/bin/env python3

import elphmod
import numpy as np
import sys

if elphmod.MPI.comm.rank != 0:
    raise SystemExit

phase = sys.argv[1] if len(sys.argv) > 1 else 't1'

thr_olap = 0.95
thr_freq = -1e-5

((q,), (D,)), amass, at, r0 = elphmod.ph.read_flfrc(
    'model/smearing+0.0150_%s_dyn.dat' % phase)[:4]

elphmod.ph.divide_by_mass(D, amass)

w02, u0 = np.linalg.eigh(D)
w0 = elphmod.ph.sgnsqrt(w02)

smearings = []
energies = []
modes = []

with open('model/smearing.dat') as lines:
    for line in lines:
        smearing, cdw, energy, overlap, freq = line.split()

        if cdw != phase:
            continue

        smearing = float(smearing)
        overlap = float(overlap)
        freq = float(freq)

        if smearing < 0.015 and (overlap <= thr_olap or freq <= thr_freq):
            continue

        ((q,), (D,)), amass, at, r = elphmod.ph.read_flfrc(
            'model/smearing%+7.4f_%s_dyn.dat' % (smearing, phase))[:4]

        elphmod.ph.divide_by_mass(D, amass)

        w2, u = np.linalg.eigh(D)
        w = elphmod.ph.sgnsqrt(w2) * 1e3 * elphmod.misc.Ry
        u = u.real

        v = (u[:, np.newaxis, :] * u0[:, 3:9, np.newaxis]).sum(axis=0).real

        weight = (abs(v) ** 2).sum(axis=0)

        w = w[weight > 0.5]
        u = u[:, weight > 0.5]
        v = v[:, weight > 0.5]

        U = (r - r0).ravel()
        U *= np.sqrt(np.repeat(amass, 3))
        U /= np.linalg.norm(U) or 1.0

        weight_Higgs = abs((u * U[:, np.newaxis]).sum(axis=0)) ** 2

        nu_Higgs = np.argmax(weight_Higgs)

        w[[-1, nu_Higgs]] = w[[nu_Higgs, -1]]
        v[:, [-1, nu_Higgs]] = v[:, [nu_Higgs, -1]]

        smearings.append(smearing)
        energies.append(w)
        modes.append(v)

smearings = np.array(smearings)
energies = np.array(energies)
modes = np.array(modes)

order = elphmod.dispersion.band_order(energies[..., :-1], modes[..., :-1])

with open('model/modes_%s.dat' % phase, 'w') as data:
    for i in range(len(smearings)):
        data.write('%6.4f' % smearings[i])

        for nu in range(len(energies[i]) - 1):
            data.write(' %6.3f' % energies[i, order[i, nu]])

        data.write(' %6.3f\n' % energies[i, -1])
