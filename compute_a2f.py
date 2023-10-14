#!/usr/bin/env python3

import elphmod
import numpy as np
import sys

comm = elphmod.MPI.comm

seedname = sys.argv[1] if len(sys.argv) > 1 else 'model/smearing+0.0050_t1'
fast = bool(int(sys.argv[2])) if len(sys.argv) > 2 else True

nk = nq = 6 if fast else 72
f = elphmod.occupations.gauss
kTel = 100.0 if fast else 20.0
kTph = 1.0 if fast else 0.2

wmin = 0.0
wmax = 40.0
nw = 101 if fast else 1001

eps = 1e-10

el = elphmod.el.Model(seedname)
el.data *= 1e3
ph = elphmod.ph.Model('%s_ifc.dat' % seedname, apply_asr_simple=True)
ph.data *= (1e3 * elphmod.misc.Ry) ** 2

el0 = elphmod.el.Model('dft/NbS2')
ph0 = elphmod.ph.Model('dft/dyn', apply_asr_simple=True)
ph0.data *= (1e3 * elphmod.misc.Ry) ** 2
elph0 = elphmod.elph.Model('dft/NbS2.epmatwp', 'dft/wigner.dat', el0, ph0)
elph0.data *= (1e3 * elphmod.misc.Ry) ** 1.5

elph = elph0.supercell(3, 3, shared_memory=True)

q = sorted(elphmod.bravais.irreducibles(nq))
weights = np.array([len(elphmod.bravais.images(q1, q2, nq)) for q1, q2 in q])
q = 2 * np.pi / nq * np.array(q)

e, U = elphmod.dispersion.dispersion_full_nosym(el.H, nk, vectors=True)
e = e[..., :9]
U = U[..., :9]

w2, u = elphmod.dispersion.dispersion(ph.D, q, vectors=True)
dangerous = np.where(w2 < eps)
w2[dangerous] = eps
w = elphmod.ph.sgnsqrt(w2)

d2 = elph.sample(q, U=U, u=u, squared=True, shared_memory=True)

omega, domega = np.linspace(wmin, wmax, nw, retstep=True)

sizes, bounds = elphmod.MPI.distribute(len(omega), bounds=True)

my_DOS = np.empty(sizes[comm.rank])
my_a2F = np.empty(sizes[comm.rank])

DOS = np.empty(len(omega))
a2F = np.empty(len(omega))

g2dd, dd = elphmod.diagrams.double_fermi_surface_average(q, e, d2, kTel, f)
g2dd[dangerous] = 0.0
g2dd /= 2 * w

for my_n, n in enumerate(range(*bounds[comm.rank:comm.rank + 2])):
    x = (w - omega[n]) / kTph

    my_DOS[my_n] = (weights[:, np.newaxis] * f.delta(x) / kTph).sum()
    my_a2F[my_n] = (weights[:, np.newaxis] * f.delta(x) / kTph * g2dd).sum()

my_DOS /= weights.sum()
my_a2F /= (weights * dd).sum()

N0 = f.delta(e / kTel).sum() / kTel / np.prod(e.shape[:-1])

my_DOS /= len(elph.cells)
my_a2F *= N0

comm.Allgatherv(my_DOS, (DOS, sizes))
comm.Allgatherv(my_a2F, (a2F, sizes))

if comm.rank == 0:
    with open('%s_a2f.dat' % seedname, 'w') as data:
        for iw in range(nw):
            data.write('%5.2f %9.6f %9.6f\n' % (omega[iw], DOS[iw], a2F[iw]))
