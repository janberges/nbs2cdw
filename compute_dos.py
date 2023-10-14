#!/usr/bin/env python3

import elphmod
import numpy as np
import sys

comm = elphmod.MPI.comm

seedname = sys.argv[1] if len(sys.argv) > 1 else 'model/smearing+0.0050_t1'
zoom = bool(int(sys.argv[2])) if len(sys.argv) > 2 else False
fast = bool(int(sys.argv[3])) if len(sys.argv) > 3 else True

nk = 40 if fast else 400

if zoom:
    wmin = -0.04
    wmax = +0.04
    nw = 81 if fast else 801
else:
    wmin = -0.40
    wmax = +1.00
    nw = 141 if fast else 1401

w = np.linspace(wmin, wmax, nw)

el = elphmod.el.Model(seedname)

e = elphmod.dispersion.dispersion_full_nosym(el.H, nk, shared_memory=True)

DOS = 0

status = elphmod.misc.StatusBar(el.size // 3, title='calculate DOS')

for n in range(el.size // 3):
    DOS += elphmod.dos.hexDOS(e[:, :, n], minimum=wmin, maximum=wmax)(w)

    status.update()

DOS /= el.size // 3

if comm.rank == 0:
    with open('%s_dos%s.dat' % (seedname, '_zoom' * zoom), 'w') as data:
        for iw in range(nw):
            data.write(('%7.4f' if zoom else '%6.3f') % w[iw])
            data.write(' %9.6f\n' % DOS[iw])
