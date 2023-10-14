#!/usr/bin/env python3

import elphmod
import numpy as np
import storylines

ph = elphmod.ph.Model('dft/dyn').supercell(3, 3)

if elphmod.MPI.comm.rank != 0:
    raise SystemExit

Margin = 0.9
margin = 0.2

thickness = 0.04

def new_plot(**kwargs):
    plot = storylines.Plot(
        style='NanoLett',

        labelpos='lt',
        labelopt='below right=0.5mm, fill=white, rounded corners',
        labelformat=lambda a: r'\textbf{%s}' % a,

        margin=margin,
        left=Margin,
        bottom=Margin,

        xmin=-40.0,
        xmax=40.0,
        xstep=30.0,
        xminorstep=10.0,
        xminormarks=True,

        **kwargs)

    plot.width = (Margin + 9 * margin - plot.double) / 5
    plot.height = plot.width

    plot.line(grid=True)

    if plot.xlabel is None:
        plot.xlabels = False
        plot.bottom = margin

    if plot.ylabel is None:
        plot.ylabels = False
        plot.left = margin

    return plot

dopings = np.linspace(-0.02, 0.02, 5)

r0 = ph.r
a = ph.a

for l, doping in enumerate(dopings):
    if np.isclose(doping, 0):
        filename = 'model/smearing+0.0050_t1'
    else:
        filename = 'model/doping%+5.2f_t1' % doping

    # panels a to e:

    r = elphmod.ph.read_flfrc('%s_dyn.dat' % filename)[3]

    plot = new_plot(
        label='abcde'[l],

        top=Margin / 2,

        ylabel=None if l else 'CDW structure',
        xyaxes=False,
        )

    plot.lines.pop()

    plot.xmin = None
    plot.xmax = None
    plot.height = 0

    plot.line(*list(zip(0 * a[0], a[0], a[0] + a[1], a[1], 0 * a[1]))[:2],
        draw='none', fill='yellow!10')

    plot.node(*(a[0] / 2)[:2], '$%s\,e$/f.u.'
        % ('%+g' % doping if doping else '0'), above=True)

    tau = np.linalg.norm(r0[1, :2] - r0[0, :2])

    r_sc = np.array([r[i] + m * a[0] + n * a[1]
        for m in [-1, 0, 1]
        for n in [-1, 0, 1]
        for i in range(len(r))])

    bonds = storylines.bonds(R1=r_sc[0::3, :2], R2=r_sc[1::3, :2],
        dmin=0.9 * tau, dmax=1.1 * tau)

    for bond in bonds:
        for bond in storylines.cut(bond, 0, a[1, 1]):
            bond = np.array([elphmod.bravais.rotate(xy, np.pi / 3)
                for xy in bond])

            for bond in storylines.cut(bond, 0, a[1, 1]):
                bond = np.array([elphmod.bravais.rotate(xy, -np.pi / 3)
                    for xy in bond])

                plot.line(*zip(*bond), color='gray')

    u = r - r0

    for i in range(len(u)):
        if i % 3 == 0:
            color = 'black'
        elif i % 3 == 1:
            color = 'gray'
        else:
            continue

        if np.any(abs(u[i, :2]) > 0.035):
            plot.line(*list(zip(r[i], r[i] + 15 * u[i]))[:2],
                thick=True, color=color,
                **{'->': True, 'shorten >': '-1mm', 'shorten <': '1mm'})

    plot.line(*r[0::3, :2].T, mark='ball', mark_size='1mm',
        ball_color='gray', only_marks=True)

    plot.line(*r[1::3, :2].T, mark='ball', mark_size='1mm',
        ball_color='yellow', only_marks=True)

    plot.line(*list(zip(0 * a[0], a[0], a[0] + a[1], a[1], 0 * a[1]))[:2])

    for line in plot.lines:
        line['y'] = -np.array(line['y'])

    plot.save('fig/doping_%s.pdf' % plot.label)

    # panels f to j:

    plot = new_plot(
        label='fghij'[l],

        ymin=0.0,
        ymax=3.0,
        ystep=1.0,
        yminorstep=0.5,
        yminormarks=True,

        ylabel=None if l else 'Density of states (1/eV)',
        )

    w, DOS = np.loadtxt('%s_dos_zoom.dat' % filename).T

    plot.line(1e3 * w, DOS, draw='none', yref=0.0, fill='cyan!50', cut=True)
    plot.line(1e3 * w, DOS, color='cyan', cut=True)

    plot.save('fig/doping_%s.pdf' % plot.label)

    # panels k to o:

    plot = new_plot(
        label='klmno'[l],

        ymin=0.0,
        ymax=0.7 + 1e-10,
        ystep=0.2,
        yminorstep=0.1,
        yminormarks=True,

        xlabel='Energy (meV)',
        ylabel=None if l else r'\hspace*{-7mm}Eliashberg spectral function',

        lpos='rt',
        lopt='below left=1mm, inner sep=2pt',
        llen='2mm',
        lbox=True,
        )

    omega, DOS, a2F = np.loadtxt('%s_a2f.dat' % filename).T

    for sgn in -1, 1:
        plot.line(sgn * omega, DOS, draw='none', yref=0.0, cut=True,
            fill='lightgray')

    plot.axes()

    for sgn in -1, 1:
        plot.line(sgn * omega, a2F, draw='magenta', yref=0.0, cut=True)

    if not l:
        plot.line(color='lightgray', line_width='2mm', line_cap='butt',
            label='DOS (1/meV)')

    plot.save('fig/doping_%s.pdf' % plot.label)

# combine panels:

storylines.combine('fig/doping.png',
    ['doping_%s' % a for a in 'abcdefghijklmno'], columns=5)
