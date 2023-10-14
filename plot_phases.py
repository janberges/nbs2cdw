#!/usr/bin/env python3

import elphmod
import numpy as np
import storylines

ph = elphmod.ph.Model('dft/dyn').supercell(3, 3)

if elphmod.MPI.comm.rank != 0:
    raise SystemExit

Margin = 0.9
margin = 0.2

def new_plot(**kwargs):
    plot = storylines.Plot(
        style='NanoLett',

        ymax=0.0,
        ymin=-4.0,
        ystep=1.0,
        yminorstep=0.25,
        yminormarks=True,

        labelpos='lt',
        labelopt='below right=0.5mm, fill=white, rounded corners',
        labelformat=lambda a: r'\textbf{%s}' % a,

        left=Margin,
        bottom=Margin,
        margin=margin,

        **kwargs)

    plot.width = (Margin + 5 * margin - plot.double) / 3
    plot.height = plot.width

    plot.line(grid=True)

    if plot.ylabel is None:
        plot.ylabels = False
        plot.left = margin

    return plot

label = dict(fill='white', rounded_corners='1pt', inner_sep='1pt',
    above_left='0.5mm')

# panel a:

r0 = ph.r
a = ph.a

for l, title in enumerate(['T1', "T1$'$", 'hexagons', 'T2', "T2$'$", 'stars']):
    xyz = 'xyz/%s.xyz' % title.replace("$'$", 'p').lower()

    atm = np.loadtxt(xyz, skiprows=2, dtype=str, usecols=0)
    r = np.loadtxt(xyz, skiprows=2, dtype=float, usecols=(1, 2, 3))

    plot = storylines.Plot(
        style='NanoLett',
        height=0,
        margin=margin,
        xyaxes=False,
        )

    plot.width = (2 * margin - new_plot().width) / 3

    plot.line(*list(zip(0 * a[0], a[0], a[0] + a[1], a[1], 0 * a[1]))[:2],
        draw='none', fill='white')

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

    plot.line(*list(zip(0 * a[0], a[0], a[0] + a[1], a[1], 0 * a[1]))[:2],
        color='gray')

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

    plot.node(0.5 * a[1, 1], -0.5 * a[1, 0], title, above=True, rotate=30)

    for line in plot.lines:
        line['x'], line['y'] = line['y'], -np.array(line['x'])

    plot.save('fig/phases_a%d.pdf' % (l + 1))

plot = new_plot(label='a', xyaxes=False)

plot.lines = []

plot.xmin = 0.0
plot.xmax = 6.0
plot.ymin = 0.0
plot.ymax = 4.0

plot.node(1, 2.7, r'\includegraphics{phases_a1}')
plot.node(3, 2.7, r'\includegraphics{phases_a2}')
plot.node(5, 2.7, r'\includegraphics{phases_a3}')
plot.node(1, 0.7, r'\includegraphics{phases_a4}')
plot.node(3, 0.7, r'\includegraphics{phases_a5}')
plot.node(5, 0.7, r'\includegraphics{phases_a6}')

plot.save('fig/phases_a.pdf')

# panel b:

style = dict(
    t1=dict(color='blue', mark='triangle*'),
    t2=dict(color='orange', mark='square*'),
    t1p=dict(color='magenta', mark='triangle*', mark_options='{rotate=180}'),
    t2p=dict(color='brown', mark='square'),
    hexagons=dict(color='green', mark='square*', mark_options='{rotate=45}'),
    stars=dict(color='red', mark='*'),
    )

thr_olap = 0.95
thr_freq = -1e-5

plot = new_plot(
    label='b',

    xmin=-0.175 - 1e-10,
    xmax=0.175 + 1e-10,
    xstep=0.1,
    xminorstep=0.025,
    xminormarks=True,

    xlabel='Doping ($e$/f.u.)',
    ylabel='Energy change (meV/f.u.)',

    lpos='ct',
    lopt='below=1mm, xshift=-4mm',
    lbox=True,
    lwid=2.5,
    lcol=3,
    )

plot.width -= Margin

data = {key: [] for key in style}

with open('model/doping.dat') as lines:
    for line in lines:
        doping, phase, energy, overlap, freq = line.split()

        if float(overlap) > thr_olap and float(freq) > thr_freq:
            data[phase].append([float(doping), float(energy)])

data = {key: np.array(value) for key, value in data.items()}

for phase in data:
    points = np.array(data[phase])

    if points.size:
        plot.line(points[:, 0], points[:, 1] * 1e3 * elphmod.misc.Ry / 9,
            cut=True, label=phase.replace('t1', 'T1').replace('t2', 'T2')
            .replace('p', r"$\smash'$"), **style[phase])

plot.node(plot.xmax, plot.ymin, r'5\,mRy', **label)

plot.save('fig/phases_b.pdf')

# panel c:

plot = new_plot(
    label='c',

    xmax=14.5,
    xstep=2,
    xminorstep=0.5,
    xminormarks=True,

    xlabel='Smearing (mRy)',
    )

plot.axes()

plot.width += Margin

data = {key: [] for key in style}

with open('model/smearing.dat') as lines:
    for line in lines:
        smearing, phase, energy, overlap, freq = line.split()

        if float(overlap) > thr_olap and float(freq) > thr_freq:
            data[phase].append([float(smearing), float(energy)])

for phase in data:
    points = np.array(data[phase])

    if points.size:
        plot.line(points[:, 0] * 1e3, points[:, 1] * 1e3 * elphmod.misc.Ry / 9,
            zindex=2, cut=True, **style[phase])

plot.node(plot.xmax, plot.ymin, r'0\,$e$/f.u.', **label)

plot.save('fig/phases_c.pdf')

# combine panels:

storylines.combine('fig/phases.png', ['phases_%s' % a for a in 'abc'],
    columns=4)
