#!/bin/bash
#SBATCH --nodes 1
#SBATCH --tasks-per-node 96
#SBATCH --partition standard96
#SBATCH --time 12:00:00

echo 'Probably this script will not work out of the box on your system.'
echo 'Please only run selected parts manually, e.g., just the plotting.'
exit

module load anaconda3 texlive intel impi

export SLURM_CPU_BIND=none

cd $SLURM_SUBMIT_DIR

nc=${SLURM_CPUS_ON_NODE:-16}
np=16
nt=$((nc/np))

cd dft

echo 'Using normconserving pseudopotentials from PseudoDojo'
echo '[1] van Setten et al., Comput. Phys. Commun. 226, 39 (2018)'
echo '[2] Hamann, Phys. Rev. B 88, 085117 (2013)'

url=http://www.pseudo-dojo.org/pseudos/nc-sr-04_pbe_standard

for pp in Nb.upf S.upf
do
    test -e $pp || (wget $url/$pp.gz && gunzip $pp)
done

source elphmodenv 1

mpirun pw.x -nk $np < scf.in > scf.out
mpirun ph.x -nk $np < ph.in > ph.out

ph2epw

mpirun pw.x -nk $np < nscf.in > nscf.out
mpirun epw.x -nk $nc < epw.in > epw.out

cd ..

mkdir -p xyz

source elphmodenv $nt

for phase in t1 t2 t1p t2p hexagons stars
do
    mpirun -n $np python3 cdw_relax.py $phase
done

mkdir -p model

mpirun -n $np python3 change_smearing.py
mpirun -n $np python3 change_doping.py
mpirun -n $np python3 change_scaling.py

source elphmodenv 1

for zoom in 0 1
do
    mpirun python3 compute_dos.py "model/scaling+0.0" $zoom 0

    for phase in t1 t1p t2p hexagons
    do
        for scaling in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9
        do
            mpirun python3 compute_dos.py "model/scaling+"$scaling"_"$phase $zoom 0
        done

        mpirun python3 compute_dos.py "model/smearing+0.0050_"$phase $zoom 0
    done
done

for phase in t1 t1p t2p hexagons
do
    mpirun python3 compute_a2f.py "model/smearing+0.0050_"$phase 0
done

for doping in -0.02 -0.01 +0.01 +0.02
do
    mpirun python3 compute_dos.py "model/doping"$doping"_t1" 1 0
    mpirun python3 compute_a2f.py "model/doping"$doping"_t1" 0
done

source elphmodenv $nc

for phase in t1 t1p t2p hexagons
do
    python3 compute_modes.py $phase
done

mkdir -p fig

for phase in t1 t1p t2p hexagons
do
    python3 plot.py $phase
done

python3 plot_doping.py
python3 plot_phases.py
