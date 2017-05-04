#!/usr/bin/python

import os, sys, getopt
import numpy as np

from mpi4py import MPI

from lammps import lammps
from lj_pylmp_hac import *

import hermeshd
from lj_hermes_hac import *

###########################################
# MPI variables
#------------------------------------------
comm = MPI.COMM_WORLD
iam = comm.Get_rank()
nprocs = comm.Get_size()

# Get the communicator for Fortran
fcomm = MPI.COMM_WORLD.py2f()
#------------------------------------------


################################################################################
# Initialize HAC variables
#-----------------------------
t_md = tstep  # from lj_pylmp_hac
t_hd = 100*t_md


#------------------------------------------------------------------------------


################################################################################
# Run
#-----------------------------
if __name__ == "__main__":

    # LAMMPS
    lmp = lammps()
    lmp = setup(lmp)
    lmp = setup_wall(lmp)
    lmp = equilibrate(lmp)

    # HERMESHD
    setup(Qio, t, dt, t1, t_start, dtout, nout, fcomm)


    while ( t < tf ):
        # LAMMPS
        lmp = run_lammps(lmp, n_md)

        # HERMES
        run_hermes()


        f_hd = get_hd_forces()
        f_md = lmp.gather_atoms("f", 1, 3)
        f_md += f_hd
        lmp.scatter_atoms("f", 1, 3, f_md)


    cleanup(t_start)
