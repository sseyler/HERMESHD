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
dt_md = 10.0    # from lj_pylmp_hac  in fs
dt_hd = 1.0e-2  # In ps

nstep_md = 10000
nstep_hd = n_md/(dt_hd*1e3/dt_md)
#------------------------------------------------------------------------------


################################################################################
# Run
#-----------------------------
if __name__ == "__main__":

    ### LAMMPS ################################
    lmp = lammps()
    lmp = setup(lmp)
    lmp = setup_wall(lmp)
    lmp = setup_buffer(lmp)
    lmp = setup_md_region(lmp)
    lmp = init_velocities(lmp)
    lmp = equilibrate(lmp, te_init, te_sim)
    ################################
    lmp.command("dump     run all dcd {} {}/md_{}.dcd".format(n_md, datdir, bname))
    lmp.command("dump_modify run pbc yes")
    ################################
    lmp.command("variable dt equal {}".format(dt_md))
    lmp.command("variable zeta  equal {}".format(0.1))
    lmp.command("variable kt    equal {}".format(1.9872036e-3*te_sim))  # gas constant in kcal/mol
    lmp.command("variable sigma equal sqrt(2*v_kt*v_zeta/v_dt)")
    lmp.command("variable ux equal 0.0")
    lmp.command("variable uy equal 0.01")
    lmp.command("variable uz equal 0.0")
    lmp.command("variable hfx atom \"-v_zeta*(vx - v_ux) + normal(0.0, v_sigma, {})\"".format(seed))
    lmp.command("variable hfy atom \"-v_zeta*(vy - v_uy) + normal(0.0, v_sigma, {})\"".format(seed))
    lmp.command("variable hfz atom \"-v_zeta*(vz - v_uz) + normal(0.0, v_sigma, {})\"".format(seed))
    ################################
    lmp.command("region  rid_right block {} {} {} {} {} {} units box".format(rb_lo, rb_hi, Lyd, Lyu, Lzd, Lzu))
    lmp.command("fix hf_right all addforce v_hfx v_hfy v_hfz every 1 region rid_right")

    ### HERMESHD ################################
    dt = np.array(0.0, dtype=float)  # temp var to initialize HERMESHD
    setup(Qio, t, dt, t1, t_start, dtout, nout, fcomm)
    dt_hd = np.array(dt_hd, dtype=float)

    for step in xrange(nstep_md):
        # LAMMPS
        lmp = run_lammps(lmp, 1, dt=dt_md, nout=100)

        # HERMES
        run_hermes(1, Qio, Q1, Q2, t, dt_hd, dtout, nout)
        # Get velocities in buffer from Qio
        # Set velocity variables using the LAMMPS command


    cleanup(t_start)
