#!/usr/bin/python

import os, sys, getopt
import numpy as np

from mpi4py import MPI

from lammps import lammps
from lj_lammps_2Dhac import *

import hermeshd
from lj_hermes_2Dhac import *

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

nstep_md = 100
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
    lmp.command("variable uy equal 0.001")
    lmp.command("variable uz equal 0.0")
    lmp.command("variable hfx atom \"-v_zeta*(vx - v_ux) + normal(0.0, v_sigma, {})\"".format(seed))
    lmp.command("variable hfy atom \"-v_zeta*(vy - v_uy) + normal(0.0, v_sigma, {})\"".format(seed))
    ################################
    lmp.command("region  rid_right block {} {} {} {} {} {} units box".format(rb_lo, rb_hi, Lyd, Lyu, Lzd, Lzu))
    lmp.command("fix hf_right all addforce v_hfx v_hfy 0 every 1 region rid_right")

    lmp.command("variable xc atom x ")
    lmp.command("variable yc atom y ")

    lmp.command("compute grid all chunk/atom bin/3d x lower {} y lower {} z lower {} ids every units box".format(dx,dy,dz))
    lmp.command("compute ctest all property/chunk grid count coord1 coord2 coord3")
    lmp.command("fix ftest all ave/chunk {} {} {} grid density/number norm sample ave one file {}".format(neve, nrep, nfre, grid_file))

    ### HERMESHD ################################
    dt = np.array(1.0e-2, dtype=float)
    hermeshd.hermeshd.setup(Qio, t, dt, t1, t_start, dtout, nout, fcomm)
    dt_hd = np.array(dt_hd, dtype=float)

    for step in xrange(nstep_md):
        # LAMMPS
        lmp = run_lammps(lmp, 10, dt=dt_md, nout=100)

        # HERMES
        run_hermes(10, Qio, Q1, Q2, t, dt_hd, t1, dtout, nout)
        # Get velocities in reservoir (from Qio)
        # Set velocity variables using the LAMMPS command


    hermeshd.hermeshd.cleanup(t_start)
