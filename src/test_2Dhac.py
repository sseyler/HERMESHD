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
fcomm = comm.py2f()       # Get the communicator for Fortran
iam = comm.Get_rank()     # Get the MPI rank of *this* process
nprocs = comm.Get_size()  # Get the total number of MPI ranks
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

    ### LAMMPS #################################################################
    lmp = lammps()
    lmp = setup_md(lmp)
    lmp = setup_wall(lmp)
    lmp = setup_buffer(lmp)
    lmp = setup_md_region(lmp)
    lmp = init_velocities(lmp)
    lmp = equilibrate(lmp, te_init, te_sim)
    ################################
    lmp.command("dump        run all dcd {} {}/md_{}.dcd".format(n_md, datdir, bname))
    lmp.command("dump_modify run pbc yes")
    ################################
    lmp.command("variable dt   equal {}".format(dt_md))
    lmp.command("variable zeta equal {}".format(0.1))
    lmp.command("variable kt   equal {}".format(1.9872036e-3*te_sim))  # gas constant in kcal/mol
    lmp.command("variable std  equal sqrt(2*v_kt*v_zeta/v_dt)")
    lmp.command("variable xc atom x")
    lmp.command("variable yc atom y")

    lmp.command("variable ux python get_Qres_var")
    lmp.command("variable uy equal 0.0")
    lmp.command("python get_Qres_var input 4 v_xc v_yc 0 1 return v_ux format f file lj_hermes_2Dhac.py")
    # lmp.command("variable ux python get_Qres_2Dhac")
    # lmp.command("python get_Qres_2Dhac input 2 v_xc v_yc return v_ux format f file lj_hermes_2Dhac.py")

    lmp.command("variable hfx atom \"-v_zeta*(vx - v_ux) + normal(0, v_std, {})\"".format(seed))
    lmp.command("variable hfy atom \"-v_zeta*(vy - v_uy) + normal(0, v_std, {})\"".format(seed))
    ################################
    lmp.command("fix hforce all addforce v_hfx v_hfy 0 every 10")
    # lmp.command("region  rid_right block {} {} {} {} {} {} units box".format(tb_lo, tb_hi, Lyd, Lyu, Lzd, Lzu))
    # lmp.command("fix hf_right all addforce v_hfx v_hfy 0 every 1 region rid_right")
    # lmp.command("compute grid all chunk/atom bin/3d x lower {} y lower {} z lower {} ids every units box".format(dx,dy,dz))
    # lmp.command("compute ctest all property/chunk grid count coord1 coord2 coord3")
    # lmp.command("fix ftest all ave/chunk {} {} {} grid density/number norm sample ave one file {}".format(neve, nrep, nfre, grid_file))


    ### HERMESHD ###############################################################
    dt_hd = np.array(dt_hd, dtype=float)
    setup_fhd(Qio, t, dt_hd, t1, t_start, dtout, nout, fcomm)
    Qres = Qio[:,:,:,:,0]  # just grab value of constant basis function
    sendbuf_Qres = Qres    # Qres defined in lj_hermes_2Dhac.py
    comm.Gather(sendbuf_Qres, recvbuf_Qres, root=0)  # Gather Qres from each rank


    ### Run ####################################################################
    for step in xrange(nstep_md):
        # LAMMPS
        lmp = run_lammps(lmp, 10, dt=dt_md, nout=100)

        # HERMES
        run_hermeshd(10, Qio, Q1, Q2, t, dt_hd, t1, dtout, nout)
        sendbuf_Qres = Qio[:,:,:,:,0]                    # Get reservoir vars from Qio
        comm.Gather(sendbuf_Qres, recvbuf_Qres, root=0)  # Gather Qres from each rank
        # Set velocity variables using the LAMMPS command
        # lmp.command("fix hforce all addforce v_hfx v_hfy 0 every 1")


    cleanup_fhd(t_start)
    ############################################################################
