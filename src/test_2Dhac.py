#!/usr/bin/python
from __future__ import print_function, division
from mpi4py import MPI

import os, sys, getopt
import numpy as np
import io, tempfile
from contextlib import contextmanager

from lammps import lammps, PyLammps
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
dxi = (nx*mpi_nx)/Lx
dyi = (ny*mpi_ny)/Ly
dzi = (nz*mpi_nz)/Lz
dx = 1./dxi
dy = 1./dyi
dz = 1./dzi

dt_md = 1.0    # from lj_pylmp_hac  in fs
dt_hd = 1.0e-2  # In ps

nstep_md = 100
nstep_hd = n_md/(dt_hd*1e3/dt_md)
#------------------------------------------------------------------------------

def get_global_idx(x, y, z):
    gid_i = int(np.ceil((x - Lxd) / dx))
    gid_j = int(np.ceil((y - Lyd) / dy))
    gid_k = int(np.ceil((z - Lzd) / dz))
    return gid_i, gid_j, gid_k

def get_mpi_idx_from_gid(gid_i, gid_j, gid_k):
    mpi_i = int( np.ceil(gid_i / nx) )
    mpi_j = int( np.ceil(gid_j / ny) )
    mpi_k = int( np.ceil(gid_k / nz) )
    return mpi_i, mpi_j, mpi_k

def get_local_idx(gid_i, gid_j, gid_k, mpi_i, mpi_j, mpi_k):
    i = int( gid_i - (mpi_i-1)*nx )
    j = int( gid_j - (mpi_j-1)*ny )
    k = int( gid_k - (mpi_k-1)*nz )
    return i, j, k

def get_mpi_rank(mpi_i, mpi_j, mpi_k):
    return mpi_ny*mpi_nx*(mpi_k-1) + mpi_nx*(mpi_j-1) + (mpi_i-1)

def xyc_to_ux(crd2D, rbuf_Qres):
    x, y = tuple(crd2D)  # assumes numpy array or list
    gid_i, gid_j, gid_k = get_global_idx(x, y, 0)
    mpi_i, mpi_j, mpi_k = get_mpi_idx_from_gid(gid_i, gid_j, gid_k)
    i, j, k = get_local_idx(gid_i, gid_j, gid_k, mpi_i, mpi_j, mpi_k)
    rank = get_mpi_rank(mpi_i, mpi_j, mpi_k)
    # print('x',x,'y',y,'mpi',mpi_i,mpi_j,mpi_k,'lid',i,j,k)
    # print(rho, m_x)
    return 0.01 * rbuf_Qres[rank,i-1,j-1,mx]/rbuf_Qres[rank,i-1,j-1,rh]  # WARN: should be 0.001

def ux_to_atomfile(natoms, rbuf_Qres, fname='ux_per_atom.dat'):
    atompos_lmp = lmp.gather_atoms("x", 1, 2)  # NOTE: Must be executed by all MPI ranks
    if iam == 0:
        if os.path.isfile(fname):
            os.unlink(fname)
            assert not os.path.exists(fname)
        tf = tempfile.NamedTemporaryFile(dir='./tmp', delete=False)
        # tf.name = fname

        xyc = np.fromiter(atompos_lmp, dtype=np.float, count=2*natoms).reshape((natoms,2)) # copy
        ux_strlist = ["{} {}\n".format(i+1,xyc_to_ux(xyc[i,:],rbuf_Qres)) for i in xrange(natoms)]
        # print("\n>>>> Here's the ux value:",xyc_to_ux(xyc[0,:],rbuf_Qres))
        tf.writelines(ux_strlist)
        tf.flush()
        fnameout = tf.name
    else:
        fnameout = None
    fnameout = comm.bcast(fnameout, root=0)
    comm.Barrier()
    return fnameout

def test_write_atomfile(natoms, fname='ux_per_atom.dat'):
    if iam == 0:
        # if os.path.isfile(fname):
        #     os.unlink(fname)
        #     assert not os.path.exists(fname)
        tf = tempfile.NamedTemporaryFile(dir='.', delete=False)
        # tf.name = fname
        strlst = ["{} {}\n".format(i+1,"<some value>") for i in xrange(natoms)]
        tf.writelines(strlst)
        tf.flush()
        fnameout = tf.name
    else:
        fnameout = None
    fnameout = comm.bcast(fnameout, root=0)
    comm.Barrier()
    return fnameout

################################################################################
# Run
#-----------------------------
if __name__ == "__main__":
    '''
    lmp.command("variable ux python lmp_py_func")
    lmp.command("python lmp_py_func input 2 0.0 0.0 return v_ux format fff file lj_hermes_2Dhac.py")
    ERROR: Cannot embed Python when also extending Python with LAMMPS (../python.cpp:156)
    Last command: python lmp_py_func input 2 0.0 0.0 return v_ux format fff file lj_hermes_2Dhac.py
    '''
    ### LAMMPS #################################################################
    lmp = lammps()  # L = PyLammps(ptr=lmp)
    setup_md(lmp)
    setup_wall(lmp)
    setup_buffer(lmp)
    setup_md_region(lmp)
    init_velocities(lmp, te_init)
    minimize(lmp, 100)
    equilibrate(lmp, te_init, te_sim, fast=False)
    #-------------------------------
    #-------------------------------
    lmp.command("variable dt   equal {}".format(dt_md))
    lmp.command("variable zeta equal {}".format(1.0))
    lmp.command("variable kt   equal {}*temp".format(1.9872036e-3))  # gas constant in kcal/mol
    lmp.command("variable std  equal sqrt(2*v_kt*v_zeta/v_dt)")
    #-------------------------------
    natoms = lmp.get_natoms()
    # atompos_lmp = lmp.gather_atoms("x", 1, 2)  # WARN: Must be executed by all MPI ranks
    # xc = None
    # L.variable("ux equal 0.1") WTF, broken: lmp.set_variable("ux", 0.1)

    ### HERMESHD ###############################################################
    Qres_shape = (nx, ny, 2)
    rbuf_Qres_shape = (nproc, nx, ny, 2)
    sbuf_Qres = np.empty(Qres_shape, dtype=np.float32)
    if iam == 0:
        rbuf_Qres = np.empty(rbuf_Qres_shape, dtype=np.float32)
    else:
        rbuf_Qres = None

    dt_hd = np.array(dt_hd, dtype=float)
    setup_fhd(Qio, t, dt_hd, t1, t_start, dtout, nout, fcomm)  # Run HERMESHD setup
    # if iam == 0: print("\n>>>>\n>>>>", 't', t, 'dtout', dtout, 'nout', nout, "\n>>>>\n")
    sbuf_Qres[:,:,:] = Qio[:,:,0,rh:my,0]        # just grab constant basis func
    # if iam == 0: print('sbuf_Qres\n',sbuf_Qres[:,:,rh],'END\n')
    # if iam == 0: print('Qio\n',Qio[:,:,:,rh,0],'END\n')
    comm.Barrier()
    comm.Gather(sbuf_Qres, rbuf_Qres, root=0)  # Gather Qres from each rank

    #-------------------------------
    ux_fname = 'ux_per_atom.dat'
    fname = ux_to_atomfile(natoms, rbuf_Qres, fname=ux_fname)
    lmp.command("variable ux atomfile {}".format(fname))
    # lmp.command("variable ux atom 0.01")
    lmp.command("variable uy atom 0.0")
    lmp.command("variable uz atom 0.0")

    lmp.command("region rid_top block {} {} {} {} {} {} units box".format(Lxd, Lxu, tb_lo, tb_hi, Lzd, Lzu))
    lmp.command("variable hfx atom \"-v_zeta*(vx - v_ux) + normal(0, v_std, {})\"".format(seed))
    lmp.command("variable hfy atom \"-v_zeta*(vy - v_uy) + normal(0, v_std, {})\"".format(seed))
    lmp.command("fix hforce all addforce v_hfx v_hfy 0 every 1 region rid_top")
    #-------------------------------

    ### Run ####################################################################
    pre, post = 'yes', 'no'
    nloop = 10000

    lmp.command("dump        dcdf all dcd {} {}/md_{}.dcd".format(100, datdir, bname))
    lmp.command("dump_modify dcdf pbc yes")
    lmp.command("dump        lmpf all custom {} {}/md_{}.dump id type x y z vx vy vz v_ux v_uy v_uz".format(100, datdir, bname))
    lmp.command("dump_modify lmpf pbc yes")

    for step in xrange(nloop):
        # LAMMPS
        if step == nloop: post = 'yes'
        run_lammps(lmp, 10, dt=dt_md, nout=100, pre=pre, post=post)
        pre = 'no'
        comm.Barrier()
        # lmp.command('print ">>> The x-velocity of the fluid is: ${ux}"')

        # HERMES
        run_hermeshd(1, Qio, Q1, Q2, t, dt_hd, t1, dtout, nout)
        # if iam == 0: print("\n>>>>\n>>>>", 't', t, 'dtout', dtout, 'nout', nout, "\n>>>>\n")
        comm.Barrier()
        sbuf_Qres[:,:,:] = Qio[:,:,0,rh:my,0]                    # Get reservoir vars from Qio
        comm.Barrier()
        comm.Gather(sbuf_Qres, rbuf_Qres, root=0)  # Gather Qres from each rank
        comm.Barrier()
        fname = ux_to_atomfile(natoms, rbuf_Qres, fname=ux_fname)
        # comm.Barrier()
        lmp.command("variable ux atomfile {}".format(fname))
        lmp.command("variable hfx atom \"-v_zeta*(vx - v_ux) + normal(0, v_std, {})\"".format(seed))
        lmp.command("variable hfy atom \"-v_zeta*(vy - v_uy) + normal(0, v_std, {})\"".format(seed))
        lmp.command("unfix hforce")
        lmp.command("fix hforce all addforce v_hfx v_hfy 0 every 1 region rid_top")


    finalize(lmp)
    cleanup_fhd(t_start)
    ############################################################################
