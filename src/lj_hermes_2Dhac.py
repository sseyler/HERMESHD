#!/usr/bin/python
from __future__ import print_function, division
from mpi4py import MPI
import numpy as np

import hermeshd

comm = MPI.COMM_WORLD
iam = comm.Get_rank()
nproc = comm.Get_size()

################################################################################
# Alias the Fortran subroutines
main_fhd  = hermeshd.hermeshd.main
setup_fhd = hermeshd.hermeshd.setup
step_fhd  = hermeshd.hermeshd.step
output_fhd  = hermeshd.hermeshd.generate_output
cleanup_fhd = hermeshd.hermeshd.cleanup

xc = hermeshd.spatial.xc
yc = hermeshd.spatial.yc
zc = hermeshd.spatial.zc

##############################################
# Instantiate some global parameters
rh = 0                     # density
mx, my, mz = 1, 2, 3       # vector momentum
en = 4                     # scalar energy
exx, eyy, ezz = 5, 6, 7    # isotropic stress
exy, exz, eyz = 8, 9, 10   # deviatoric stress

mpi_nx, mpi_ny, mpi_nz = 4, 4, 1
nx, ny, nz = 2, 2, 1
nQ, nB = 11, 8


##############################################
# Field arrays
Qio = np.empty((nx, ny, nz, nQ, nB), order='F', dtype=np.float32)
Q1  = np.empty((nx, ny, nz, nQ, nB), order='F', dtype=np.float32)
Q2  = np.empty((nx, ny, nz, nQ, nB), order='F', dtype=np.float32)

Qres_shape = (nx, ny, nz, nQ)
Qres = np.empty(Qres_shape, order='F', dtype=np.float32)

sendbuf_Qres = None
if iam == 0:
    recvbuf_Qres = np.empty((nproc,)+Qres_shape, order='F', dtype=np.float32)
else:
    recvbuf_Qres = None
################################################################################

# Time variables
t  = np.array(0.0,   dtype=float)  # works w/ np.float32 and None
tf = np.array(1.0e2, dtype=float)

# Timing and output variables
t1      = np.array(0.0, dtype=float)
t_start = np.array(0.0, dtype=float)
dtout   = tf/1000
nout    = np.array(0,   dtype=int)  # works w/ np.int32 and None


def run_hermeshd(nstep_hd, Qio, Q1, Q2, t, dt, t1, dtout, nout):
    for i in xrange(nstep_hd):
        step_fhd(Qio, Q1, Q2, t, dt)
        output_fhd(Qio, t, dt, t1, dtout, nout)


def get_Qres_2Dhac(x, y):
    global recvbuf_Qres
    gid_i, gid_j, gid_k = get_global_idx(x, y, 0)
    mpi_i, mpi_j, mpi_k = get_mpi_idx_from_gid(gid_i, gid_j, gid_k)
    i, j, k = get_local_idx(gid_i, gid_j, gid_k, mpi_i, mpi_j, mpi_k)
    rank = get_mpi_rank(mpi_i, mpi_j, mpi_k)
    return recvbuf_Qres[rank,i-1,j-1,k-1,1]

def get_Qres_var(x, y, z, fv):
    global recvbuf_Qres
    gid_i, gid_j, gid_k = get_global_idx(x, y, z)
    mpi_i, mpi_j, mpi_k = get_mpi_idx_from_gid(gid_i, gid_j, gid_k)
    i, j, k = get_local_idx(gid_i, gid_j, gid_k, mpi_i, mpi_j, mpi_k)
    rank = get_mpi_rank(mpi_i, mpi_j, mpi_k)
    return recvbuf_Qres[rank,i-1,j-1,k-1,fv]


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
