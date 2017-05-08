#!/usr/bin/python

from mpi4py import MPI

# Alias the Fortran subroutines
main = hermeshd.hermeshd.main
setup = hermeshd.hermeshd.setup
step = hermeshd.hermeshd.step
generate_output = hermeshd.hermeshd.generate_output
cleanup = hermeshd.hermeshd.cleanup

# Instantiate some global parameters
nx, ny, nz = 5, 5, 5
nQ, nB = 11, 8

# Field arrays
Qio = np.empty((nx, ny, nz, nQ, nB), order='F', dtype=np.float32)
Q1  = np.empty((nx, ny, nz, nQ, nB), order='F', dtype=np.float32)
Q2  = np.empty((nx, ny, nz, nQ, nB), order='F', dtype=np.float32)

# Time variables
t  = np.array(0.0,   dtype=float)  # works w/ np.float32 and None
tf = np.array(1.0e2, dtype=float)

# Timing and output variables
t1      = np.array(0.0, dtype=float)
t_start = np.array(0.0, dtype=float)
dtout   = tf/1000
nout    = np.array(0,   dtype=int)  # works w/ np.int32 and None


def run_hermes(nstep_hd, Qio, Q1, Q2, t, dt, t1, dtout, nout):
    for i in xrange(nstep_hd):
        step(Qio, Q1, Q2, t, dt)
        generate_output(Qio, t, dt, t1, dtout, nout)
