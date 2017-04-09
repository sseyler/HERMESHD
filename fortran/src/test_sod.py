from mpi4py import MPI
import numpy as np

import hermeshd
# import hermes

###########################################
# This works
#------------------------------------------
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
# print "Proc %d out of %d procs" % (rank,size)

# fcomm = MPI.COMM_WORLD.py2f()
# hermeshd.hermeshd.main(fcomm)
#------------------------------------------


###########################################
# Current testing
#------------------------------------------

# temp = hermeshd.hermeshd.temp
#
# a, b, c = np.array(10), np.array(10), np.array(10)
#
# if rank == 0: print "Before:  a = {}  b = {}  c = {}".format(a,b,c)
# stuff = temp(a=a, b=b, c=c)
# if rank == 0: print stuff
# if rank == 0: print "After:   a = {}  b = {}  c = {}".format(a,b,c)


# Get the communicator for Fortran
fcomm = MPI.COMM_WORLD.py2f()

# Alias the Fortran subroutines
main = hermeshd.hermeshd.main
setup = hermeshd.hermeshd.setup
step = hermeshd.hermeshd.step
generate_output = hermeshd.hermeshd.generate_output
cleanup = hermeshd.hermeshd.cleanup

# Instantiate some global parameters
nx, ny, nz = 50, 1, 1
nQ, nB = 11, 8

# Field arrays
Qio = np.empty((nx, ny, nz, nQ, nB), order='F', dtype=float)
Q1  = np.empty((nx, ny, nz, nQ, nB), order='F', dtype=float)
Q2  = np.empty((nx, ny, nz, nQ, nB), order='F', dtype=float)

# Time variables
t  = np.array(0.0, dtype=float)  # works w/ np.float32 and None
tf = np.array(7.0e-4, dtype=float)
dt = np.array(0.0, dtype=float)

# Timing and output variables
t1      = np.array(0.0, dtype=float)
t_start = np.array(0.0, dtype=float)
dtout   = np.array(0.0, dtype=float)
nout    = np.array(0, dtype=int)  # works w/ np.int32 and None


##############################
# I. SETUP
#-----------------------------
# Qio, t, dt, t1, t_start, dtout, nout, fcomm = setup(Qio, t, dt, t1, t_start, dtout, nout, fcomm)
setup(Qio, t, dt, t1, t_start, dtout, nout, fcomm)

# if rank == 0: print "t = {}   dt = {}   nout = {}".format(t, dt, nout)
# print np.ascontiguousarray(Qio)

##############################
# I. SIMULATION
#-----------------------------
while ( t < tf ):
    step(Qio, Q1, Q2, t, dt)
    generate_output(Qio, t, dt, t1, dtout, nout)
    # if rank == 0: print "t = {}   dt = {}   nout = {}".format(t, dt, nout)

##############################
# I. CLEANUP
#-----------------------------
cleanup(t_start)
