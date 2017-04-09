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
Qio = np.empty((nx, ny, nz, nQ, nB), order='F', dtype=np.float32)
Q1  = np.empty((nx, ny, nz, nQ, nB), order='F', dtype=np.float32)
Q2  = np.empty((nx, ny, nz, nQ, nB), order='F', dtype=np.float32)

# Time variables
t  = 0.0
tf = 7.0e-4
dt = 0.0

# Timing and output variables
t1      = 0.0
t_start = 0.0
dtout = 0.0
nout  = 0



temp = hermeshd.hermeshd.temp

a, b, c, d = 0, 10, 20, 30

if rank == 0: print "Before:  a = {}  b = {}  c = {}  d = {}".format(a,b,c,d)
stuff = temp(a, b, c)
if rank == 0: print stuff
if rank == 0: print "After:   a = {}  b = {}  c = {}  d = {}".format(a,b,c,d)

# ##############################
# # I. SETUP
# #-----------------------------
# # Qio, t, dt, t1, t_start, dtout, nout, fcomm = setup(Qio, t, dt, t1, t_start, dtout, nout, fcomm)
# stuff = setup(Qio, t, dt, t1, t_start, dtout, nout, fcomm)
#
# Qio, t, dt, t1, t_start, dtout, nout, fcomm = stuff
# print "t = {}   dt = {}   nout = {}".format(t, dt, nout)
# # print np.ascontiguousarray(Qio)
#
# ##############################
# # I. SIMULATION
# #-----------------------------
# while ( t < tf ):
#     Qio, Q1, Q2, t, dt = step(Qio, Q1, Q2, t, dt)
#     generate_output(Qio, t, dt, t1, nout)
#     # print "t = {}   dt = {}   nout = {}".format(t, dt, nout)
#
# ##############################
# # I. CLEANUP
# #-----------------------------
# cleaup(t_start)
