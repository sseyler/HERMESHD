#!/usr/bin/python
from mpi4py import MPI

import os, sys, getopt
import numpy as np
from lammps import lammps

comm = MPI.COMM_WORLD
iam = comm.Get_rank()


bname  = "lj_pylmp_hac"  # base name of simulation files
datdir = "data/{}_test".format(bname)    # name of output data directory

te_init = 30.0
te_sim  = 94.4
pr_init = 1.0
pr_sim  = 1.0

nmin   = 1000      # number of emin steps
tstep  = 10.0      # timestep in fs
nsteps = 10000     # number of timesteps
nout   = 100      # output frequency

mass = 39.948
La = 6.0       # lattice spacing in A
L  = 34.68     #78.45  # length of single box dimension in A
Lx = L*1.50
Ly = L
Lz = L

x  = 1.0/La
y  = 1.0/La
z  = 1.0/La
xx = Lx*x
yy = Ly*y
zz = Lz*z

bdx = Lx/5.0  # buffer dx (width of each buffer slab)
# x-direction
llo = 0.0
lhi = llo + bdx
rhi = Lx
rlo = rhi - bdx

lje = 0.23748
ljs = 3.4
ljcut = 12.0

# frequencies for taking averages for buffer region particles
neve, nrep, nfre = 2, 5, 10   # avg over every 2 steps, 5 times (over 10 total steps)

# Output files of particles in buffer region for
#  DENSITY
rh_l_file = "{}/rh.lbuffer".format(datdir)
rh_r_file = "{}/rh.rbuffer".format(datdir)
#  VELOCITY
ux_l_file = "{}/ux.lbuffer".format(datdir)
ux_r_file = "{}/ux.rbuffer".format(datdir)
uy_l_file = "{}/uy.lbuffer".format(datdir)
uy_r_file = "{}/uy.rbuffer".format(datdir)
uz_l_file = "{}/uz.lbuffer".format(datdir)
uz_r_file = "{}/uz.rbuffer".format(datdir)
#  TEMPERATURE
te_l_file = "{}/te.lbuffer".format(datdir)
te_r_file = "{}/te.rbuffer".format(datdir)


def setup(lmp):
    print_mpi("Setting up simulation...\n")

    inputfile = read(sys.argv[1:])
    if (inputfile): lmp.file(inputfile)

    lmp.command("units real")
    lmp.command("newton on")

    lmp.command("dimension    3")
    lmp.command("boundary     f p p")
    lmp.command("atom_style   atomic")
    lmp.command("atom_modify  map hash")
    lmp.command("lattice      fcc {}".format(La))

    lmp.command("region mybox block 0 {} 0 {} 0 {}".format(xx, yy, zz))
    lmp.command("create_box   1 mybox")
    lmp.command("create_atoms 1 region mybox units box")

    lmp.command("mass  1 {}".format(mass))

    lmp.command("pair_style  lj/cut {}".format(ljcut))
    lmp.command("pair_coeff  1 1 {} {} {}".format(lje, ljs, ljcut))
    lmp.command("pair_modify shift yes")
    lmp.command("neighbor     3.0 bin")
    lmp.command("neigh_modify delay 0 every 20 check no")

    lmp.command("velocity all create {} 87287 loop geom".format(te_init))
    lmp.command("thermo_style multi")

    # WARN: LAMMPS claims the system must be init before write_dump can be used...
    lmp.command("write_dump all xyz {}/init_{}.xyz".format(datdir, bname))
    xyz_to_pdb("{}/init_{}.xyz".format(datdir, bname))
    lmp.command("restart {} {}_a.res {}_b.res".format(1000, bname, bname))
    return lmp


def setup_buffer(lmp):
    # STEP 1: Define a "chunk" of atoms with an implicit buffer region
    lmp.command("compute cid_left  all chunk/atom bin/1d x {} {} discard yes bound x 0.0 0.2 units reduced".format("lower", 0.2))
    lmp.command("compute cid_right all chunk/atom bin/1d x {} {} discard yes bound x 0.8 1.0 units reduced".format(0.8, 0.2))

    # STEP 2: Use the pre-defined "chunk" from step 1 to compute an average
    # DENSITY
    lmp.command("fix rh_left  all ave/chunk {} {} {} cid_left  density/mass norm sample ave one file {}".format(neve, nrep, nfre, rh_l_file))
    lmp.command("fix rh_right all ave/chunk {} {} {} cid_right density/mass norm sample ave one file {}".format(neve, nrep, nfre, rh_r_file))
    # VELOCITIES
    lmp.command("fix ux_left  all ave/chunk {} {} {} cid_left  vx norm sample ave one file {}".format(neve, nrep, nfre, ux_l_file))
    lmp.command("fix ux_right all ave/chunk {} {} {} cid_right vx norm sample ave one file {}".format(neve, nrep, nfre, ux_r_file))
    lmp.command("fix uy_left  all ave/chunk {} {} {} cid_left  vy norm sample ave one file {}".format(neve, nrep, nfre, uy_l_file))
    lmp.command("fix uy_right all ave/chunk {} {} {} cid_right vy norm sample ave one file {}".format(neve, nrep, nfre, uy_r_file))
    lmp.command("fix uz_left  all ave/chunk {} {} {} cid_left  vz norm sample ave one file {}".format(neve, nrep, nfre, uz_l_file))
    lmp.command("fix uz_right all ave/chunk {} {} {} cid_right vz norm sample ave one file {}".format(neve, nrep, nfre, uz_r_file))
    # TEMPERATURE (OR PRESSURE)
    lmp.command("compute te_left_t  all temp/chunk cid_left  temp com yes")
    lmp.command("compute te_right_t all temp/chunk cid_right temp com yes")
    lmp.command("fix te_left  all ave/time {} {} {} c_te_left_t  ave one file {}".format(neve, nrep, nfre, te_l_file))
    lmp.command("fix te_right all ave/time {} {} {} c_te_right_t ave one file {}".format(neve, nrep, nfre, te_r_file))

    return lmp


def setup_wall(lmp):
    # Create region that's the size of the whole simulation domain
    # lmp.command("region rid_wall block {} {} {} {} {} {} units box".format(0.0, Lx, 0.0, Ly, 0.0, Lz))
    # Setup boundary using the wall/region fix
    # lmp.command("fix wall all wall/region mybox lj126 {} {} 8.0".format(lje, ljs))
    lmp.command("fix wall_xlo all wall/lj1043 xlo {} {} {} {} units box".format(-0.25*ljs, 0.1*lje, 0.5*ljs, 0.5*ljs))
    lmp.command("fix wall_xhi all wall/lj1043 xhi {} {} {} {} units box".format(Lx+0.25*ljs, 0.1*lje, 0.5*ljs, 0.5*ljs))
    return lmp


def finalize(lmp):
    lmp.close()
    MPI.Finalize()


def minimize(lmp):
    print_mpi("Minimizing for {} steps...".format(nmin))
    lmp.command("thermo     10")
    lmp.command("dump       emin all dcd {} {}/em_{}.dcd".format(10, datdir, bname))

    lmp.command("minimize   0.0 0.0 {} {}".format(nmin, 10*nmin))
    lmp.command("write_dump all xyz {}/em_{}.xyz".format(datdir, bname))
    xyz_to_pdb("{}/em_{}.xyz".format(datdir, bname))

    lmp.command("undump emin")
    return lmp


def equilibrate(lmp):
    print_mpi("NVT equilibration for 10000 steps...")
    lmp.command("thermo   100")
    lmp.command("timestep 1.0")
    lmp.command("fix      1 all nvt temp {} {} 100.0 tchain 1".format(te_init, te_sim))
    lmp.command("dump     eq1 all dcd {} {}/eq1_{}.dcd".format(100, datdir, bname))

    lmp.command("run      10000")
    lmp.command("write_dump all xyz {}/eq1_{}.xyz".format(datdir, bname))
    xyz_to_pdb("{}/eq1_{}.xyz".format(datdir, bname))

    lmp.command("unfix 1")
    lmp.command("undump eq1")

    print_mpi("NVT equilibration for 10000 steps...")
    lmp.command("thermo   100")
    lmp.command("timestep 5.0")
    lmp.command("fix      1 all nvt temp {} {} 100.0 tchain 1".format(te_sim, te_sim))
    lmp.command("dump     eq2 all dcd {} {}/eq2_{}.dcd".format(100, datdir, bname))

    lmp.command("run      10000")
    lmp.command("write_dump all xyz {}/eq2_{}.xyz".format(datdir, bname))
    xyz_to_pdb("{}/eq2_{}.xyz".format(datdir, bname))

    lmp.command("unfix 1")
    lmp.command("undump eq2")
    return lmp


def run_lammps(lmp, n_md):
    print_mpi("Running NVE simulation for {} steps...".format(n_md))
    lmp.command("thermo {}".format(n_md/10))
    lmp.command("timestep {}".format(tstep))
    lmp.command("fix      1 all nve")
    lmp.command("dump     run all dcd {} {}/md_{}.dcd".format(nout, datdir, bname))

    lmp.command("run      {}".format(n_md))
    lmp.command("write_dump all xyz {}/md_{}.xyz".format(datdir, bname))
    xyz_to_pdb("{}/md_{}.xyz".format(datdir, bname))
    lmp.command("write_restart {}.res".format(bname))
    return lmp


def read(argv):
    inputfile = ''
    outputfile = ''
    try:
        opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
    except getopt.GetoptError:
        print 'test_lj.py -i <inputfile> -o <outputfile>'
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print 'test_lj.py -i <inputfile> -o <outputfile>'
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg
    return inputfile


def xyz_to_pdb(xyzfile):
    import MDAnalysis as mda
    u = mda.Universe(xyzfile)
    u.dimensions = np.array([Lx, Ly, Lz, 90.00, 90.00, 90.00])
    u.atoms.write(os.path.splitext(xyzfile)[0] + '.pdb')


def print_mpi(msg, iam=iam, print_id=0):
    if (iam == print_id): print(msg)



def add_position(x_md, atom_ids, dx, dy, dz):
    for aid in atom_ids:
        x_md[aid][0] += dx
        x_md[aid][1] += dy
        x_md[aid][2] += dz

def add_velocity(v_md, atom_ids, dvx, dvy, dvz):
    for aid in atom_ids:
        v_md[aid][0] += dvx
        v_md[aid][1] += dvy
        v_md[aid][2] += dvz

def add_force(f_md, atom_ids, dfx, dfy, dfz):
    for aid in atom_ids:
        f_md[aid][0] += dfx
        f_md[aid][1] += dfy
        f_md[aid][2] += dfz


if __name__ == "__main__":
    # from ctypes import *

    lmp = lammps()
    lmp = setup(lmp)
    lmp = setup_wall(lmp)

    natoms = lmp.get_natoms()
    # natoms = lmp.extract_global("natoms",0)

    x_md = lmp.extract_atom("x", 3)  # Pointer to underlying array in LAMMPS
    v_md = lmp.extract_atom("v", 3)
    f_md = lmp.extract_atom("f", 3)

    print_mpi('>>> x_md[0][0]: {}'.format(x_md[0][0]))
    print_mpi('>>> x_md[0][1]: {}'.format(x_md[0][1]))
    print_mpi('>>> x_md[0][2]: {}'.format(x_md[0][2]))


    ### Run #######################
    lmp = run_lammps(lmp, 10)

    print_mpi('>>> x_md[0][0]: {}'.format(x_md[0][0]))
    print_mpi('>>> x_md[0][1]: {}'.format(x_md[0][1]))
    print_mpi('>>> x_md[0][2]: {}'.format(x_md[0][2]))

    epsilon = 0.5
    add_position(x_md, [0], epsilon, epsilon, epsilon)

    print_mpi('>>> x_md[0][0]: {}'.format(x_md[0][0]))
    print_mpi('>>> x_md[0][1]: {}'.format(x_md[0][1]))
    print_mpi('>>> x_md[0][2]: {}'.format(x_md[0][2]))







    # f_md = lmp.gather_atoms("f", 1, 3)
    # print('>>> ', len(f_md))
    # print("Global coords from gather_atoms =",f_md[0],f_md[1],f_md[31])

    # natoms = lmp.get_natoms()
    # n3 = 3*natoms
    # x = (n3*c_double)()
    # x[0] = x coord of atom with ID 1
    # x[1] = y coord of atom with ID 1
    # x[2] = z coord of atom with ID 1
    # x[3] = x coord of atom with ID 2
    # print('>>> ', x[1296])
    # x_md = np.asarray(x)
    # print('>>> ', x_md[1296])
    #
    #
    # x_md = lmp.gather_atoms("x", 1, 3)
    # print('>>> Global number of atomic coordinates: ', len(x_md))
    #
    # x_md = lmp.gather_atoms("x", 1, 3)
    # print('>>> Global number of atomic coordinates: ', type(x_md))
    #
    # x_md = lmp.gather_atoms("x", 1, 3)
    # print('>>> x_md[0]: ', x_md[0])
    # x_md = lmp.gather_atoms("x", 1, 3)
    # print('>>> x_md[0,1]: ', x_md[0][1])
    # x_md = lmp.gather_atoms("x", 1, 3)
    # print('>>> x_md[1296]: ', x_md[1296])
    # x_md = lmp.gather_atoms("x", 1, 3)
    # print('>>> x_md[3888]: ', x_md[3888])
