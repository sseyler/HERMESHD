#!/usr/bin/python
from __future__ import print_function, division
from mpi4py import MPI

import os, sys, getopt

import numpy as np
from lammps import lammps, PyLammps

comm = MPI.COMM_WORLD
iam = comm.Get_rank()
nproc = comm.Get_size()

seed = 12345

bname  = "lj_pylmp_2Dhac"  # base name of simulation files
datdir = "data/{}_test".format(bname)    # name of output data directory

### Vars #######################
eta = 2.084e-3  # in Poise (= 0.1 Pascal-seconds = 0.1 Pa s = g/cm/s)
g = 45.5  # lattice geometry factor (from Giupponi, et al. JCP 2007)
dx = 10.0  # HD cell size in Angstroms
zeta_bare = 1.0  # bare friction coefficient
zeta_eff = 1./zeta_bare + 1/(g*eta*dx)

te_init = 50.0
te_sim  = 100.0
pr_init = 1.0
pr_sim  = 1.0

dt_md  = 10.0      # timestep in fs
n_md = 10
t_md = n_md*dt_md

################################
nmin   = 200      # number of emin steps
nsteps = 10000     # number of timesteps
nout   = 100      # output frequency

mass = 39.948
La  = 6.0          # lattice spacing in A
bdy = 10.0         # buffer dz (height of each reservoir)
a   = bdy/1.0      # grid cell dimensions in A

L  = 80.0              # length of single box dimension in A  # 34.68 or 78.45
Lx = Ly = L
Lz = L/12  # + 2*bdz   # z direction length is 4*L/3 + total height of top/bottom reservoirs
dy = dx            # grid cell dimensions in A
dz = Lz

Lxu =  Lx/2.0
Lxd = -Lx/2.0
Lyu =  Ly/2.0
Lyd = -Ly/2.0
Lzu =  Lz/2.0
Lzd = -Lz/2.0

xxu = Lxu/La
yyu = Lyu/La
zzu = Lzu/La
xxd = Lxd/La
yyd = Lyd/La
zzd = Lzd/La

# x-direction
bb_lo = Lyd
bb_hi = Lyd + bdy
tb_lo = Lyu - bdy
tb_hi = Lyu

LJe = 0.23748
LJs = 3.4
LJc = 12.0

LJWe = 0.75*LJe
LJWs = 0.75*LJs
LJWc = bdy
################################

# frequencies for taking averages for buffer region particles
neve, nrep, nfre = 2, 5, 10   # avg over every 2 steps, 5 times (over 10 total steps)

# Output files of particles in buffer region for
#  DENSITY
rh_s_file = "{}/rh.sim".format(datdir)
grid_file = "{}/grid.sim".format(datdir)

rh_b_file = "{}/rh.bbuffer".format(datdir)
rh_t_file = "{}/rh.tbuffer".format(datdir)
#  VELOCITY
ux_b_file = "{}/ux.bbuffer".format(datdir)
ux_t_file = "{}/ux.tbuffer".format(datdir)
uy_b_file = "{}/uy.bbuffer".format(datdir)
uy_t_file = "{}/uy.tbuffer".format(datdir)
#  TEMPERATURE
te_b_file = "{}/te.bbuffer".format(datdir)
te_t_file = "{}/te.tbuffer".format(datdir)


def setup(lmp):
    print_mpi("Setting up simulation...\n")

    inputfile = read(sys.argv[1:])
    if (inputfile): lmp.file(inputfile)

    lmp.command("units real")
    lmp.command("newton on")

    lmp = create_box(lmp)

    lmp.command("pair_style  lj/cut {}".format(LJc))
    lmp.command("pair_coeff  1 1 {} {} {}".format(LJe, LJs, LJc))
    lmp.command("pair_modify shift yes")
    lmp.command("neighbor     3.0 bin")
    lmp.command("neigh_modify delay 0 every 20 check no")
    lmp.command("thermo_style multi")

    # WARN: LAMMPS claims the system must be init before write_dump can be used...
    lmp.command("write_dump all xyz {}/init_{}.xyz".format(datdir, bname))
    xyz_to_pdb("{}/init_{}.xyz".format(datdir, bname))
    lmp.command("restart {} {}_a.res {}_b.res".format(1000, bname, bname))
    return lmp

def create_box(lmp):
    lmp.command("dimension    2")
    lmp.command("boundary     p f p")
    lmp.command("atom_style   atomic")
    lmp.command("atom_modify  map hash")
    lmp.command("lattice      sq2 {} origin 0.0 0.25 0.0".format(La))

    lmp.command("region simreg block {} {} {} {} {} {} units box".format(Lxd,Lxu,Lyd,Lyu,Lzd,Lzu))
    lmp.command("region latreg block {} {} {} {} {} {} units lattice".format(xxd,xxu,yyd+LJWs/10,yyu-LJWs/10,zzd,zzu))
    lmp.command("create_box   1 simreg")
    lmp.command("create_atoms 1 region latreg units box")
    lmp.command("mass  1 {}".format(mass))
    return lmp


def init_velocities(lmp):
    lmp.command("velocity all create {} 87287 loop geom".format(te_init))
    return lmp


def setup_buffer(lmp):
    # STEP 1: Define a "chunk" of atoms with an implicit buffer region
    lmp.command("compute cid_bot all chunk/atom bin/1d y {} {} discard yes bound y {} {} units box".format(bb_lo, bdy, bb_lo, bb_hi))
    lmp.command("compute cid_top all chunk/atom bin/1d y {} {} discard yes bound y {} {} units box".format(tb_lo, bdy, tb_lo, tb_hi))
    # STEP 2a: Use the pre-defined "chunk" from step 1 to compute an average DENSITY
    lmp.command("fix rh_bot all ave/chunk {} {} {} cid_bot density/mass norm sample ave one file {}".format(neve, nrep, nfre, rh_b_file))
    lmp.command("fix rh_top all ave/chunk {} {} {} cid_top density/mass norm sample ave one file {}".format(neve, nrep, nfre, rh_t_file))
    # STEP 2b: Use the pre-defined "chunk" from step 1 to compute average VELOCITIES
    lmp.command("fix ux_bot all ave/chunk {} {} {} cid_bot vx norm sample ave one file {}".format(neve, nrep, nfre, ux_b_file))
    lmp.command("fix ux_top all ave/chunk {} {} {} cid_top vx norm sample ave one file {}".format(neve, nrep, nfre, ux_t_file))
    lmp.command("fix uy_bot all ave/chunk {} {} {} cid_bot vy norm sample ave one file {}".format(neve, nrep, nfre, uy_b_file))
    lmp.command("fix uy_top all ave/chunk {} {} {} cid_top vy norm sample ave one file {}".format(neve, nrep, nfre, uy_t_file))
    # STEP 2c: Use the pre-defined "chunk" from step 1 to compute an average TEMPERATURE (OR PRESSURE)
    lmp.command("compute te_bot_t all temp/chunk cid_bot temp com yes")
    lmp.command("compute te_top_t all temp/chunk cid_top temp com yes")
    lmp.command("fix te_bot all ave/time {} {} {} c_te_bot_t ave one file {}".format(neve, nrep, nfre, te_b_file))
    lmp.command("fix te_top all ave/time {} {} {} c_te_top_t ave one file {}".format(neve, nrep, nfre, te_t_file))
    return lmp


def setup_md_region(lmp):
    lmp.command("compute cid_sim all chunk/atom bin/1d y {} {} discard yes bound y {} {} units box".format(bb_hi, tb_lo-bb_hi, bb_hi, tb_lo))
    lmp.command("fix rh_sim  all ave/chunk {} {} {} cid_sim density/mass norm sample ave one file {}".format(neve, nrep, nfre, rh_s_file))
    return lmp


def setup_wall(lmp):
    lmp.command("fix wall_ylo all wall/lj126 ylo {} {} {} {} units box".format(Lyd, LJWe, LJWs, LJWc))
    lmp.command("fix wall_yhi all wall/lj126 yhi {} {} {} {} units box".format(Lyu, LJWe, LJWs, LJWc))
    return lmp


def finalize(lmp):
    lmp.close()
    MPI.Finalize()


def minimize(lmp, style='cg'):
    print_mpi(">>> Minimizing for {} steps...".format(nmin))
    lmp.command("thermo     100")
    lmp.command("dump       emin all dcd {} {}/em_{}.dcd".format(10, datdir, bname))

    lmp.command("min_style {}".format(style))
    lmp.command("minimize   0.0 0.0 {} {}".format(nmin, 100*nmin))
    lmp.command("write_dump all xyz {}/em_{}.xyz".format(datdir, bname))
    xyz_to_pdb("{}/em_{}.xyz".format(datdir, bname))
    lmp.command("undump emin")
    return lmp


def equilibrate(lmp, te_i, te_f, fast=False):
    nsteps_eq = 1000 if fast else 25000
    print_mpi(">>> NVT equilibration for {} steps...".format(nsteps_eq))
    lmp.command("thermo   5000")
    lmp.command("timestep 1.0")
    lmp.command("fix  1   all nvt temp {} {} 100.0 tchain 1".format(te_i, te_f))
    lmp.command("fix  e2d all enforce2d")
    lmp.command("dump eq1 all dcd {} {}/eq1_{}.dcd".format(100, datdir, bname))

    lmp.command("run      {}".format(nsteps_eq))
    lmp.command("write_dump all xyz {}/eq1_{}.xyz".format(datdir, bname))
    xyz_to_pdb("{}/eq1_{}.xyz".format(datdir, bname))
    lmp.command("unfix 1")
    lmp.command("undump eq1")
    return lmp


def run_lammps(lmp, nstep, dt=dt_md, nout=nout):
    print_mpi(">>> Running NVE simulation for {} steps...".format(nstep))
    lmp.command("thermo   {}".format(nout))
    lmp.command("timestep {}".format(dt))
    lmp.command("fix  1   all nve")
    lmp.command("fix  e2d all enforce2d")

    lmp.command("run      {}".format(nstep))
    # lmp.command("write_dump all xyz {}/md_{}.xyz".format(datdir, bname))
    # xyz_to_pdb("{}/md_{}.xyz".format(datdir, bname))
    lmp.command("write_restart {}.res".format(bname))
    return lmp


def read(argv):
    inputfile = ''
    outputfile = ''
    try:
        opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
    except getopt.GetoptError:
        print('test_lj.py -i <inputfile> -o <outputfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('test_lj.py -i <inputfile> -o <outputfile>')
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



mpi_nx, mpi_ny, mpi_nz = 4, 4, 1
nx, ny, nz = 2, 2, 1
nQ, nB = 11, 8

def get_gid(x, y, z):
    gid_i = int(np.ceil((x - Lxd) / dx))
    gid_j = int(np.ceil((y - Lyd) / dy))
    gid_k = int(np.ceil((z - Lzd) / dz))
    return gid_i, gid_j, gid_k

# def get_local_cell_idx(x, y, z):
#     lid_i = int(np.ceil((x - loc_Lxd) / dx))
#     lid_j = int(np.ceil((y - loc_Lyd) / dy))
#     lid_k = int(np.ceil((z - loc_Lzd) / dz))
#     return lid_i, lid_j, lid_k

def get_mpi_idx_from_gid(gid_i, gid_j, gid_k):
    mpi_i = int( np.ceil(gid_i / nx) )
    mpi_j = int( np.ceil(gid_j / ny) )
    mpi_k = int( np.ceil(gid_k / nz) )
    return mpi_i, mpi_j, mpi_k

def get_lid_from_gid(gid_i, gid_j, gid_k, mpi_i, mpi_j, mpi_k):
    i = int( gid_i - (mpi_i-1)*nx )
    j = int( gid_j - (mpi_j-1)*ny )
    k = int( gid_k - (mpi_k-1)*nz )
    return i, j, k

def get_mpi_idx_and_lid(gid_i, gid_j, gid_k):
    mpi_i, mpi_j, mpi_k = get_mpi_idx_from_gid(gid_i, gid_j, gid_k)
    i = int( gid_i - (mpi_i-1)*nx )
    j = int( gid_j - (mpi_j-1)*ny )
    k = int( gid_k - (mpi_k-1)*nz )
    return [(mpi_i, i), (mpi_j, j), (mpi_k, k)]

def get_mpi_rank(mpi_i, mpi_j, mpi_k):
    return mpi_ny*mpi_nx*(mpi_k-1) + mpi_nx*(mpi_j-1) + (mpi_i-1)


if __name__ == "__main__":
    lmp = lammps()

    lmp = setup(lmp)
    lmp = setup_wall(lmp)
    lmp = setup_buffer(lmp)
    lmp = setup_md_region(lmp)
    lmp = init_velocities(lmp)
    # lmp = minimize(lmp)
    lmp = equilibrate(lmp, te_init, te_sim, fast=True)

    ################################
    lmp.command("dump     run all dcd {} {}/md_{}.dcd".format(100, datdir, bname))
    lmp.command("dump_modify run pbc yes")

    lmp.command("variable dt equal {}".format(dt_md))
    lmp.command("variable zeta  equal {}".format(0.1))
    lmp.command("variable kt    equal {}".format(1.9872036e-3*te_sim))  # gas constant in kcal/mol
    lmp.command("variable sigma equal sqrt(2*v_kt*v_zeta/v_dt)")
    lmp.command("variable ux equal 0.0")
    lmp.command("variable uy equal 0.01")
    lmp.command("variable hfx atom \"-v_zeta*(vx - v_ux) + normal(0.0, v_sigma, {})\"".format(seed))
    lmp.command("variable hfy atom \"-v_zeta*(vy - v_uy) + normal(0.0, v_sigma, {})\"".format(seed))




    ################################################################################
    from ctypes import *
    natoms = lmp.get_natoms()
    atom_pos = None

    xc_md = lmp.gather_atoms("x", 1, 2)  # WARN: Must be executed by all MPI ranks
    # x_md_ptr = lmp.extract_atom("x", 3)
    # lmp.command("variable xc atom x")
    # lmp.command("variable yc atom y")
    # x_md = lmp.extract_variable("xc", "all", 1)
    # y_md = lmp.extract_variable("yc", "all", 1)
    # if iam == 0:
    #     # atom_pos = np.zeros((natoms,2))
    #     # for aid in xrange(natoms):
    #     #     atom_pos[aid,0] = xc_md[2*aid]
    #     #     atom_pos[aid,1] = xc_md[2*aid+1]
    #     atom_pos = np.fromiter(xc_md, dtype=np.float, count=2*natoms).reshape((natoms,2)) # copy
    #     for aid in xrange(natoms):
    #         xc = atom_pos[aid,0]
    #         yc = atom_pos[aid,1]
    #         print( aid+1, xc, yc, ' | ', get_gid(xc, yc, 0), ' | ', get_mpi_idx_and_lid(*get_gid(xc, yc, 0)))

    Qres_shape = (nx, ny, nz, nQ, 1)
    # Qres = np.empty(Qres_shape, order='F', dtype=np.float32)

    Qres = np.zeros(Qres_shape, dtype=np.float32) + iam

    sendbuf = Qres
    recvbuf = None

    if iam == 0:
        recvbuf = np.empty((nproc,)+Qres_shape, dtype=np.float32)

    # Gather Qres from each rank
    comm.Gather(sendbuf, recvbuf, root=0)

    if iam == 0:
        atom_pos = np.fromiter(xc_md, dtype=np.float, count=2*natoms).reshape((natoms,2)) # copy

        for aid in xrange(natoms):
            xc = atom_pos[aid,0]
            yc = atom_pos[aid,1]
            zc = 0

            gid_i, gid_j, gid_k = get_gid(xc, yc, zc)
            mpi_i, mpi_j, mpi_k = get_mpi_idx_from_gid(gid_i, gid_j, gid_k)

            i, j, k = get_lid_from_gid(gid_i, gid_j, gid_k, mpi_i, mpi_j, mpi_k)
            rank = get_mpi_rank(mpi_i, mpi_j, mpi_k)

            print("Hello, I'm atom",aid,"in rank",rank,"and I have a Qres value of",recvbuf[rank,i-1,j-1,k-1,0,:]) # WARNING: Sean, your indices start at 1, idiot



    ################################################################################






    ### Run #######################
    # lmp.command("region rid_bot block {} {} {} {} {} {} units box".format(Lzd, Lzu, Lyd, Lyu, bb_lo, bb_hi))
    # lmp.command("region rid_top block {} {} {} {} {} {} units box".format(Lzd, Lzu, Lyd, Lyu, tb_lo, tb_hi))
    # lmp.command("group lbuff dynamic all region rid_left  every {}".format(n_md))
    # lmp.command("group rbuff dynamic all region rid_right every {}".format(n_md))
    # lmp.command("fix hf_top all addforce v_hfx v_hfy v_hfz every 1 region rid_top")
    lmp.command("compute grid all chunk/atom bin/2d x lower {} y lower {} ids every units box".format(dx,dy))
    lmp.command("compute ctest all property/chunk grid count coord1 coord2")
    lmp.command("fix ftest all ave/chunk {} {} {} grid density/number norm sample ave one file {}".format(neve, nrep, nfre, grid_file))


    # lmp.command("fix hforce all addforce v_hfx v_hfy 0 every 1")

    for i in xrange(1):
        lmp = run_lammps(lmp, 10000, dt=dt_md, nout=500)
