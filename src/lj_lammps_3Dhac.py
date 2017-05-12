#!/usr/bin/python

from mpi4py import MPI

import os, sys, getopt

import numpy as np
from lammps import lammps

comm = MPI.COMM_WORLD
iam = comm.Get_rank()

seed = 12345

bname  = "lj_pylmp_hac"  # base name of simulation files
datdir = "data/{}_test".format(bname)    # name of output data directory

### Vars #######################
eta = 2.084e-3  # in Poise (= 0.1 Pascal-seconds = 0.1 Pa s = g/cm/s)
g = 45.5  # lattice geometry factor (from Giupponi, et al. JCP 2007)
dx = 10.0  # HD cell size in Angstroms
zeta_bare = 1.0  # bare friction coefficient
zeta_eff = 1./zeta_bare + 1/(g*eta*dx)

te_init = 50.0
te_sim  = 94.4
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
bdz = 10.0         # buffer dz (height of each reservoir)
a   = bdz/1.0      # grid cell dimensions in A
dx  = dy = dz = a  # grid cell dimensions in A

L  = 30.0                # length of single box dimension in A  # 34.68 or 78.45
Lz = (4./3.)*L + 2*bdz   # z direction length is 4*L/3 + total height of top/bottom reservoirs
Lx = Ly = L

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
bb_lo = Lzd
bb_hi = Lzd + bdz
tb_lo = Lzu - bdz
tb_hi = Lzu

LJe = 0.23748
LJs = 3.4
LJc = 12.0

LJWe = 0.75*LJe
LJWs = 0.75*LJs
LJWc = bdz
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
uz_b_file = "{}/uz.bbuffer".format(datdir)
uz_t_file = "{}/uz.tbuffer".format(datdir)
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
    lmp.command("dimension    3")
    lmp.command("boundary     p p f")
    lmp.command("atom_style   atomic")
    lmp.command("atom_modify  map hash")
    lmp.command("lattice      fcc {}".format(La))

    lmp.command("region simreg block {} {} {} {} {} {} units box".format(Lxd,Lxu,Lyd,Lyu,Lzd,Lzu))
    lmp.command("region latreg block {} {} {} {} {} {} units lattice".format(xxd,xxu,yyd,yyu,zzd+LJWs/La*1.5,zzu-LJWs/La*1.5))
    lmp.command("create_box   1 simreg")
    lmp.command("create_atoms 1 region latreg units box")
    lmp.command("mass  1 {}".format(mass))
    return lmp


def init_velocities(lmp):
    lmp.command("velocity all create {} 87287 loop geom".format(te_init))
    return lmp


def setup_buffer(lmp):
    # STEP 1: Define a "chunk" of atoms with an implicit buffer region
    lmp.command("compute cid_bot all chunk/atom bin/1d z {} {} discard yes bound z {} {} units box".format(bb_lo, bdz, bb_lo, bb_hi))
    lmp.command("compute cid_top all chunk/atom bin/1d z {} {} discard yes bound z {} {} units box".format(tb_lo, bdz, tb_lo, tb_hi))
    # STEP 2a: Use the pre-defined "chunk" from step 1 to compute an average DENSITY
    lmp.command("fix rh_bot all ave/chunk {} {} {} cid_bot density/mass norm sample ave one file {}".format(neve, nrep, nfre, rh_b_file))
    lmp.command("fix rh_top all ave/chunk {} {} {} cid_top density/mass norm sample ave one file {}".format(neve, nrep, nfre, rh_t_file))
    # STEP 2b: Use the pre-defined "chunk" from step 1 to compute average VELOCITIES
    lmp.command("fix ux_bot all ave/chunk {} {} {} cid_bot vx norm sample ave one file {}".format(neve, nrep, nfre, ux_b_file))
    lmp.command("fix ux_top all ave/chunk {} {} {} cid_top vx norm sample ave one file {}".format(neve, nrep, nfre, ux_t_file))
    lmp.command("fix uy_bot all ave/chunk {} {} {} cid_bot vy norm sample ave one file {}".format(neve, nrep, nfre, uy_b_file))
    lmp.command("fix uy_top all ave/chunk {} {} {} cid_top vy norm sample ave one file {}".format(neve, nrep, nfre, uy_t_file))
    lmp.command("fix uz_bot all ave/chunk {} {} {} cid_bot vz norm sample ave one file {}".format(neve, nrep, nfre, uz_b_file))
    lmp.command("fix uz_top all ave/chunk {} {} {} cid_top vz norm sample ave one file {}".format(neve, nrep, nfre, uz_t_file))
    # STEP 2c: Use the pre-defined "chunk" from step 1 to compute an average TEMPERATURE (OR PRESSURE)
    lmp.command("compute te_bot_t all temp/chunk cid_bot temp com yes")
    lmp.command("compute te_top_t all temp/chunk cid_top temp com yes")
    lmp.command("fix te_bot all ave/time {} {} {} c_te_bot_t ave one file {}".format(neve, nrep, nfre, te_b_file))
    lmp.command("fix te_top all ave/time {} {} {} c_te_top_t ave one file {}".format(neve, nrep, nfre, te_t_file))
    return lmp


def setup_md_region(lmp):
    lmp.command("compute cid_sim all chunk/atom bin/1d z {} {} discard yes bound z {} {} units box".format(bb_hi, 4*bdz, bb_hi, tb_lo))
    lmp.command("fix rh_sim  all ave/chunk {} {} {} cid_sim density/mass norm sample ave one file {}".format(neve, nrep, nfre, rh_s_file))
    return lmp


def setup_wall(lmp):
    lmp.command("fix wall_zlo all wall/lj126 zlo {} {} {} {} units box".format(Lzd, LJWe, LJWs, LJWc))
    lmp.command("fix wall_zhi all wall/lj126 zhi {} {} {} {} units box".format(Lzu, LJWe, LJWs, LJWc))
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


def equilibrate(lmp, te_i, te_f):
    print_mpi(">>> NVT equilibration for 10000 steps...")
    lmp.command("thermo   100")
    lmp.command("timestep 1.0")
    lmp.command("fix      1 all nvt temp {} {} 100.0 tchain 1".format(te_i, te_f))
    lmp.command("dump     eq1 all dcd {} {}/eq1_{}.dcd".format(100, datdir, bname))

    lmp.command("run      10000")
    lmp.command("write_dump all xyz {}/eq1_{}.xyz".format(datdir, bname))
    xyz_to_pdb("{}/eq1_{}.xyz".format(datdir, bname))
    lmp.command("unfix 1")
    lmp.command("undump eq1")
    return lmp


def run_lammps(lmp, nstep, dt=dt_md, nout=nout):
    print_mpi(">>> Running NVE simulation for {} steps...".format(nstep))
    lmp.command("thermo   {}".format(nout))
    lmp.command("timestep {}".format(dt))
    lmp.command("fix      1 all nve")

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


if __name__ == "__main__":
    lmp = lammps()
    lmp = setup(lmp)
    lmp = setup_wall(lmp)
    lmp = setup_buffer(lmp)
    lmp = setup_md_region(lmp)
    lmp = init_velocities(lmp)
    # lmp = minimize(lmp)
    lmp = equilibrate(lmp, te_init, te_sim)

    ################################
    lmp.command("dump     run all dcd {} {}/md_{}.dcd".format(n_md, datdir, bname))
    lmp.command("dump_modify run pbc yes")

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

    ### Run #######################
    # lmp.command("region rid_bot block {} {} {} {} {} {} units box".format(Lzd, Lzu, Lyd, Lyu, bb_lo, bb_hi))
    # lmp.command("region rid_top block {} {} {} {} {} {} units box".format(Lzd, Lzu, Lyd, Lyu, tb_lo, tb_hi))
    # lmp.command("group lbuff dynamic all region rid_left  every {}".format(n_md))
    # lmp.command("group rbuff dynamic all region rid_right every {}".format(n_md))
    # lmp.command("fix hf_top all addforce v_hfx v_hfy v_hfz every 1 region rid_top")
    lmp.command("compute grid all chunk/atom bin/3d x lower {} y lower {} z lower {} ids every units box".format(dx,dy,dz))
    lmp.command("compute ctest all property/chunk grid count coord1 coord2 coord3")
    lmp.command("fix ftest all ave/chunk {} {} {} grid density/number norm sample ave one file {}".format(neve, nrep, nfre, grid_file))

    lmp.command("fix hforce all addforce v_hfx v_hfy v_hfz every 1")

    for i in xrange(1):
        lmp = run_lammps(lmp, 10000, dt=dt_md, nout=100)
