import os
import random
import argparse

"""Calls Gnl to genrate a synthetic benchmark and translates it to blif.

Parameters
----------
size : int
    Number of LUTs in the circuit.
p : float
    Rent's exponent.
fixed_o : int
    Fixed number of outputs of the circuit.
filename : str
    Name of the file in which to store the resulting blif.
batch : Optional[bool], default = False
    Generates all benchmarks prescribed by default in this file.

Returns
-------
None
"""

parser = argparse.ArgumentParser()
parser.add_argument("--size")
parser.add_argument("--p")
parser.add_argument("--fixed_o")
parser.add_argument("--filename")
parser.add_argument("--batch")

args = parser.parse_args()

size = int(args.size)

BATCH = False
try:
    BATCH = int(args.batch)
except:
    pass

p = 1.0
if not BATCH:
    p = float(args.p)

filename = args.filename

perimeter_p = 0.36 
#From Christie and Stroobandt, "The Interpretation and Application of Rent's Rule", 2000,
#last paragraph of Section II---region II Rent's exponent of an X86 CPU.

lut_size_distribution = {6 : 0.35, 5 : 0.25, 4 : 0.17, 3 : 0.07, 2 : 0.16}
#From Fig. 9 (balanced) of Hutton et al., "Improving FPGA Performance and Area Using an Adaptive Logic Module", 2004.

ff_to_lut_ratio = 0.7
#From ISPD'16 design [3].
#Also conforms to Figure 4 of Lewis et al., "Architectural Enhancements in Stratix V", 2013.

region2_size = int(0.25 * size)
I = 50
#Arbitrarily chosen.
O = -1
try:
    O = int(args.fixed_o)
except:
    pass

##########################################################################
def export_gnl_input():
    """Exports the Gnl input file.

    Parameters
    ----------
    None
    
    Returns
    -------
    None
    """
    
    circ_name = filename
    try:
        circ_name = os.path.basename(circ_name)
    except:
        pass
    circ_name = circ_name.rsplit('.', 1)[0]

    ff_count = int(ff_to_lut_ratio * size)
    lut_counts = [int(lut_size_distribution[K] * size) for K\
                  in sorted(lut_size_distribution, reverse = True)]
    real_size = sum(lut_counts) + ff_count

    txt = "[library]\n"
    txt += "name=lib\n"
    txt += "latch=dff 1 1\n"
    txt += "gate=lut6 6 1\n"
    txt += "gate=lut5 5 1\n"
    txt += "gate=lut4 4 1\n"
    txt += "gate=lut3 3 1\n"
    txt += "gate=lut2 2 1\n"
    txt += "[circuit]\n"
    txt += "name=%s\n" % circ_name
    txt += "libraries=lib\n"
    txt += "distribution=%s\n" % (' '.join([str(c) for c in [ff_count] + lut_counts]))
    txt += "size=%d\n" % (real_size - region2_size)
    txt += "  p=%.2f\n" % p
    txt += "size=%d\n" % real_size
    if O < 0:
        txt += "  p=%.2f\n" % perimeter_p
    txt += "  I=%d" % I
    if O >= 0:
        txt += "\n  O=%d" % O

    with open(filename.rsplit('.', 1)[0] + ".gnl", "w") as outf:
        outf.write(txt)
##########################################################################

##########################################################################
def call_gnl(seed):
    """Calls Gnl to generate the netlist.

    Parameters
    ----------
    seed : int
        Random number generator seed.

    Returns
    -------
    None
    """

    circ_name = filename
    try:
        circ_name = os.path.basename(circ_name)
    except:
        pass
    circ_name = circ_name.rsplit('.', 1)[0]

    container = "keen_jackson"
    cp_to_container = "docker cp %s.gnl %s:/" % (circ_name, container)
    call = "docker exec -it %s /gnl -seed %d /%s.gnl" % (container, seed, circ_name)
    cp_from_container = "docker cp %s:/%s.hnl ./" % (container, circ_name)

    os.system(cp_to_container)
    os.system(call)
    os.system(cp_from_container)
##########################################################################

##########################################################################
def gnl_to_blif():
    """Converts the netlist output by Gnl to blif.

    Parameters
    ----------
    None

    Returns
    -------
    None
    """

    #------------------------------------------------------------------------#
    def replace_ff(line):
        """Replaces an FF instantiation with the appropriate blif latch.

        Parameters
        ----------
        line : str
            A Gnl FF instantiation.

        Returns
        -------
        str
            A blif latch instantiation.
        """

        ff, d, q = line.split()
       
        return ".latch %s %s re clk 2\n" % (d, q)
    #------------------------------------------------------------------------#

    #------------------------------------------------------------------------#
    def replace_lut(line):
        """Replaces a LUT instantiation with the appropriate blif names.

        Parameters
        ----------
        line : str
            A Gnl LUT instantiation.

        Returns
        -------
        str
            A blif LUT instantiation.
        """

        ins = line.split()[1:-1]
        o = line.split()[-1]

        return ".names %s %s\n%s\n" % (' '.join(ins), o, '1' * len(ins) + " 1")
    #------------------------------------------------------------------------#

    circ_name = filename
    try:
        circ_name = os.path.basename(circ_name)
    except:
        pass
    circ_name = circ_name.rsplit('.', 1)[0]

    with open(circ_name + ".hnl", "r") as inf:
        lines = inf.readlines()

    txt = ".model %s\n" % circ_name
    rd = False
    for line in lines:
        if line.startswith("circuit %s" % circ_name):
            rd = True
            continue
        if not rd:
            continue

        if line.startswith("input"):
            txt += line.replace("input", ".inputs").strip() + " clk\n"
        elif line.startswith("output"):
            txt += line.replace("output", ".outputs")
        elif line.startswith("dff"):
            txt += replace_ff(line)
        elif line.startswith("lut"):
            txt += replace_lut(line)
        elif line.startswith("end"):
            txt += ".end"

    with open("%s.blif" % circ_name, "w") as outf:
        outf.write(txt)
##########################################################################

##########################################################################
def batch():
    """Generates a number of benchmarks.

    Parameters
    ----------
    None

    Returns
    -------
    None
    """

    if not BATCH:
        return

    exponents = [0.7]
    name_template = "synth_p%.2f_size%d_seed%d.blif"

    sample_size = 10
    min_seed, max_seed = 0, 100000
    random.seed(0)
    seeds = [random.randint(min_seed, max_seed) for s in range(0, sample_size)]

    global p
    global filename

    for p in exponents:
        for seed in seeds:
            filename = name_template % (p, size, seed)
            export_gnl_input()
            call_gnl(seed)
            gnl_to_blif()

    os.system("mkdir gnl_benchmarks/")
    os.system("mv synth* gnl_benchmarks/")
    exit(0)
##########################################################################

batch()
export_gnl_input()
call_gnl()
gnl_to_blif()
