"""Measures the delay of global (intercluster) wires.

Parameters
----------
K : int
    LUT size.
N : int
    Cluster size.
tech : float
    Technology node (16, 7, 5, 4, 3.0, 3.1).
    3.0 corresponds to F3a in the paper and 3.1 to F3b.
insert_rep : Optional[bool], default = False
    Specifies that buffer insertion should be performed.
fixed_length : Optional[int], default = 0
    Specifies the logical length of the wire the delay
    of which should be measured. Otherwise, the maximum
    length of a wire such that its delay is 1% smaller
    than the delay of two wires of the previous length
    joined together through a switch-block multiplexer.
    We enumerate only the lengths that are powers of 2.
    Note that here it is assumed that the wires are aligned
    with one another and with the drivers. In reality where
    the wire may be offset from the driver and/or the sink mux,
    the actual delay will be larger. Hence 1% seems meaningful as
    a threshold.
horizontal : Optional[bool], default = False
    Specifies that the measurement is performed on a horizontal wire.
    otherwise, vertical is assumed.

Returns
-------
Tuple[int]
    Buffer sizes
float
    Measured delay

Notes
-----
Buffer sizes determined here are later used in the more precise model,
where less optimistic alignment between the wire and its driver and the
sink multiplexers is assumed. Thus obtained delay is used for the final
architecture annotation.
"""

import os
import math
import argparse
import sys
sys.path.insert(0,'..')

import setenv
import tech

parser = argparse.ArgumentParser()
parser.add_argument("--K")
parser.add_argument("--N")
parser.add_argument("--tech")
parser.add_argument("--insert_rep")
parser.add_argument("--fixed_length")
parser.add_argument("--horizontal")
args = parser.parse_args()

K = int(args.K)
N = int(args.N)

try:
    tech_node = int(args.tech)
except:
    tech_node = float(args.tech)

INSERT_REP = False
try:
    INSERT_REP = int(args.insert_rep)
except:
    pass

FIXED_LENGTH = 0
try:
    FIXED_LENGTH = int(args.fixed_length)
except:
    pass

HORIZONTAL = False
try:
    HORIZONTAL = int(args.horizontal)
except:
    pass

node_index = tech.nodes.index(tech_node)
device_node_index = tech.nodes.index(int(tech_node))

MxR = tech.MxR[node_index]
MxC = tech.MxC[node_index] * 1e-15

MyR = tech.MyR[node_index]
MyC = tech.MyC[node_index] * 1e-15
GP = tech.GP[device_node_index]
FP = tech.FP[device_node_index]
GL = tech.GL[device_node_index]
vdd = tech.vdd[device_node_index]
Rvia = tech.stacked_via[node_index]

base_tile_height = 2 ** (K - 4) * 2 * (4 + 2) * N * GP
if HORIZONTAL:
    base_tile_height = 2 * (2 ** 4) * 10 * FP

Cw = base_tile_height / 1000.0 * MyC
orig_Cw = Cw

Rw = base_tile_height / 1000.0 * MyR
orig_Rw = Rw

ptm_path = "\"/home/snikolic/FPGA21/ptm/%dnm.l\" %dNM_FINFET_HP\n"

##########################################################################
def update_rc(rep_no, D0, D1, WL):
    """Conservatively scales the tile height due to repeater insertion.

    Parameters
    ----------
    rep_no : int
        Repeater count.
    D0 : int
        Drive strength of the first inverter.
    D1 : int
        Drive strength of subsequent inverters.
    WL : length of the wire.

    Returns
    -------
    None
    """


    inv_inc = lambda d : math.ceil((2 * d + 1) / 25.0)
    #We assume that 25 fin pitches are provided by the mux width.
    #(two SRAMS + row_cnt = 4, for 16:1) +1 is due to the well separation.

    overlapped_wires = WL
    rep_freq = rep_no / float(WL)
    reps_per_tile = overlapped_wires * rep_freq

    total_inv_inc = inv_inc(D0) + inv_inc(D1) if 2 * (D0 + D1 + 1) > 25 else 1
    #We may be able to pack the whole buffer in one gate pitch
    total_inv_inc *= N * GP * (1 + reps_per_tile)
    #There is one wire per LUT.

    scale_factor = (base_tile_height + total_inv_inc) / float(base_tile_height)

    global Cw
    global Rw
    Cw = orig_Cw * scale_factor
    Rw = orig_Rw * scale_factor
##########################################################################

##########################################################################
def gen_header(D0, D1):
    """Outputs the parameters, subcircuits, and voltage sources.

    Parameters
    ----------
    D0 : int
        Drive strength of the first inverter.
    D1 : int
        Drive strength of subsequent inverters.

    Returns
    -------
    str
        SPICE file header.
    """

    txt = ".TITLE MAX_WL\n\n"
    txt += ".LIB %s\n" % (ptm_path % (int(tech_node), int(tech_node)))
    txt += ".TRAN 1p 16n\n.OPTIONS BRIEF=1\n\n"

    txt += ".PARAM Cw=%g\n" % Cw
    txt += ".PARAM Rw=%g\n" % Rw
    txt += ".PARAM Rvia=%g\n\n" % Rvia
    txt += ".PARAM D0=%d\n" % D0
    txt += ".PARAM D1=%d\n\n" % D1
    txt += ".PARAM GL=%dn\n" % GL
    txt += ".PARAM supply_v=%.2f\n\n" % vdd

    #------------------------------------------------------------------------#
    def add_wire_subckt(stages = 1):
        """Adds a single tile length wire subcircuit.

        Parameters
        ----------
        stages : Optional[int], default = 1
            Number of Pi stages in the wire model.
        
        Returns
        -------
        str
            Subcircuit description.
        """

        txt = ".SUBCKT wire n_in n_out\n"

        stage_template = "C%d_1 %s gnd" +  " C='Cw/(2*%d)'\n" % stages
        stage_template += "R%d_1 %s %s" + " R='Rw/%d'\n" % stages
        stage_template += "C%d_2 %s gnd" + " C='Cw/(2*%d)'\n" % stages


        if stages == 1:
            txt += stage_template % (1, "n_in", 1, "n_in", "n_out", 1, "n_out")
        else:
            n_in = "n_in"
            for stage in range(0, stages - 1):
                n_out = "n_out%d" % stage
                txt += stage_template % (stage, n_in, stage, n_in, n_out, stage, n_out)
                n_in = n_out
            stage += 1
            txt += stage_template % (stage, n_out, stage, n_out, "n_out", stage, "n_out")

        txt += ".ENDS\n\n"

        return txt
    #------------------------------------------------------------------------#

    txt += add_wire_subckt(N * 2 ** (K - 4))

    txt += ".SUBCKT buf n_in n_out vdd\n"
    txt += "Rin n_in n_in_via R=Rvia\n"
    txt += "MN1 n_mid n_in_via gnd gnd nmos L=GL nfin=D0\n"
    txt += "MP1 n_mid n_in_via vdd vdd pmos L=GL nfin=D0\n"
    txt += "MN2 n_out_pre_via n_mid gnd gnd nmos L=GL nfin=D1\n"
    txt += "MP2 n_out_pre_via n_mid vdd vdd pmos L=GL nfin=D1\n"
    txt += "Rout n_out_pre_via n_out R=Rvia\n"
    txt += ".ENDS\n\n"

    txt += ".SUBCKT rep n_in n_out vdd\n"
    txt += "Rin n_in n_in_via R=Rvia\n"
    txt += "MN1 n_out_pre_via n_in_via gnd gnd nmos L=GL nfin=D1\n"
    txt += "MP1 n_out_pre_via n_in_via vdd vdd pmos L=GL nfin=D1\n"
    txt += "Rout n_out_pre_via n_out R=Rvia\n"
    txt += ".ENDS\n\n"

    txt += ".SUBCKT tg_on n_in n_out vdd\n"
    txt += "MN1 n_in vdd n_out gnd nmos L=GL nfin=1\n"
    txt += "MP1 n_in gnd n_out vdd pmos L=GL nfin=1\n"
    txt += ".ENDS\n\n"

    txt += ".SUBCKT tg_off n_in n_out vdd\n"
    txt += "MN1 n_in gnd n_out gnd nmos L=GL nfin=1\n"
    txt += "MP1 n_in vdd n_out vdd pmos L=GL nfin=1\n"
    txt += ".ENDS\n\n"


    #------------------------------------------------------------------------#
    def add_mux_subckt(row_cnt, col_cnt, state):
        """Constructs a multiplexer out of minimum-sized transmission gates.

        Parameters
        ----------
        row_cnt : int:
            Number of rows.
        col_cnt : int
            Number of columns.
        state : str
            on, partial, or off.
        
        Returns
        -------
        str
            Description of the multiplexer.
        int
            Advanced index.
        """
   
        ind = 0 

        txt = ".SUBCKT mux_%d_%d_%s n_in n_out vdd\n" % (row_cnt, col_cnt, state) 
        for r in range(0, row_cnt):
            row_node = "n_r_%d" % r
            for c in range(0, col_cnt):
                i = "n_in" if c == 0  and r == 0 else "n_dummy_%d" % (r * col_cnt + c)
                tg_model = "tg_on" if c == 0 and state != "off" else "tg_off"
                txt += "Xtg_%d %s %s vdd %s\n" % (ind, i, row_node, tg_model)
                ind += 1
            tg_model = "tg_on" if r == 0 and state == "on" else "tg_off"
            txt += "Xtg_%d %s %s vdd %s\n" % (ind, row_node, "n_out", tg_model)
            ind += 1

        return txt + ".ENDS\n\n"
    #------------------------------------------------------------------------#

    txt += add_mux_subckt(4, 4, "on")

    txt += "Vps vdd gnd supply_v\n"
    txt += "Vin n_in gnd PULSE (0 supply_v 0 0 0 2n 4n)\n\n"

    txt += "Xmux1 n_in n_mux_out vdd mux_4_4_on\n"
    txt += "Xmux_buf n_mux_out n_buf_out vdd buf\n\n"

    return txt
##########################################################################

##########################################################################
def gen_wire(hops):
    """Generates the entire wire model.

    Parameters
    ----------
    hops : str
        Specifies whether a buffer exist at a given tile (1) or not (0).

    Returns
    -------
    str
        Description of the entire wire.
    """

    n_in = "n_buf_out"
    txt = ""
    for h_cnt, h in enumerate(hops):
        n_out = "n_out%d" % h_cnt
        if int(h):
            txt += "Xrep%d %s %s_rep vdd rep\n" % (h_cnt, n_in, n_out)
            txt += "Xwire%d %s_rep %s wire\n" % (h_cnt, n_out, n_out)
        else:
            txt += "Xwire%d %s %s wire\n" % (h_cnt, n_in, n_out)
        n_in = n_out

    #------------------------------------------------------------------------#
    def insert_load(mux_cnt):
        """Inserts the load into the netlist.

        Parameters
        ----------
        mux_cnt : int
            Number of multiplexers that will form the load.
        
        Returns
        -------
        str
            Load netlist.
        """

        txt = "Rload_via %s n_load_via R=Rvia\n" % n_out
        for i in range(0, mux_cnt):
            txt += "Xmux_load%d n_load_via n_out_load%d vdd mux_4_4_on\n\n" % (i, i)
        
        return txt
    #------------------------------------------------------------------------#

    txt += insert_load(6)

    fall_trig = "RISE=2" if len([h for h in hops if h == '1']) % 2 else "FALL=2"
    rise_trig = "FALL=2" if len([h for h in hops if h == '1']) % 2 else "RISE=2"

    trig_node = "n_in"

    txt += ".MEASURE tfall TRIG V(%s) VAL='supply_v/2' %s\n" % (trig_node, fall_trig)
    txt += "+                  TARG V(n_load_via) VAL supply_v/2 FALL=2\n\n"
    txt += ".MEASURE trise TRIG V(%s) VAL='supply_v/2' %s\n" % (trig_node, rise_trig)
    txt += "+                  TARG V(n_load_via) VAL supply_v/2 RISE=2\n\n"

    txt += ".MEASURE tfall_mux TRIG V(n_in) VAL='supply_v/2' FALL=2\n"
    txt += "+                  TARG V(n_mux_out) VAL supply_v/2 FALL=2\n\n"
    txt += ".MEASURE trise_mux TRIG V(n_in) VAL='supply_v/2' RISE=2\n"
    txt += "+                  TARG V(n_mux_out) VAL supply_v/2 RISE=2\n\n"

    txt += ".MEASURE tfall_buf TRIG V(n_mux_out) VAL='supply_v/2' FALL=2\n"
    txt += "+                  TARG V(n_buf_out) VAL supply_v/2 FALL=2\n\n"
    txt += ".MEASURE trise_buf TRIG V(n_mux_out) VAL='supply_v/2' RISE=2\n"
    txt += "+                  TARG V(n_buf_out) VAL supply_v/2 RISE=2\n\n"

    txt += ".MEASURE tfall_wire TRIG V(n_buf_out) VAL='supply_v/2' FALL=2\n"
    txt += "+                  TARG V(%s) VAL supply_v/2 FALL=2\n\n" % n_out
    txt += ".MEASURE trise_wire TRIG V(n_buf_out) VAL='supply_v/2' RISE=2\n"
    txt += "+                  TARG V(%s) VAL supply_v/2 RISE=2\n\n" % n_out

    txt += ".END"

    return txt
##########################################################################

##########################################################################
def measure(D0, D1, L, rep_no):
    """Constructs the SPICE netlist and measures the delay of a given buffering.

    Parameters
    ----------
    D0 : int
        Drive strength of the first inverter.
    D1 : int
        Drive strength of subsequent inverters.
    L : int
        Length of the wire in tile lengths.
    rep_no : int
        Number of repeaters.
    
    Returns
    -------
    float
        Average between the rise and the fall times.
    """

    update_rc(rep_no, D0, D1, L)

    txt = gen_header(D0, D1)

    if rep_no:
        rep_space = L / (rep_no + 1)
        hops = "0"
        for i in range(1, L):
            hops += '1' if i % rep_space == 0 else '0'
    else:
        hops = '0' * L

    print hops

    txt += gen_wire(hops)

    netlist_filename = "sim_wl_K%dN%dT%s%s.sp" % (K, N, str(tech_node), ("L%d" % L if FIXED_LENGTH else ""))
    hspice_dump = "hspice_K%dN%dT%s%s.dump" % (K, N, str(tech_node), ("L%d" % L if FIXED_LENGTH else ""))
 
    with open(netlist_filename, "w") as outf:
        outf.write(txt)

    hspice_call = os.environ["HSPICE"] + " %s > %s" % (netlist_filename, hspice_dump)
    os.system(hspice_call)

    scale_dict = {'f' : 1e-15, 'p' : 1e-12, 'n' : 1e-9}

    with open(hspice_dump, "r") as inf:
        lines = inf.readlines()

    os.system("rm %s" % hspice_dump)

    for line in lines:
        if "tfall=" in line:
            tfall = float(line.split()[1][:-1]) * scale_dict[line.split()[1][-1]]
        elif "trise=" in line:
            trise = float(line.split()[1][:-1]) * scale_dict[line.split()[1][-1]]

    return (trise + tfall) / 2
##########################################################################

##########################################################################
def find_optimum(L):
    """Finds the delay of the optimally buffered wire of L hops.

    Parameters
    ----------
    L : int
        Length of the wire in tile lengths.

    Returns
    -------
    float
        Delay of the optimal buffering
    """

    max_D0 = 5
    max_D1_over_D0 = 5

    rep_nos = [0] if not INSERT_REP else [2 ** i - 1 for i in range(int(math.log(L, 2)), -1, -1)]

    min_td = float("inf")
    for D0 in range(max_D0, 0, -1):
        for D1_over_D0 in range(max_D1_over_D0, 0, -1):
            D1 = D0 * D1_over_D0
            for rep_no in rep_nos:
                td = measure(D0, D1, L, rep_no)
                print D0, D1, rep_no, td
                if td < min_td:
                    min_td = td
                    best_D0 = D0
                    best_D1 = D1
                    best_rep_no = rep_no

    return min_td, best_D0, best_D1, best_rep_no
##########################################################################

##########################################################################
def print_results(buf_dict, td_dict, L):
    """Prints the results.

    Parameters
    ----------
    buf_dict : Dict[int, Tuple[int]]
        Dictionary of buffer sizes of the explored wire lengths.
    td_dict : Dict[int, float]
        Dictionary of the measured delays of the explored wire lengths.
    L : int
        Maximum (or fixed) wire length.

    Returns
    -------
    None
    """

    txt = "Maximum length = %d\n" % L

    txt += "Buffer sizes per length:\n"
    for l in sorted(buf_dict):
        txt += "%d: %d %d\n" % (l, buf_dict[l][0], buf_dict[l][1])

    txt += "Delays per length:\n"
    for l in sorted(td_dict):
        txt += "%d: %.2f ps\n" % (l, 1e12 * td_dict[l])

    filename = "buf_cache/%sK%dN%dT%s.log" % ('H' if HORIZONTAL else '', K, N, args.tech)
    with open(filename, "w") as outf:
        outf.write(txt)

    print(txt[:-1])
##########################################################################


if FIXED_LENGTH:
    td, best_D0, best_D1, rep_no = find_optimum(FIXED_LENGTH)
    td_dict = {FIXED_LENGTH : td}
    buf_dict = {FIXED_LENGTH : (best_D0, best_D1)}
    print_results(buf_dict, td_dict, FIXED_LENGTH)
    exit(0)

THR = 1

td_dict = {}
buf_dict = {}

td1, D01, D11, rep_no1 = find_optimum(1)
td_dict.update({1 : td1})
buf_dict.update({1 : (D01, D11)})

td2, D02, D12, rep_no2 = find_optimum(2)
td_dict.update({2 : td2})
buf_dict.update({2 : (D02, D12)})

composition = 2 * td1
single = td2

print composition, single

last_len = 1
if (composition - single) / composition * 100 < THR:
    print_results(buf_dict, td_dict, 1)
    exit(0)

last_len = 2
L = last_len

while((composition - single) / composition * 100 >= THR):
    last_len = L
    composition = 2 * single
    L *= 2
    single, best_D0, best_D1, best_rep_no = find_optimum(L)
    print L, single
    td_dict.update({L : single})
    buf_dict.update({L : (best_D0, best_D1)})

print_results(buf_dict, td_dict, last_len)

os.system("rm -f *.sp *.st0 *.ic0 *.mt0")
