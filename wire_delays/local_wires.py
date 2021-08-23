"""Measures the delay of a local (intracluster) wire (LAB line in Stratix architectures).

Parameters
----------
K : int
    LUT size.
N : int
    Cluster size.
tech : float
    Technology node (16, 7, 5, 4, 3.0, 3.1).
    3.0 corresponds to F3a in the paper and 3.1 to F3b.
density : float
    Crossbar density [0.0, 1.0]
use_my : Optional[bool], default = False
    Specifies that the long vertical wire (LAB line itself)
    should be raised from Mx to My.
insert_rep : Optional[bool], default = False
    Specifies that buffer insertion should be performed.
rep_at_lut_x : Optional[bool], default = False
    Specifies that repeater location on the x-axis should be at LUT output.
    Otherwise, they are optimistically placed inline with the vertical wire.
meas_cb : Optional[Tuple[float]], default = False
    Specifies that the delay from connection-block input to LUT input should be measured.
    Otherise, measurement is performed from the output of one LUT to the input of another,
    inside the same cluster. If the argument is passed, it should hold the x-offset of the
    vertical wire (LAB line), x-offset of the connection-block multiplexer, and the size
    of the connection-block multiplexer, respectively.
all_endpoints : Optional[bool], default = False
    Specifies that delays to all reachable (respecting crossbar sparsity) LUT inputs
    should be measured. Otherwise, only one delay, from the bottom-most LUT output, to
    a representative input in the middle of the cluster height is returned.
rebuffer : Optional[bool], default = False
    Specifies if the buffer sweep is to be performed.
    If not, the previously determined sizes are loaded from a cache.
wd : Optional[str], default = .
    Specify the working directory.

Returns
-------
Tuple[int]
    Buffer sizes
float
    Measured delay
""" 

import math
import os
import argparse
import copy
import networkx as nx
import sys
sys.path.insert(0,'..')

import setenv
import tech

parser = argparse.ArgumentParser()
parser.add_argument("--K")
parser.add_argument("--N")
parser.add_argument("--tech")
parser.add_argument("--density")
parser.add_argument("--use_my")
parser.add_argument("--insert_rep")
parser.add_argument("--rep_at_lut_x")
parser.add_argument("--meas_cb")
parser.add_argument("--all_endpoints")
parser.add_argument("--rebuffer")
parser.add_argument("--wd")
args = parser.parse_args()

K = int(args.K)
N = int(args.N)
density = float(args.density)

try:
    tech_node = int(args.tech)
except:
    tech_node = float(args.tech)

USE_MY = False
try:
    USE_MY = int(args.use_my)
except:
    pass

INSERT_REP = False
try:
    INSERT_REP = int(args.insert_rep)
except:
    pass

REP_AT_LUT_X = False
try:
    REP_AT_LUT_X = int(args.rep_at_lut_x)
except:
    pass

node_index = tech.nodes.index(tech_node)
node_device_index = tech.node_names.index(int(tech_node))

MxR = tech.MxR[node_index]
MxC = tech.MxC[node_index] * 1e-15

MyR = tech.MyR[node_index]
MyC = tech.MyC[node_index] * 1e-15
GP = tech.GP[node_device_index]
FP = tech.FP[node_device_index]
GL = tech.GL[node_device_index]
vdd = tech.vdd[node_device_index]

Rvia = tech.MxMx_via[node_index]
if USE_MY:
    Rvia = tech.stacked_via[node_index]

D0 = 1
D1 = 1
DO_REBUFFER = False
buf_log_filename = "buf_cache/K%dN%dD%.2fR%dX%dY%dT%s.log"\
                 % (K, N, density, 1 if INSERT_REP else 0,\
                    1 if REP_AT_LUT_X else 0, 1 if USE_MY else 0, args.tech)
try:
    DO_REBUFFER = int(args.rebuffer)
except:
    try:
        with open(buf_log_filename, "r") as inf:
            lines = inf.readlines()
            D0 = int(lines[-2].split()[0])
            D1 = int(lines[-2].split()[1])
    except:
        DO_REBUFFER = True

SIM_ALL = False
try:
    SIM_ALL = int(args.all_endpoints)
except:
    pass

Cw = 1.0 / 1000.0 * MxC
Rw = 1.0 / 1000.0 * MxR
Cwy = 1.0 / 1000.0 * MyC
Rwy = 1.0 / 1000.0 * MyR

spice_model_path = "\"%s/" % os.path.abspath("../spice_models/") + "%dnm.l\" %dNM_FINFET_HP\n"


get_w = lambda row_cnt : 21 + row_cnt
#Width in fin pitches

get_h = lambda row_cnt, col_cnt : 2 * max(row_cnt, col_cnt)
#Height in gate pitches

WD = None
try:
    WD = args.wd
except:
    pass

##########################################################################
def get_mux_dimensions(I, max_height):
    """Finds the row and column count for the given mux.

    Parameters
    ----------
    I : int
        Input count.
    max_height : int
        Maximum height in gate pitches.

    Returns
    -------
    Tuple[int]
        Row and column count.
    """
       
    col_cnt = int(min(max_height / 2.0, math.ceil(I ** 0.5)))
    row_cnt = int(math.ceil(I / float(col_cnt)))

    return row_cnt, col_cnt
##########################################################################

p = 0.8
#Desired Rent's exponent.

cluster_inputs = int(math.ceil((K * (N ** p)) / N) * N)
crossbar_mux_size = density * (N + cluster_inputs)

lut_height = 2 ** (K - 4) * 2 * (4 + 2)
crossbar_mux_rcnt, crossbar_mux_ccnt = get_mux_dimensions(crossbar_mux_size, lut_height)

crossbar_mux_phys_width = get_w(crossbar_mux_rcnt)
crossbar_mux_phys_height = get_h(crossbar_mux_rcnt, crossbar_mux_ccnt)

cur_height = 0
cur_x_offset = -1
mux_input_coords = []
for m in range(0, K):
    coords = ((0.5 + cur_x_offset) * crossbar_mux_phys_width, cur_height * crossbar_mux_phys_height)
    mux_input_coords.append(coords)
    if (cur_height + 1) * crossbar_mux_phys_height > lut_height:
        cur_height = 0
        cur_x_offset -= 1
    else:
        cur_height += 1

if density != 1:
    #Now sparsify the crossbar:
    to_remove = int(math.floor(K * density))
    rem_space = int(math.ceil(1 / (1.0 - density)))
    rem_list = []
    for i in range(0, K):
        if len(rem_list) < to_remove and i % rem_space == 0:
            rem_list.append(mux_input_coords[i])
    for mux in rem_list:
        mux_input_coords.remove(mux)

min_x = min(c[0] for c in mux_input_coords)
max_x = max(c[0] for c in mux_input_coords)

lut_x = 2 ** 4 * 10 / 2.0
lut_y = lut_height / (2 ** (K - 4) / 2.0)

meas_cb = False
try:
    wire_x = float(args.meas_cb.split()[0])
    cb_x = float(args.meas_cb.split()[1])
    cb_size = int(args.meas_cb.split()[2])
    vertical_x = 0.5 * (max_x + cb_x)
    meas_cb = True
except:
    vertical_x = 0.5 * (min_x + lut_x)

all_mux_coords = copy.deepcopy(mux_input_coords)

for i in range(1, N):
    for m in mux_input_coords:
        all_mux_coords.append((m[0], i * lut_height + m[1])) 

net = nx.Graph()

if meas_cb:
    net.add_node("wire_t", coords = (wire_x, lut_y), mux = False)
    net.add_node('s', coords = (cb_x, lut_y), mux = True, size = cb_size)
    net.add_edge("wire_t", "s", spine = False)
else: 
    net.add_node('s', coords = (lut_x, lut_y), mux = False)

ncnt = 0
prev_cnt = None
joint_dict = {}
for y in sorted(set([m[1] for m in all_mux_coords] + [lut_y])):
    net.add_node(ncnt, coords = (vertical_x, y), mux = False)
    if prev_cnt is not None:
        net.add_edge(prev_cnt, ncnt, spine = True)
    if y == lut_y:
        net.add_edge('s', ncnt)
    joint_dict.update({y : ncnt})
    prev_cnt = ncnt
    ncnt += 1

muxes_per_height = {}
for m in all_mux_coords:
    try:
        muxes_per_height[m[1]].append(m[0])
    except:
        muxes_per_height.update({m[1] : [m[0]]})

ncnt += 1
for h in muxes_per_height:
    muxes_per_height[h].sort()
    prev_cnt = None
    prev_x = -1
    joint_added = False
    for m in muxes_per_height[h]:
        net.add_node(ncnt, coords = (m, h), mux = True)
        if (m > vertical_x and prev_x < vertical_x):
            if prev_cnt != None:
                net.add_edge(prev_cnt, joint_dict[h])
            net.add_edge(joint_dict[h], ncnt)
            joint_added = True
        elif prev_cnt is not None:
            net.add_edge(prev_cnt, ncnt)
        prev_cnt = ncnt
        prev_x = m
        ncnt += 1
    if not joint_added:
        net.add_edge(prev_cnt, joint_dict[h])

#seen = set()
#for u in sorted(net.nodes(), key = lambda n : net.node[n]["coords"]):
#    for v in sorted(net.nodes(), key = lambda m : net.node[m]["coords"]):
#        if net.has_edge(u, v) and not tuple(sorted([u, v])) in seen:
#            print u, v, net.node[u]["coords"], net.node[v]["coords"]
#            seen.add(tuple(sorted([u, v])))
#exit(0)

#Compute the multiplexer states
##########################################################################
on_fanout = 2
#Corresponds roughly to the typical fanout of a net in the MCNC circuits.

max_fanout = K * N * density

prob_partial = 1.0 / crossbar_mux_ccnt
#Probability that the first stage of the mux is turned on.

partial_fanout = int(math.ceil(prob_partial * max_fanout))

mux_nodes = [u for u, attrs in net.nodes(data = True) if attrs["mux"] and u != 's']
mux_nodes.sort(key = lambda m : (net.node[m]["coords"][1], net.node[m]["coords"][0]))

partial_space = len(mux_nodes) / partial_fanout
for i, m in enumerate(mux_nodes):
    if i % partial_space == 0:
        net.node[m]["state"] = "partial"
    else:
        net.node[m]["state"] = "off"

net.node[mux_nodes[-1]]["state"] = "on"
target_mux = mux_nodes[len(mux_nodes) / 2] if not INSERT_REP else mux_nodes[-1]
print net.node[target_mux]["coords"]
net.node[target_mux]["state"] = "on"
net.add_node('t', coords = (0, net.node[target_mux]["coords"][1]), mux = False)
net.add_edge(target_mux, 't')
if meas_cb:
    net.node['s']["state"] = "on"

##########################################################################
def conv_nx_to_spice(net, invert_trig = False):
    """Converts the net to a spice netlist.

    Parameters
    ----------
    net : nx.Graph
        The net graph.
    invert_trig : Optional[bool], default = False
        Specifies if the trigger signal should be inverted or not.

    Returns
    -------
    str
        SPICE netlist description.
    """

    txt = ".TITLE LOCAL_WIRE_MEAS\n\n"
    txt += ".LIB %s\n" % (spice_model_path % (int(tech_node), int(tech_node)))
    txt += ".TRAN 1p 16n\n.OPTIONS BRIEF=1\n\n"

    txt += ".PARAM Cw=%g\n" % Cw
    txt += ".PARAM Rw=%g\n" % Rw
    txt += ".PARAM Cwy=%g\n" % Cwy
    txt += ".PARAM Rwy=%g\n" % Rwy
    txt += ".PARAM Rvia=%g\n\n" % Rvia
    txt += ".PARAM D0=%d\n" % D0
    txt += ".PARAM D1=%d\n\n" % D1
    txt += ".PARAM GL=%dn\n" % GL
    txt += ".PARAM supply_v=%.2f\n\n" % vdd

    #------------------------------------------------------------------------#
    def add_wire_subckt(stages = 1, my_wire = False):
        """Adds a single tile length wire subcircuit.

        Parameters
        ----------
        stages : Optional[int], default = 1
            Number of Pi stages in the wire model.
        my_wire : Optional[bool], default = false
            Specifies if the wire should be traced on My layer.            

        Returns
        -------
        str
            Subcircuit description.
        """

        txt = ".SUBCKT %swire n_in n_out l=1\n" % ("my_" if my_wire else '')

        stage_template = "C%d_1 %s gnd" +  " C='l*Cw%s/(2*%d)'\n" % ('y' if my_wire else '', stages)
        stage_template += "R%d_1 %s %s" + " R='l*Rw%s/%d'\n" % ('y' if my_wire else '', stages)
        stage_template += "C%d_2 %s gnd" + " C='l*Cw%s/(2*%d)'\n" % ('y' if my_wire else '', stages)


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

    txt += add_wire_subckt(1)
    txt += add_wire_subckt(1, my_wire = True)
    #NOTE: Because we break the routing tee into many small segments,
    #due to the numerous connected multiplexers, one RC stage provides
    #a reasonable approximation.

    txt += ".SUBCKT buf n_in n_out vdd STRENGTH0=D0 STRENGTH1=D1\n"
    txt += "Rin n_in n_in_via R=Rvia\n"
    txt += "MN1 n_mid n_in_via gnd gnd nmos L=GL nfin=STRENGTH0\n"
    txt += "MP1 n_mid n_in_via vdd vdd pmos L=GL nfin=STRENGTH0\n"
    txt += "MN2 n_out_pre_via n_mid gnd gnd nmos L=GL nfin=STRENGTH1\n"
    txt += "MP2 n_out_pre_via n_mid vdd vdd pmos L=GL nfin=STRENGTH1\n"
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
        txt += "R_in_via n_in n_in_post_via R=Rvia\n" 
        for r in range(0, row_cnt):
            row_node = "n_r_%d" % r
            for c in range(0, col_cnt):
                i = "n_in_post_via" if c == 0  and r == 0 else "n_dummy_%d" % (r * col_cnt + c)
                tg_model = "tg_on" if c == 0 and state != "off" else "tg_off"
                txt += "Xtg_%d %s %s vdd %s\n" % (ind, i, row_node, tg_model)
                ind += 1
            tg_model = "tg_on" if r == 0 and state == "on" else "tg_off"
            txt += "Xtg_%d %s %s vdd %s\n" % (ind, row_node, "n_out", tg_model)
            ind += 1

        return txt + ".ENDS\n\n"
    #------------------------------------------------------------------------#

    txt += add_mux_subckt(crossbar_mux_rcnt, crossbar_mux_ccnt, "on")
    txt += add_mux_subckt(crossbar_mux_rcnt, crossbar_mux_ccnt, "off")
    txt += add_mux_subckt(crossbar_mux_rcnt, crossbar_mux_ccnt, "partial")

    if meas_cb:
        cb_rcnt, cb_ccnt = get_mux_dimensions(cb_size, lut_height)
        if cb_rcnt != crossbar_mux_rcnt or cb_ccnt != crossbar_mux_ccnt:
            txt += add_mux_subckt(cb_rcnt, cb_ccnt, "on")

    txt += "Vps vdd gnd supply_v\n"
    if meas_cb:
        txt += "Vin n_in gnd PULSE (0 supply_v 0 0 0 2n 4n)\n\n"
        txt += "Rvia_in_cb n_in n_wire_t R=%f\n\n" % tech.stacked_via[node_index]
    else:
        txt += "Vin n_pre_mux gnd PULSE (0 supply_v 0 0 0 2n 4n)\n\n"
        txt += "Xtg_drv_mux_2_1 n_pre_mux n_in vdd tg_on\n"
        txt += "Xbuf_drv n_in n_s vdd buf\n\n"

    ecnt = -1
    mux_cnt = -1
    for u, v in net.edges():
        ecnt += 1
        u_coords = net.node[u]["coords"]
        v_coords = net.node[v]["coords"]

        if u_coords[0] == v_coords[0]:
            L = abs(u_coords[1] - v_coords[1]) * GP
        elif u_coords[1] == v_coords[1]:
            L = abs(u_coords[0] - v_coords[0]) * FP
        else:
            print "Diagonal wire!"
            raise ValueError

        is_spine = net[u][v].get("spine", False)
        model = "my_wire" if is_spine and USE_MY else "wire"

        if meas_cb and tuple(sorted([u, v])) == ('s', "wire_t"):
            txt += "Xwire%d n_wire_t n_s_pre_mux %s L=%.2f\n" % (ecnt, model, L)
        else:
            if u != 't' and v != 't':
                if INSERT_REP and net[u][v].get("buf", False):
                    if REP_AT_LUT_X:
                        L_rep_to_wire = abs(lut_x - vertical_x) * FP
                        print "inserting between", u, v
                        txt += "Xwire%d_wire_to_rep n_%s n_%s_wire_to_rep %s L=%.2f\n" % (ecnt, str(u), str(u), model, L_rep_to_wire)
                        txt += "Xrep%d n_%s_wire_to_rep n_%s_rep_out vdd rep\n" % (ecnt, str(u), str(u))
                        txt += "Xwire%d_rep_to_wire n_%s_rep_out n_%s_rep_to_wire %s L=%.2f\n" % (ecnt, str(u), str(u), model, L_rep_to_wire)
                        txt += "Xwire%d n_%s_rep_to_wire n_%s %s L=%.2f\n" % (ecnt, str(u), str(v), model, L)
                    else:
                        txt += "Xrep%d n_%s n_%s_rep_out vdd rep\n" % (ecnt, str(u), str(u))
                        txt += "Xwire%d n_%s_rep_out n_%s %s L=%.2f\n" % (ecnt, str(u), str(v), model, L)
                else:
                    txt += "Xwire%d n_%s n_%s %s L=%.2f\n" % (ecnt, str(u), str(v), model, L)
            else:
                txt += "Xwire%d n_t_drv n_t %s L=%.2f\n" % (ecnt, model, L)
                txt += "Xbuf_load n_t n_out vdd buf STRENGTH0=1 STRENGTH1=2\n\n"
    
        mux_node = None
        if net.node[v]["mux"]:
            mux_node = v
        elif net.node[u]["mux"]:
            mux_node = u
        if mux_node is not None:
            if meas_cb and mux_node == 's':
                if u == "wire_t" or v == "wire_t":
                    txt += "Xmux_cb n_s_pre_mux n_cb_buf_in vdd mux_%d_%d_on\n" % (cb_rcnt, cb_ccnt)
                    txt += "Xbuf_drv n_cb_buf_in n_s vdd buf\n\n"
            else:
                mux_cnt += 1
                n_out = "n_dummy_mux_out%d" % mux_cnt
                if 't' in net[mux_node]:
                    n_out = "n_t_drv"
                txt += "Xmux%d n_%s %s vdd mux_%d_%d_%s\n\n" % (mux_cnt, str(mux_node), n_out,\
                       crossbar_mux_rcnt, crossbar_mux_ccnt, net.node[mux_node]["state"])

    txt += ".MEASURE tfall TRIG V(n_in) VAL='supply_v/2' %s=2\n" % ("RISE" if invert_trig else "FALL")
    txt += "+                  TARG V(n_t) VAL supply_v/2 FALL=2\n\n"
    txt += ".MEASURE trise TRIG V(n_in) VAL='supply_v/2' %s=2\n" % ("FALL" if invert_trig else "RISE")
    txt += "+                  TARG V(n_t) VAL supply_v/2 RISE=2\n\n"

    txt += ".END"

    return txt
##########################################################################

netlist_filename = "sim_local_K%dN%dT%s.sp" % (K, N, str(tech_node))
hspice_dump = "hspice_K%dN%dT%s.dump" % (K, N, str(tech_node))

##########################################################################
def measure(invert_trig = False):
    """Calls HSPICE to obtain the delay.
    
    Parameters
    ----------
    invert_trig : Optional[bool], default = False
        Specifies if the trigger signal should be inverted or not.

    Returns
    -------
    float
        Delay.
    """

    to_write = conv_nx_to_spice(net, invert_trig = invert_trig)

    wd = os.getcwd()
    if WD is not None:
        os.chdir(WD)

    with open(netlist_filename, "w") as outf:
        outf.write(to_write)
   
    hspice_call = os.environ["HSPICE"] + " %s > %s" % (netlist_filename, hspice_dump)
    os.system(hspice_call)
   
    scale_dict = {'f' : 1e-15, 'p' : 1e-12, 'n' : 1e-9}
   
    with open(hspice_dump, "r") as inf:
        lines = inf.readlines()

    os.system("rm " + hspice_dump)
    os.system("rm -f *.sp *.st0 *.ic0 *.mt0")
    
    if WD is not None:
        os.chdir(wd)

    for line in lines:
        if "tfall=" in line:
            tfall = float(line.split()[1][:-1]) * scale_dict[line.split()[1][-1]]
        elif "trise=" in line:
            trise = float(line.split()[1][:-1]) * scale_dict[line.split()[1][-1]]
    if trise < 0 or tfall < 0:
        print("Negative time!")
        if not DO_REBUFFER:
            raise ValueError
        else:
            return float('inf')
   
    return (trise + tfall) / 2
##########################################################################

##########################################################################
def rebuffer(invert_trig = False):
    """Finds another optimal buffer, given that we now know the precise load.

    Parameters
    ----------
    invert_trig : Optional[bool], default = False
        Specifies if the trigger signal should be inverted or not.

    Returns
    -------
    Tuple[int, float]
        Buffer size and delay.
    """

    max_D0 = 5
    max_D1_over_D0 = 5
    
    min_td = float("inf")
    for d0 in range(max_D0, 0, -1):
        for d1_over_D0 in range(max_D1_over_D0, 0, -1):
            d1 = d0 * d1_over_D0
            global D0
            D0 = d0
            global D1
            D1 = d1
            td = measure(invert_trig = invert_trig)
            print D0, D1, td
            if td > 0 and td < min_td:
                min_td = td
                best_D0 = D0
                best_D1 = D1

    return min_td, best_D0, best_D1
##########################################################################

##########################################################################
def sim_all_endpoints():
    """Finds the delays to all endpoints, without rebuffering.

    Parameters
    ----------
    None

    Returns
    -------
    List[float]
        A sorted list of all delays.
    """

    tds = []

    for target_mux in mux_nodes:
        net.remove_node('t')
        for i, m in enumerate(mux_nodes):
            if i % partial_space == 0:
                net.node[m]["state"] = "partial"
            else:
                net.node[m]["state"] = "off"

        net.node[mux_nodes[-1]]["state"] = "on"
        print net.node[target_mux]["coords"]
        net.node[target_mux]["state"] = "on"
        net.add_node('t', coords = (0, net.node[target_mux]["coords"][1]), mux = False)
        net.add_edge(target_mux, 't')

        sp = nx.shortest_path(net, 's', target_mux)
        inv_cnt = 0
        for i, u in enumerate(sp[:-1]):
            if net[u][sp[i + 1]].get("buf", False):
                inv_cnt += 1

        td = measure(invert_trig = inv_cnt % 2)
        tds.append(td)

    return sorted(tds)
##########################################################################

##########################################################################
def insert_reps(fixed_comb = None, sim_all = False):
    """Inserts repeaters in the spine.

    Parameters
    ----------
    fixed_comb : Optional[List[int]], default = None
        A single predetermined insertion combination to be used.
    sim_all : Optional[bool]
        Determines if all endpoint delays should be found.

    Returns
    -------
    float
        Delay with quasi optimal insertion.
    int
        D0
    int
        D1
    List[int]
        Repeater presence list
    """

    if not INSERT_REP:
        return
    
    spine = [(u, v) for u, v, attrs in net.edges(data = True)
             if attrs.get("spine", False) and net.node[u]["coords"][1] > lut_y\
                                          and net.node[v]["coords"][1] > lut_y]
    spine.sort(key = lambda e : min(net.node[e[0]]["coords"][1], net.node[e[1]]["coords"][1]))

    if fixed_comb:
        combs = [fixed_comb]
    else:
        combs = []
        for space in range(1, len(spine) / 2 + 1):
            combs.append([0])
            for i in range(1, len(spine)):
                if i % space == 0:
                    combs[-1].append(1)
                else:
                    combs[-1].append(0)
        combs.append([0 for i in range(0, len(spine))])

    all_min_td = float("inf")
    all_best_D0 = 0
    all_best_D1 = 0
    all_best_comb = []
    
    for comb in reversed(combs):
        print comb
        invert_trig = len([i for i in comb if i]) % 2
        for i, e in enumerate(spine):
            u, v = e
            net[u][v]["buf"] = comb[i]
        min_td, best_D0, best_D1 = rebuffer(invert_trig = invert_trig)
        print best_D0, best_D1
        print min_td
        if min_td < all_min_td:
            all_min_td = min_td
            all_best_D0 = best_D0
            all_best_D1 = best_D1
            all_best_comb = comb

    if sim_all:
        comb = all_best_comb
        invert_trig = len([i for i in comb if i]) % 2
        for i, e in enumerate(spine):
            u, v = e
            net[u][v]["buf"] = comb[i]

        global D0
        D0 = all_best_D0
        global D1
        D1 = all_best_D1
        
        return sim_all_endpoints(), all_best_D0, all_best_D1, all_best_comb

    return all_min_td, all_best_D0, all_best_D1, all_best_comb
##########################################################################

##########################################################################
def log_results(d0, d1, td, filename):
    """Logs the buffering results.

    Parameters
    ----------
    d0 : int
        Size of the first buffer stage.
    d1 : int
        Size of the second buffer stage.
    td : float
        Delay in picoseconds.
    filename : str
        Name of the log file.

    Returns
    -------
    None
    """

    with open(filename, "w") as outf:
        outf.write("%d %d\n" % (d0, d1))
        outf.write("%.2f" % td)
##########################################################################

if INSERT_REP:
    min_td, d0, d1, comb = insert_reps()
    min_td *= 1e12
    print(min_td)
    #NOTE: No need to cache the buffers here, as these are one-time experiments.
    if SIM_ALL:
        global D0
        D0 = d0
        global D1
        D1 = d1
        print(insert_reps(fixed_comb = comb, sim_all = True)[0])   
elif DO_REBUFFER:
    min_td, best_D0, best_D1 = rebuffer()
    min_td *= 1e12
    print(min_td)
    log_results(best_D0, best_D1, min_td, buf_log_filename)
    if SIM_ALL:
        global D0
        D0 = best_D0
        global D1
        D1 = best_D1
        print(sim_all_endpoints())
elif SIM_ALL:
    print(sim_all_endpoints())
else:
    print(measure())
    #We read the sizes from the cache, hence no need to log again.

os.system("rm -f *.sp *.st0 *.ic0 *.mt0")
