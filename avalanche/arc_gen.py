"""Generates the architecture xml description and the routing-resource-graph.

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
grid_w : int
    Width of the FPGA grid.
grid_h : int
    Height of the FPGA grid.
wire_file : str
    Name of the file containing the channel composition description.
arc_name : str
    Name of the exported architecture.
change_grid_dimensions : Optional[str], default = None
    Generates another architecture from the specified file, merely
    changing the grid dimensions.
physical_square : Optional[bool], default = True
    Changes logical grid dimensions after tile dimensions are known,
    so as to obtain a grid that is closer to a square, maintaining the tile count.
only_pad : Optional[bool], default = False
    Instructs the script to only complete the LEN-1 wire
    padding to obtain the tile area statistics, without
    spicing the wires or generating the architecture files.
import_padding : Optional[str], default = None
    Specifies a padding log from which to inherit padding.
    It is assumed that the log has been created for N=operating_point
    (defined later in the code, usually 8).
make_sb_clique : Optional[bool], default = False
    Specifies that dummy, high-cost edges should be created
    between all switch-block inputs and outputs, broken by nodes
    with an assigned cost, amenable to avalanche updates.
first_clique_iter : Optional[bool], default = False
    Specifies that this is the first iteration of clique search and that the
    VPR routing log should not be parsed.
no_stored_unused_edge_increase : Optional[bool], default = False
    Instructs the script not to increase the unused iteration count for the potential
    edges. Useful when calling the resize routines during wafer generation.
robustness_level : Optional[int], default = 2
    Sets the measurement robustness level:
    0 -> measurement is performed on the first found wire of the given type.
    1 -> a single representative wire is picked for measurement,
         based on the multiplexer width and fanout (median product).
    2 -> average over all wires is determined for each type, but the
         target is fixed at the fanout at average distance from the wire.
    3 -> all targets are measured for all wires. This usually brings little
         change (a fraction of a ps at 4nm, for example).
    Note that in all cases, the horizontal wire is assumed to be at the
    middle height of its LUT and vertical wire in the middle of the tile.
    What changes is the location of the driving multiplexer.
insert_empty_ring : Optional[bool], default = False
    Inserts a ring of empty clusters adjacent to the I/Os, to cover the fringe effects.
    The grid dimensions are increased accordingly.
greedy_switch_search : Optional[bool], default = False
    Instructs the script to configure avalanche parameters for pure greedy search.
    (starts at zero costs, no updates, no delays).
load_pattern_from_file : Optional[str], default = None
    Instructs the script to load a (partial) SB-switch-pattern from a file.
pattern_is_init_only : Optional[bool], default = False
    Tells the script that the loaded pattern is only for initialization.
pick_top_only : Optional[bool], default = False
    Instructs the script to only pick the most used switch, breaking ties using delay from the previous architecture.
adoption_threshold : Optional[float], default = 1.2
    Threshold for utilization drop above which all switches are adopted into the pattern.
load_mux_stack_order : Optional[str], default = None
    Specifies that the multiplexer stacking order should be read from a file.

Returns
-------
None
""" 

import time
import os
import copy
import ast
import networkx as nx
import math
import argparse
import numpy as np
import random
import sys
sys.path.insert(0,'..')

import setenv
import tech

import parse_routing

parser = argparse.ArgumentParser()
parser.add_argument("--K")
parser.add_argument("--N")
parser.add_argument("--tech")
parser.add_argument("--density")
parser.add_argument("--grid_w")
parser.add_argument("--grid_h")
parser.add_argument("--wire_file")
parser.add_argument("--arc_name")
parser.add_argument("--change_grid_dimensions")
parser.add_argument("--physical_square")
parser.add_argument("--only_pad")
parser.add_argument("--import_padding")
parser.add_argument("--make_sb_clique")
parser.add_argument("--first_clique_iter")
parser.add_argument("--no_stored_unused_edge_increase")
parser.add_argument("--robustness_level")
parser.add_argument("--insert_empty_ring")
parser.add_argument("--greedy_switch_search")
parser.add_argument("--load_pattern_from_file")
parser.add_argument("--pattern_is_init_only")
parser.add_argument("--pick_top_only")
parser.add_argument("--adoption_threshold")
parser.add_argument("--load_mux_stack_order")

args = parser.parse_args()
K = int(args.K)
N = int(args.N)
density = float(args.density)

try:
    tech_node = int(args.tech)
except:
    tech_node = float(args.tech)

node_index = tech.nodes.index(tech_node)
node_device_index = tech.node_names.index(int(tech_node))

MxR = tech.MxR[node_index]
MxC = tech.MxC[node_index] * 1e-15

MyR = tech.MyR[node_index]
MyC = tech.MyC[node_index] * 1e-15
MyP = tech.MyP[node_index]

GP = tech.GP[node_device_index]
FP = tech.FP[node_device_index]
GL = tech.GL[node_device_index]
vdd = tech.vdd[node_device_index]

ABS_MIN_HEIGHT = ABS_MIN_WIDTH = 7
#NOTE: VPR computes the router lookahead by using the coordinates (3, 3)--(5, 5) as starts.
#Hence, if we have a grid smaller than that, it is going to crash.

grid_w = int(args.grid_w)
if grid_w < ABS_MIN_WIDTH:
    grid_w = ABS_MIN_WIDTH
grid_h = int(args.grid_h)
if grid_h < ABS_MIN_HEIGHT:
    grid_h = ABS_MIN_HEIGHT

default_delay = float("inf")
#Proxy delay replaced by measurements.
operating_point = 8
#Cluster size for which the optimum segmentation was found.

tap_phi = 1
#Spacing between taps.
tap_M = max(1, operating_point / N)
#Total number of taps.

Rvia = tech.stacked_via[node_index]

local_driver = (1, 1)
local_buf_filename = "../wire_delays/buf_cache/K%dN%dD%.2fR%dX%dY%dT%s.log"\
                   % (K, N, density, 0, 0, 0, args.tech)
with open(local_buf_filename, "r") as inf:
    lines = inf.readlines()
    local_driver = (int(lines[-2].split()[0]), int(lines[-2].split()[1]))
    default_cb_delay = 1e-12 * float(lines[-1].strip())
    local_feedback_delay = default_cb_delay

Cw = 1.0 / 1000.0 * MxC
Rw = 1.0 / 1000.0 * MxR
Cwy = 1.0 / 1000.0 * MyC
Rwy = 1.0 / 1000.0 * MyR

io_crop = 1
#Number of columns/rows on the periphery belonging to I/O.

mult_col = {"start" : (-1, -1), "freq" : 0, "height" : 4}
mem_col = {"start" : (-1, -1), "freq" : 0, "height" : 6}

block_ids = {"empty" : 0, "io" : 1, "clb" : 2, "empty_clb" : 3, "mult" : 4, "mem" : 5} 
mux_ids = {"__vpr_delayless_switch__" : 0, "cb" : 1, "sb" : 2}
seg_ids = {}

cut_corners = True
#Specifies if the corners should be left empty.
TOP_BOTTOM_IO = False
#Specifies that the I/O pads should be located only at top and bottom of the chip.
#This is not uncommon (Versal, Agilex) and allows for keeping the same I/O capacity
#accross different cluster sizes.
IO_CAPACITY = 8 if TOP_BOTTOM_IO else N
#Capacity of a single I/O pad.

indent = "    "
#Default indentation.

p = 0.8
#Target Rent's exponent for determining cluster input count.
cluster_inputs = int(math.ceil((K * (N ** p)) / N) * N)
#cluster_inputs = 4 * N
#NOTE: To minimize noise, it is wiser to always put 4 X N as the number of inputs.
#Note that this holds for N in {4, 8, 16}, but only for N=2 it is smaller. On the
#other hand, the Betz and Rose formula would give 3 X 3 = 9 inputs for N=2, which is
#almost equal to the 4 X 2 = 8.
O = 2
#Outputs per LUT
crossbar_mux_size = density * (O * N + cluster_inputs)

lut_height = 2 ** (K - 4) * 2 * (4 + 2)
lut_width = 2 ** 4 * 10

H = []
V = []
 
K6N8_LUT4 = 8 * 4
KN_LUT4 = N * 2 ** (K - 4)
print("Scaling ratio from K6N8: %d" % int(math.ceil((K6N8_LUT4 / float(KN_LUT4)))))

with open(args.wire_file, "r") as inf:
    lines = inf.readlines()
for line in lines:
    if line[0] == 'H':
        H.append((int(line.split()[1]), int(line.split()[2])))
    elif line[0] == 'V':
        L = int(line.split()[1])
        L = int(math.ceil((L * K6N8_LUT4 / float(KN_LUT4))))
        if L < 1:
            L = 1
        V.append((L, int(line.split()[2])))

ONLY_PAD = False
try:
    ONLY_PAD = int(args.only_pad)
except:
    pass

PHYSICAL_SQUARE = False
try:
    PHYSICAL_SQUARE = int(args.physical_square)
except:
    pass

MAKE_SB_CLIQUE = False
try:
    MAKE_SB_CLIQUE = int(args.make_sb_clique)
except:
    pass

FINALIZE_CLIQUE = False
FIRST_CLIQUE_ITER = False
if MAKE_SB_CLIQUE:
    try:
        FIRST_CLIQUE_ITER = int(args.first_clique_iter)
    except:
        pass

PATTERN_IS_INIT_ONLY = False
try:
    PATTERN_IS_INIT_ONLY = int(args.pattern_is_init_only)
except:
    pass

PICK_TOP_ONLY = False
try:
    PICK_TOP_ONLY = int(args.pick_top_only)
except:
    pass

ADOPTION_THRESHOLD = 1.2
try:
    ADOPTION_THRESHOLD = float(args.adoption_threshold)
except:
    pass

NO_STORED_UNUSED_EDGE_INCREASE = False
try:
    NO_STORED_UNUSED_EDGE_INCREASE = int(args.no_stored_unused_edge_increase)
except:
    pass

print "NO_STORED_UNUSED_EDGE_INCREASE", NO_STORED_UNUSED_EDGE_INCREASE

NEUTRAL_BLE = N / 2
MAX_BLE_SPAN = 1
#Specifies the maximum BLE distance between the connected nodes, when constructing the clique.
assert ((NEUTRAL_BLE - MAX_BLE_SPAN) >= 0) and ((NEUTRAL_BLE + MAX_BLE_SPAN) <= (N - 1)),\
       "Neutral BLE and maximum BLE span incompatible!"

ROBUSTNESS_LEVEL = 2
try:
    ROBUSTNESS_LEVEL = int(args.robustness_level)
except:
    pass

INSERT_EMPTY_RING = False
try:
    INSERT_EMPTY_RING = int(args.insert_empty_ring)
except:
    pass

GREEDY_SWITCH_SEARCH = False
try:
    GREEDY_SWITCH_SEARCH = int(args.greedy_switch_search)
except:
    pass

SB_MUX_ORDER = None
##########################################################################
def read_buffer_cache(tech_name):
    """Reads the buffer sizes from the cache.
    
    Parameters
    ----------
    tech_name : str
        Name of the technology node (args.tech).

    Returns
    -------
    Dict[str, Tuple[int]]
        A dictionary of buffer sizes per wire type.
    """

    H_drivers = {}
    V_drivers = {}

    filename = "../wire_delays/buf_cache/%s" + "K6N8T%s.log" % tech_name
    
    try:
        with open(filename % 'H', "r") as inf:
            lines = inf.readlines()
    except:
        print("Buffer cache for technology %s not found (%s).\n\
               Please run ../wire_delays/global_wire.py to produce it."\
             % (tech_name, filename % 'H'))
        raise FileNotFoundError

    rd = False
    for line in lines:
        if not rd:
            if line.startswith("Buffer sizes per length"):
                rd = True
            continue
        if line.startswith("Delays per length"):
            break
        L = int(line.split(':')[0])
        d0 = int(line.split()[1])
        d1 = int(line.split()[2])
        H_drivers.update({L : (d0, d1)})

    try:
        with open(filename % '', "r") as inf:
            lines = inf.readlines()
    except:
        print("Buffer cache for technology %s not found (%s).\n\
               Please run ../wire_delays/global_wire.py to produce it."\
             % (tech_name, filename % ''))
        raise FileNotFoundError

    rd = False
    for line in lines:
        if not rd:
            if line.startswith("Buffer sizes per length"):
                rd = True
            continue
        if line.startswith("Delays per length"):
            break
        L = int(line.split(':')[0])
        L = max(1, (L * K6N8_LUT4) / KN_LUT4)
        d0 = int(line.split()[1])
        d1 = int(line.split()[2])
        V_drivers.update({L : (d0, d1)})


    #FIXME: Run a sweep here too.
    H_drivers.update({6 : H_drivers[4]})

    return H_drivers, V_drivers
##########################################################################
 
H_drivers, V_drivers = read_buffer_cache(args.tech)
spice_model_path = "\"%s/" % os.path.abspath("../spice_models/") + "%dnm.l\" %dNM_FINFET_HP\n"

COMPRESS_RR = True
#Instructs the script to compress the produced RR-graph to save storage space
#that can otherwise become problematic over many architectures.

HUMAN_READABLE = False
#Specifies if the node names and edges should be integer-based in the exported RR-graph
#(Needed by VPR), or strings that correspond to wire and pin identifiers.

D0 = D1 = 1
#Default driver sizes.

#Switch pattern parameters:
###########################

DISJOINT_SB = True
#Specifies a disjoint pattern instead of the default full.
#(always bounded to pins at the same LUT height)
#By disjoint, we mean that each wire drives only one wire of each group,
#but does drive all groups. In case all wires at one LUT height are of different
#type, full and disjoint patterns are the same.
DISJOINT_CB = True
#Specifies a disjoint pattern instead of the default full.
#(always bounded to pins at the same LUT height)
MAX_LUT_FANOUT = float('inf')
#Maximum LUT fanout.
ADD_LEN_1_TWISTS = True
#Allow connections among LEN-1 wires accross neighboring LUT heights.
ONLY_CONTINUATION_TWISTS =  True
#Perform cross-LUT twisting only on continuation wires.
CUT_CROSS_CLB_TWISTS = True
#Prevent the connections between neighboring LUT-height LEN-1 wires
#from crossing the tile boundaries. In principle, this can be done,
#to further improve routability, but causes VPR to deem the RR-graphs illegal.

SEPARATE_TAPS = False
#Specifies if the taps should be separated as different nodes in the RR-graph.
#This allows for less conservative modeling of tap delays but may compromise
#router lookahead effectiveness and cause routability issues.

get_index = lambda wire : int(wire.split("_tap")[0].split('_')[-1])
#Fetch the wire index within its group.

##########################################################################
def add_io_pins(G):
    """Adds the IO pin nodes to the routing-resource graph.
        
    Parameters
    ----------
    G : nx.MultiDiGraph
        The routing-resource graph.

    Returns
    -------
    None
    """

    p = -1
    for i in range(0, IO_CAPACITY):
        p += 1
        u = "io_%d_opad_in" % i
        G.add_node(u, node_type = "io_opad_in", p = p)
        p += 1
        u = "io_%d_ipad_out" % i
        G.add_node(u, node_type = "io_ipad_out", p = p)
        p += 1
        u = "io_%d_clk_in" % i
        G.add_node(u, node_type = "io_clk_in", p = p)
##########################################################################

##########################################################################
def add_cluster_pins(G):
    """Adds the cluster pin nodes to the routing-resource graph.
        
    Parameters
    ----------
    G : nx.MultiDiGraph
        The routing-resource graph.

    Returns
    -------
    None
    """

    p = -1
    for i in range(0, cluster_inputs):
        p += 1
        ble_cnt = i / (cluster_inputs / N)
        ble_i_cnt = i % (cluster_inputs / N)
        u = "ble_%d_cb_out_%d" % (ble_cnt, ble_i_cnt)
        G.add_node(u, node_type = "cb_out", p = p)

    for n in range(0, N):
        for o in range(0, O):
            p += 1
            u = "ble_%d_o_%d" % (n, o)
            G.add_node(u, node_type = "clb_out", p = p)
    p += 1
    G.add_node("ble_clk",  node_type = "clb_clk", p = p)
##########################################################################

##########################################################################
def compose_channels(G, H, V):
    """Composes the vertical and the horizontal channels
    from the high-level specification.

    Parameters
    ----------
    G : nx.MultiDiGraph
        The routing-resource graph.
    H : List[Tuple[int]]
        The list of track lengths and occurrences per BLE of each,
        in the horizontal channel.
    V : List[Tuple[int]]
        The list of track lengths and occurrences per BLE of each,
        in the vertical channel.

    Returns
    -------
    None

    Notes
    -----
    For now we assume that all wires are unidirectional. Hence the occurence
    of each track is effectively multiplied by two.
    """

    p_h = -1
    p_v = -1
    for n in range(0, N):
        for h in H:
            for i in range(0, h[1]):
                u = "ble_%d_H%d_L_%d" % (n, h[0], i)
                p_h += 1
                G.add_node(u, node_type = "h_track", p = p_h)
                u = "ble_%d_H%d_R_%d" % (n, h[0], i)
                p_h += 1
                G.add_node(u, node_type = "h_track", p = p_h)
        for v in V:
            for i in range(0, v[1]):
                for tap in range(0, min(v[0], tap_M) if SEPARATE_TAPS else 1):
                    up_template = "ble_%d_V%d_U_%d_tap_%d"
                    u = up_template % (n, v[0], i, tap)
                    p_v += 1
                    G.add_node(u, node_type = "v_track", p = p_v)
                    if tap > 0:
                        mux_type = "V%d_tap_%d" % (v[0], tap)
                        offset = (0, v[0] if not SEPARATE_TAPS else\
                                  tap_phi if tap > 1 else v[0] - (tap_M - 1))
                        G.add_edge(up_template % (n, v[0], i, tap - 1), u,\
                                   mux_type = mux_type, offset = offset)
                    down_template = "ble_%d_V%d_D_%d_tap_%d"
                    u = down_template % (n, v[0], i, tap)
                    p_v += 1
                    G.add_node(u, node_type = "v_track", p = p_v)
                    if tap > 0:
                        mux_type = "V%d_tap_%d" % (v[0], tap)
                        offset = (0, -1 * v[0] if not SEPARATE_TAPS else\
                                  -1 * (tap_phi if tap > 1 else v[0] - (tap_M - 1)))
                        G.add_edge(down_template % (n, v[0], i, tap - 1), u,\
                                   mux_type = mux_type, offset = offset)
##########################################################################

##########################################################################
def add_taps(G):
    """Adds taps into the >>tap_M<< connection blocks from the end of each wire,
    spaced at >>tamp_phi<<, for vertical wires, and only the last block for the
    horizontal ones (due to the vertically stacked BLE layout assumption).

    Parameters
    ----------
    G : nx.MultiDiGraph
        The routing-resource graph.

    Returns
    -------
    None
    """

    cb_iter = 0
    for u, attrs in G.nodes(data = True):
        if attrs["node_type"] == "h_track":
            target_ble = u.split('H')[0]
            L = int(u.split('_')[2][1:])
            d = u.split('_')[-2]
            offset = (-1 * L if d == 'L' else L, 0)
            type_count = [h for h in H if h[0] == L][0][1]
            con_freq = 1
            if DISJOINT_CB:
                for i in range(0, cluster_inputs / N):
                    if type_count < cluster_inputs / N:
                        sought_index = get_index(u)
                        if sought_index == type_count:
                            sought_index = 0
                        if i % type_count != sought_index:
                            continue
                    else:
                        sought_index = i
                        if get_index(u) % (cluster_inputs / N) != sought_index:
                            continue
                    G.add_edge(u, target_ble + "cb_out_%d" % i, offset = offset, mux_type = "cb")
        elif attrs["node_type"] == "v_track":
            tap = int(u.split('_')[-1])
            if SEPARATE_TAPS:
                tap_list = [u]
            else:
                tap_list = []
                for sweep_tap in range(0, tap_M):
                    tap_list.append(u.replace("_tap_%d" % tap, "_tap_%d" % sweep_tap))
            for sweep_u in tap_list:
                tap = int(sweep_u.split('_')[-1])
                target_ble = sweep_u.split('V')[0]
                L = int(sweep_u.split('_')[2][1:])
                type_count = [v for v in V if v[0] == L][0][1]
                d = sweep_u.split("_tap")[0].split('_')[-2]
                if SEPARATE_TAPS:
                    y_offset = tap_phi if tap > 0 else L - (tap_M * tap_phi - 1)
                else:
                    y_offset = L - (tap_M * tap_phi - 1) + tap
                offset = (0, -1 * y_offset if d == 'D' else y_offset)
                attach_u = sweep_u if SEPARATE_TAPS else sweep_u.replace("_tap_%d" % tap, "_tap_0")
                if DISJOINT_CB:
                    for i in range(0, cluster_inputs / N):
                        if type_count < cluster_inputs / N:
                            sought_index = get_index(u)
                            if sought_index == type_count:
                                sought_index = 0
                            if i % type_count != sought_index:
                                continue
                        else:
                            sought_index = i
                            if get_index(u) % (cluster_inputs / N) != sought_index:
                                continue
                        G.add_edge(u, target_ble + "cb_out_%d" % i, offset = offset, mux_type = "cb", tap = tap)
##########################################################################

##########################################################################
def add_clb_to_sb(G):
    """Adds CLB output drivers to the SB muxes. We asume simply that each BLE
    output drives those and only those SB muxes that are at its height. This
    is what likely happens in Agilex, 7-Series (See Fig. 19 of Morten's report),
    and what makes a lot of sense given the pitch and resistance issues, as otherwise
    we would need even more tracks to bring the output of a BLE to the appropriate
    SB mux.

    Parameters
    ----------
    G : nx.MultiDiGraph
        The routing-resource graph.

    Returns
    -------
    None
    """

    for u, attrs in G.nodes(data = True):
        if attrs["node_type"] == "h_track":
            mux_type = u.split('_')[2]
            source_ble = u.split('H')[0]
            for o in range(0, O):
                G.add_edge(source_ble + "o_%d" % o, u, mux_type = mux_type)
        elif attrs["node_type"] == "v_track":
            if "_tap" in u:
                tap = int(u.split('_')[-1])
                if tap != 0 and tap_M > 1:
                    continue
            mux_type = u.split('_')[2] + "_tap_0"
            source_ble = u.split('V')[0]
            for o in range(0, O):
                G.add_edge(source_ble + "o_%d" % o, u, mux_type = mux_type)
##########################################################################

##########################################################################
def read_sb_pattern_from_file():
    """Reads an SB-switch pattern from a file, generated externally.
    Adopted and banned edges are both read. Any other potential edge will
    be inserted later depending on whether a search is being performed or not.

    Parameters
    ----------
    None

    Returns
    -------
    List[Tuple[str]]
        A list of static-graph edges, specifying the switch-block switches.
    List[Tuple[str]]
        A list of static-graph edges, specifying the switch-block switches
        no longer in consideration.
    """

    if args.load_pattern_from_file is None:
        return None, None

    with open(args.load_pattern_from_file, "r") as inf:
        lines = inf.readlines()

    used = []
    banned = []
    for line in lines:
        if line.startswith('~'):
            banned.append(line.strip())
        else:
            used.append(line.strip())

    return used, banned
##########################################################################
  
##########################################################################
def read_sb_pattern_from_vpr(vpr_log_filename, rr_graph_filename):
    """Reads a switch-block pattern by extracting the used potential edges
    from a VPR routing log.

    Parameters
    ----------
    vpr_log_filename : str
        Name of the VPR routing log file.
    rr_graph_filename : str
        Name of the rr-graph file used to produce the routing log. 

    Returns
    -------
    List[Tuple[str]]
        A list of static-graph edges, specifying the switch-block switches.
    List[Tuple[str]]
        A list of static-graph edges, specifying the switch-block switches
        no longer in consideration.
    """

    #------------------------------------------------------------------------#
    def get_costs(lines):
        """Reads the final costs of the edge-splitting nodes from the VPR log,
        if it contains them.

        Parameters
        ----------
        lines : List[str]
            List of lines read from the VPR log.

        Returns
        -------
        Dict[str, float]
            Dictionary of costs, indexed by edge-splitter names.
        """

        cost_dict = {}
        for line in lines:
            if " -> " in line:
                e = line.split('(', 1)[1].split(')', 1)[0]
                c = float(line.split(" -> ")[1].split('(')[1].split(')')[0])
                #FIXME: In the final printout, the stored and the computed base costs differ.
                #If printed during iterations, all seems fine.
                cost_dict.update({e : c})
        
        return cost_dict
    #------------------------------------------------------------------------#

    #Retrieve the stored accepted/banned edges.
    store_file = "stored_edges.save"
    stored_edges = []
    stored_banned_edges = []
    try:
        with open(store_file, "r") as inf:
            for line in inf.readlines():
                if line[0] != "~":
                    stored_edges.append(line.strip())
                else:
                    stored_banned_edges.append(line.strip()[1:])
    except:
        pass

    #Now read the current usage statistics.
    lines = []
    try:
        with open(vpr_log_filename, "r") as inf:
            lines = inf.readlines()
    except:
        pass

    costs = get_costs(lines)
    free = [e for e in costs if costs[e] == 0 and not e in stored_edges]

    if GREEDY_SWITCH_SEARCH:
        free = []
    
    rd = 0
    usages = {}
    banned_edges = []
    for line in lines:
        if "Total routing area:" in line:
            rd = 1
            continue
        if rd == 1 and "Edge-splitter costs:" in line:
            rd += 1
            continue
        if rd < 2:
            continue
    
        try:
            e = line.split('(', 1)[1].split(')', 1)[0]
            if e in stored_edges:
                continue
            usage = int(line.split(',', 1)[0].rsplit('(', 1)[1])
            if usage == 0:
                banned_edges.append(e)
            usages.update({e : usage})
        except:
            break

    try:
        max_usage = max(usages.values())
    except:
        max_usage = 0.0
        global FINALIZE_CLIQUE
        if NO_STORED_UNUSED_EDGE_INCREASE:
            with open("base_costs.dump", "r") as inf:
                lines = inf.readlines()
    
                if lines[2][0] == '0':
                    FINALIZE_CLIQUE = True
        else:
            FINALIZE_CLIQUE = True
        print "V1: No new switches added to the pattern. Finalizing the clique."
        #raw_input()

    usage_thr = float(max_usage) / ADOPTION_THRESHOLD

    if max_usage <= 1e-6:
        used_potential_edges = []
        global FINALIZE_CLIQUE
        if NO_STORED_UNUSED_EDGE_INCREASE:
            with open("base_costs.dump", "r") as inf:
                lines = inf.readlines()
    
                if lines[2][0] == '0':
                    FINALIZE_CLIQUE = True
        else:
            FINALIZE_CLIQUE = True
        print "V2: No new switches added to the pattern. Finalizing the clique."
        #raw_input()
    elif not NO_STORED_UNUSED_EDGE_INCREASE:
        used_potential_edges = [e for e in usages if usages[e] >= usage_thr]
        if PICK_TOP_ONLY:
            used_potential_edges = [e for e in usages if usages[e] == max_usage]
            if len(used_potential_edges) > 1:
                td_dict = {}
                delay_file = "arc_delays.log"
                with open(delay_file, "r") as inf:
                    lines = inf.readlines()
                for line in lines:
                    td_dict.update({line.split()[0] : float(line.split()[1])})
                used_potential_edges = [min(used_potential_edges, key = lambda s : td_dict[s])]
 
    #TODO: Choose this by a variable.
    if free:
        used_potential_edges = free

    if NO_STORED_UNUSED_EDGE_INCREASE:
        #NOTE: This means that we are in the wafer-creation segment, and that all the accepted edges have
        #already been stored. Hence, we should not add any new ones.
        used_potential_edges = []

    print("Used edges: %d\n" % len(used_potential_edges))

    ban_thr = 3500
    ban_counters = {e.split()[0] : int(e.split()[1]) for e in stored_banned_edges\
                    if (e.split()[0] in banned_edges) or (int(e.split()[1]) >= ban_thr)}
    for e in banned_edges:
        if e in ban_counters:
            if not NO_STORED_UNUSED_EDGE_INCREASE:
                ban_counters[e] += 1
        else:
            ban_counters.update({e : 1})

    used_potential_edges = set(used_potential_edges + stored_edges)
    #banned_edges = set([names[seg_id] for seg_id in banned_edges] + stored_banned_edges)
    banned_edges = set([e for e in ban_counters if ban_counters[e] >= ban_thr])

    txt = ""
    for e in sorted(used_potential_edges):
        txt += "%s\n" % e
    for e in sorted(ban_counters):
        txt += "~%s %d\n" % (e, ban_counters[e])
    #for e in sorted(banned_edges):
    #    txt += "~%s\n" % e


    with open(store_file, "w") as outf:
        outf.write(txt[:-1])

    #with open("ban.log", "a") as outf:
    #    outf.write("%s\n\n" % str(banned_edges))

    return used_potential_edges, banned_edges
##########################################################################

##########################################################################
def lut_canonical_potential_edge(potential_edge):
    """Returns a canonical name of a potential edge, with respect to LUT height.

    Parameters
    ----------
    potential_edge : str
        Instantiated name of the potential edge to be canonicized.

    Returns
    -------
    str
        A canonical potential edge.
    """

    prefix, u, v = potential_edge.split("__")

    lut_span = int(v.split('_')[1]) - int(u.split('_')[1])
    if lut_span < 0:
        offset_str = "lutm%d_" % abs(lut_span)
    else:
        offset_str = "lutp%d_" % abs(lut_span)

    canonical = "__".join([prefix, '_'.join(u.split('_')[2:]),\
                offset_str + '_'.join(v.split('_')[2:])])

    return canonical    
##########################################################################

##########################################################################
def add_sb_to_sb(G, make_clique = False):
    """Assigns track drivers to the SB muxes. Each mux gets driven from
    all the tracks entering the same BLE section, apart from those
    coming from the direction of the wire it drives. This creates a
    reasonably flexible and very straight-forward to implement switch,
    akin to subset, with full turn and span flexibility.

    Parameters
    ----------
    G : nx.MultiDiGraph
        The routing-resource graph.
    make_clique : Optional[bool], default = False
        Local equivalent to args.make_sb_clique.
        It does not make sense to construct the clique
        while padding is in progress.    


    Returns
    -------
    None

    Notes
    -----
    We can optionally balance the horizontal and vertical channel mux widths
    and maybe even snap them to the nearest larger SRAM-optimal size. Then
    the extra wires can be used to expand connectivity from the subsets, which
    is seemingly what Agilex does.

    VPR uses the following coordinate layout scheme:

    ------           -----          
    ..... |          ..... |
    .   . |          .   . |
    .   . |          .   . |
    .   . |          .   . |
    .   . |          .   . |
    ..... |          ..... |


    ------           -----          
    ..... |          ..... |
    .   . |          .   . |
    .   . |          .   . |
    .   . |          .   . |
    .   . |          .   . |
    ..... |          ..... |

    
    Hence U-L turns go with (0, L - 1) offset,
          U-R turns go with (1, L - 1) offset,
          D-L turns go with (0, -L) offset,
          D-R turns go with (1, -L) offset,
          R-U turns go with (L - 1, 1) offset,
          R-D turns go with (L - 1, 0) offset,
          L-U turns go with (-L, 1) offset,
          L-D turns go with (-L, 0) offset,
          U-U straights go with (0, L) offset,
          D-D straights go with (0, -L) offset,
          R-R straights go with (L, 0) offset,
          L-L straights go with (-L, 0) offset.    
    """

    offset_dict = {\
                   ('U', 'L') : lambda L : (0, L - 1),\
                   ('U', 'R') : lambda L : (1, L - 1),\
                   ('D', 'L') : lambda L : (0, -1 * L),\
                   ('D', 'R') : lambda L : (1, -1 * L),\
                   ('R', 'U') : lambda L : (L - 1, 1),\
                   ('R', 'D') : lambda L : (L - 1, 0),\
                   ('L', 'U') : lambda L : (-1 * L, 1),\
                   ('L', 'D') : lambda L : (-1 * L, 0),\
                   ('U', 'U') : lambda L : (0, L),\
                   ('D', 'D') : lambda L : (0, -1 * L),\
                   ('R', 'R') : lambda L : (L, 0),\
                   ('L', 'L') : lambda L : (-1 * L, 0)\
                  }

    wire_dict = {}
    for u, attrs in G.nodes(data = True):
        if attrs["node_type"] in ("h_track", "v_track"):
            ble = '_'.join(u.split('_')[:2])
            d = u.split('_')[-2]
            if "_tap" in u:
                tap = int(u.split('_')[-1])
                if tap != 0 and tap != tap_M - 1:
                    continue
                d = u.split("_tap")[0].split('_')[-2]
            try:
                wire_dict[ble][d].append(u)
            except:
                try:
                    wire_dict[ble].update({d : [u]})
                except:
                    wire_dict.update({ble : {d : [u]}})

    is_loopback = lambda d_in, d_out : True if set([d_in, d_out])\
                  in (set(['L', 'R']), set(['U', 'D'])) else False


    if args.load_pattern_from_file is not None:
            existing_edges, banned_edges = read_sb_pattern_from_file()
    elif make_clique:
        if FIRST_CLIQUE_ITER:
            existing_edges, banned_edges = [], []
        else:
            existing_edges, banned_edges = read_sb_pattern_from_vpr("vpr_stdout.log",\
                                                                    args.arc_name.rsplit(".xml", 1)[0] + "_rr.xml")
    if make_clique or args.load_pattern_from_file is not None:
        added = 0
        nodes = list(sorted(G.nodes()))
        for wire_target in nodes:
            if not G.node[wire_target]["node_type"] in ("h_track", "v_track"):
                continue
            d_target = wire_target.split('_')[-2]
            if "_tap" in wire_target:
                tap = int(wire_target.split('_')[-1])
                if tap != 0 and tap_M > 1:
                    continue
                d_target = wire_target.split("_tap")[0].split('_')[-2]
            mux_type = wire_target.split('_')[2] + ("_tap_0" if "tap" in wire_target else '')
            target_ble = int(wire_target.split('_')[1])
            for wire_source in nodes:
                if not G.node[wire_source]["node_type"] in ("h_track", "v_track"):
                    continue
                L = int(wire_source.split('_')[2][1:])
                d_source = wire_source.split('_')[-2]
                if "_tap" in wire_source:
                    tap = int(wire_source.split('_')[-1])
                    if tap == 0 and tap_M > 1 and SEPARATE_TAPS:
                        continue
                    if tap_M > 1 and SEPARATE_TAPS:
                        L = tap_phi
                    d_source = wire_source.split("_tap")[0].split('_')[-2]
                if is_loopback(d_source, d_target):
                    #FIXME: Actually, we want to be general and consider U-turns. Allow this!
                    continue
                source_ble = int(wire_source.split('_')[1])
                if abs(source_ble - target_ble) > MAX_BLE_SPAN:
                    continue
                offset = offset_dict[(d_source, d_target)](L)
                cost = -1
                break_node = "potential_edge__%s__%s" % (wire_source, wire_target)
                if lut_canonical_potential_edge(break_node) in banned_edges:
                    #Edges that were previously not used at all can be skipped in the next iteration,
                    #without much loss, in general.
                    continue
                if lut_canonical_potential_edge(break_node) in existing_edges:
                    #print "Adding (%s -> %s)" % (wire_source, wire_target)
                    cost = 0
                    G.add_edge(wire_source, wire_target, mux_type = mux_type, offset = offset, tap = -1)
                elif FINALIZE_CLIQUE or (args.load_pattern_from_file is not None and not PATTERN_IS_INIT_ONLY):
                    continue
                else:
                    G.add_node(break_node, node_type = "potential_edge_" + G.node[wire_target]["node_type"],\
                               cost = cost, seg = lut_canonical_potential_edge(break_node))
                    G.add_edge(wire_source, break_node, mux_type = lut_canonical_potential_edge(break_node),\
                               offset = offset, tap = -1)
                    #NOTE: Maybe it would be advisable to switch the order of dummy-real, to make the lookahead more effective.
                    G.add_edge(break_node, wire_target, mux_type = mux_type, offset = (0, 0), tap = -1)
                    added += 1
        
        print("Potential edges: %d" % added)

        #Now add the dummy switches to hanging wires in the fixed pattern, to be able to measure the delay.
        for wire_target in nodes:
            if not G.node[wire_target]["node_type"] in ("h_track", "v_track"):
                continue
            real_fanout = [c for c in G[wire_target] if not "potential_edge" in c and G.node[c]["node_type"] in ("h_track", "v_track")]
            if real_fanout:
                continue
            d_target = wire_target.split('_')[-2]
            if "_tap" in wire_target:
                tap = int(wire_target.split('_')[-1])
                if tap != 0 and tap_M > 1:
                    continue
                d_target = wire_target.split("_tap")[0].split('_')[-2]
            mux_type = wire_target.split('_')[2] + ("_tap_0" if "tap" in wire_target else '')
            target_ble = int(wire_target.split('_')[1])
            for wire_source in nodes:
                if wire_source != wire_target:
                    continue
                L = int(wire_source.split('_')[2][1:])
                d_source = wire_source.split('_')[-2]
                if "_tap" in wire_source:
                    tap = int(wire_source.split('_')[-1])
                    if tap == 0 and tap_M > 1 and SEPARATE_TAPS:
                        continue
                    if tap_M > 1 and SEPARATE_TAPS:
                        L = tap_phi
                    d_source = wire_source.split("_tap")[0].split('_')[-2]
                if is_loopback(d_source, d_target):
                    #FIXME: Actually, we want to be general and consider U-turns. Allow this!
                    continue
                source_ble = int(wire_source.split('_')[1])
                if abs(source_ble - target_ble) > MAX_BLE_SPAN:
                    continue
                offset = offset_dict[(d_source, d_target)](L)
                G.add_edge(wire_source, wire_target, mux_type = mux_type, offset = offset, tap = -1,\
                           delay_operating_point_only = True)

        return

    for ble in wire_dict:
        for d_target in wire_dict[ble]:
            for d_source in wire_dict[ble]:
                if is_loopback(d_target, d_source):
                    continue
                for wire_target in wire_dict[ble][d_target]:
                    tap = int(wire_target.split('_')[-1])
                    if "_tap" in wire_target and tap != 0 and tap_M > 1:
                        continue
                    mux_type = wire_target.split('_')[2] + ("_tap_0" if "tap" in wire_target else '')
                    if wire_target.split('_')[2][0] == 'H':
                        target_type_count = [h for h in H if h[0] == int(wire_target.split('_')[2][1:])][0][1]
                    else:
                        target_type_count = [v for v in V if v[0] == int(wire_target.split('_')[2][1:])][0][1]
                    for wire_source in wire_dict[ble][d_source]:
                        if wire_source.split('_')[2][0] == 'H':
                            source_type_count = [h for h in H if h[0] == int(wire_source.split('_')[2][1:])][0][1]
                        else:
                            source_type_count = [v for v in V if v[0] == int(wire_source.split('_')[2][1:])][0][1]
                        L = int(wire_source.split('_')[2][1:])
                        if "_tap" in wire_source:
                            tap = int(wire_source.split('_')[-1])
                            if tap == 0 and tap_M > 1 and SEPARATE_TAPS:
                                continue
                            if tap_M > 1 and SEPARATE_TAPS:
                                L = tap_phi
                        offset = offset_dict[(d_source, d_target)](L)
                        if DISJOINT_SB:
                            if source_type_count < target_type_count:
                                sought_index = get_index(wire_source) + 1
                                if sought_index == source_type_count:
                                    sought_index = 0
                                if get_index(wire_target) % source_type_count != sought_index:
                                    continue
                            else:
                                sought_index = get_index(wire_target) + 1
                                if sought_index == target_type_count:
                                    sought_index = 0
                                if get_index(wire_source) % target_type_count != sought_index:
                                    continue
                            #sought_index = get_index(wire_target) + 1
                            #if sought_index == target_type_count:
                            #    sought_index = 0
                            #if get_index(wire_source) != sought_index:
                            #    continue
                        G.add_edge(wire_source, wire_target, mux_type = mux_type, offset = offset, tap = -1)
                        #if "V4" in wire_source and "V4" in wire_target:
                        #    print wire_source, wire_target, L, offset

    #-------------------------------------------------------------------------#
    def is_twist_wire(wire):
        """Checks whetehr the wire is a twist wire or not.

        Parameters
        ----------
        wire : str
            Wire being checked.

        Returns
        --------
        bool
            True if yes, else False.
        """

        v_twist_len = max(1, K6N8_LUT4 / KN_LUT4)
        h_twist_len = 1
        
        L = int(wire.split('_')[2][1:])
        if 'V' in wire and L != v_twist_len:
            return False

        if 'H' in wire and L != h_twist_len:
            return False
   
        return True
    #-------------------------------------------------------------------------#

    if ADD_LEN_1_TWISTS:
        for ble in wire_dict:
            ble_ind = int(ble.split('_')[-1])
            up_offset = 0
            up_ind = ble_ind + 1
            if up_ind == N:
                up_ind = 0
                up_offset = 1
            down_offset = 0
            down_ind = ble_ind - 1
            if down_ind < 0:
                down_ind = N - 1
                down_offset = -1
            for d_target in wire_dict[ble]:
                for d_source in wire_dict[ble.replace(str(ble_ind), str(up_ind))]:
                    if is_loopback(d_target, d_source):
                        continue
                    if d_source != d_target and ONLY_CONTINUATION_TWISTS:
                        continue
                    for wire_target in wire_dict[ble][d_target]:
                        if not is_twist_wire(wire_target):
                            continue
                        tap = int(wire_target.split('_')[-1])
                        if "_tap" in wire_target and tap != 0 and tap_M > 1:
                            continue
                        mux_type = wire_target.split('_')[2] + ("_tap_0" if "tap" in wire_target else '')
                        if wire_target.split('_')[2][0] == 'H':
                            target_type_count = [h for h in H if h[0] == int(wire_target.split('_')[2][1:])][0][1]
                        else:
                            target_type_count = [v for v in V if v[0] == int(wire_target.split('_')[2][1:])][0][1]
                        for wire_source in wire_dict[ble.replace(str(ble_ind), str(up_ind))][d_source]:
                            L = int(wire_source.split('_')[2][1:])
                            if not is_twist_wire(wire_source):
                                continue
                            if "_tap" in wire_source:
                                tap = int(wire_source.split('_')[-1])
                                if tap == 0 and tap_M > 1 and SEPARATE_TAPS:
                                    continue
                                if tap_M > 1 and SEPARATE_TAPS:
                                    L = tap_phi
                            offset = offset_dict[(d_source, d_target)](L)
                            if CUT_CROSS_CLB_TWISTS and up_offset:
                                continue
                            offset = (offset[0], offset[1] + up_offset)
                            sought_index = get_index(wire_target) + 1
                            if sought_index == target_type_count:
                                sought_index = 0
                            if DISJOINT_SB and sought_index != get_index(wire_source):
                                continue
                            G.add_edge(wire_source, wire_target, mux_type = mux_type, offset = offset, tap = -1)
                            #if "V4" in wire_source and "V4" in wire_target:
                            #    print "up_ble: ", wire_source, wire_target, L, offset
                for d_source in wire_dict[ble.replace(str(ble_ind), str(down_ind))]:
                    if is_loopback(d_target, d_source):
                        continue
                    if d_source != d_target and ONLY_CONTINUATION_TWISTS:
                        continue
                    for wire_target in wire_dict[ble][d_target]:
                        if not is_twist_wire(wire_target): 
                            continue
                        tap = int(wire_target.split('_')[-1])
                        if "_tap" in wire_target and tap != 0 and tap_M > 1:
                            continue
                        mux_type = wire_target.split('_')[2] + ("_tap_0" if "tap" in wire_target else '')
                        if wire_target.split('_')[2][0] == 'H':
                            target_type_count = [h for h in H if h[0] == int(wire_target.split('_')[2][1:])][0][1]
                        else:
                            target_type_count = [v for v in V if v[0] == int(wire_target.split('_')[2][1:])][0][1]
                        for wire_source in wire_dict[ble.replace(str(ble_ind), str(down_ind))][d_source]:
                            L = int(wire_source.split('_')[2][1:])
                            if not is_twist_wire(wire_source):
                                continue
                            if "_tap" in wire_source:
                                tap = int(wire_source.split('_')[-1])
                                if tap == 0 and tap_M > 1 and SEPARATE_TAPS:
                                    continue
                                if tap_M > 1 and SEPARATE_TAPS:
                                    L = tap_phi
                            offset = offset_dict[(d_source, d_target)](L)
                            if CUT_CROSS_CLB_TWISTS and down_offset:
                                continue
                            offset = (offset[0], offset[1] + down_offset)
                            sought_index = get_index(wire_target) + 1
                            if sought_index == target_type_count:
                                sought_index = 0
                            if DISJOINT_SB and sought_index != get_index(wire_source):
                                continue
                            G.add_edge(wire_source, wire_target, mux_type = mux_type, offset = offset, tap = -1)
                            #if "V4" in wire_source and "V4" in wire_target:
                            #    print "down_ble: ", wire_source, wire_target, L, offset
##########################################################################

##########################################################################
def add_cluster_sinks_and_sources(G):
    """Adds the sinks and the sources for the cluster ports, used
    to represent the logical equivalence classes.

    Parameters
    ----------
    G : nx.MultiDiGraph
        The routing-resource graph.

    Returns
    -------
    None
    """

    for u, attrs in G.nodes(data = True):
        if attrs["node_type"] == "cb_out":
            G.add_edge(u, "I_SINK", mux_type = "__vpr_delayless_switch__")
        elif attrs["node_type"] == "clb_out":
            G.add_edge("O_SOURCE", u, mux_type = "__vpr_delayless_switch__")
        elif attrs["node_type"] == "clb_clk":
            G.add_edge(u, "CLK_SINK", mux_type = "__vpr_delayless_switch__")
    
    G.node["I_SINK"]["node_type"] = G.node["O_SOURCE"]["node_type"] = G.node["CLK_SINK"]["node_type"] = "dummy"
    G.node["I_SINK"]['p'] = 0
    G.node["O_SOURCE"]['p'] = 1
    G.node["CLK_SINK"]['p'] = 2
##########################################################################

##########################################################################
def add_io_sinks_and_sources(G):
    """Adds the sinks and the sources for the IO ports.

    Parameters
    ----------
    G : nx.MultiDiGraph
        The routing-resource graph.

    Returns
    -------
    None
    """

    for u, attrs in G.nodes(data = True):
        if attrs["node_type"] == "io_opad_in":
            io_cnt = int(u.split('_')[1])
            v = "IO_%d_OPAD_SINK" % io_cnt
            G.add_edge(u, v, mux_type = "__vpr_delayless_switch__")
            G.node[v]["node_type"] = "dummy"
            G.node[v]['p'] = 0
        elif attrs["node_type"] == "io_ipad_out":
            io_cnt = int(u.split('_')[1])
            v = "IO_%d_IPAD_SOURCE" % io_cnt
            G.add_edge(v, u, mux_type = "__vpr_delayless_switch__")
            G.node[v]["node_type"] = "dummy"
            G.node[v]['p'] = 1
        elif attrs["node_type"] == "io_clk_in":
            io_cnt = int(u.split('_')[1])
            v = "IO_%d_CLK_SINK" % io_cnt
            G.add_edge(u, v, mux_type = "__vpr_delayless_switch__")
            G.node[v]["node_type"] = "dummy"
            G.node[v]['p'] = 2
##########################################################################

##########################################################################
def get_chan_width(C):
    """Returns the channel width.

    Parameters
    ----------
    C : List[Tuple[int]]
        The list of track lengths and occurrences per BLE of each,
        in the given channel (V or H).

    Returns
    -------
    int
        Channnel width.
    """

    return 2 * N * sum([t[0] * t[1] for t in C])
##########################################################################

##########################################################################
def generate_grid(grid_w, grid_h):
    """Generates the grid layout according to the specification.

    Parameters
    ----------
    grid_w : int
        Number of columns (vertical channels).
    grid_h : int
        Number of rows (horizontal channels).

    Returns
    -------
    Dict[Tuple[int], str]
        A dictionary of block types, indexed by the grid coordinates.
    """

    grid = {}
    for x in range(0, grid_w):
        nomem = False
        nomult = False
        for y in range(0, grid_h):
            grid.update({(x, y) : "clb"})
            if x < io_crop or x > grid_w - 1 - io_crop or y < io_crop or y > grid_h - 1 - io_crop:
                grid[(x, y)] = "io"
     
    if cut_corners:
        grid[(0, 0)] = "empty"
        grid[(0, grid_h - 1)] = "empty"
        grid[(grid_w - 1, 0)] = "empty"
        grid[(grid_w - 1, grid_h - 1)] = "empty"
    if TOP_BOTTOM_IO:
        for y in range(0, grid_h):
            grid[(0, y)] = "empty"
            grid[(grid_w - 1, y)] = "empty"

    return grid
##########################################################################

##########################################################################
def export_chan_tags(H, V, grid_w, grid_h, potential_edge_no = 0):
    """Exports the channel tags in the VTR8 RR-graph format.

    Parameters
    ----------    
    H : List[Tuple[int]]
        The list of track lengths and occurrences per BLE of each,
        in the horizontal channel.
    V : List[Tuple[int]]
        The list of track lengths and occurrences per BLE of each,
        in the vertical channel.
    grid_w : int
        Number of columns (vertical channels).
    grid_h : int
        Number of rows (horizontal channels).
    potential_edge_no : Optional[int], default = 0
        Number of potential edges whose dummy nodes increase the formal channel width.
        This is necessary for VPR lookahead computation not to stop too early.

    Returns
    -------
    str
        Text of the tags.
    """

    h_width = get_chan_width(H) + potential_edge_no 
    v_width = get_chan_width(V) + potential_edge_no

    txt = indent + "<channels>\n"
    txt += 2 * indent + "<channel chan_width_max=\"%d\" x_min=\"%d\" y_min=\"%d\" x_max=\"%d\" y_max=\"%d\"/>\n"\
        % (max(h_width, v_width), h_width, v_width, h_width, v_width)

    for x in range(0, grid_h):
        txt += 3 * indent + "<x_list index=\"%d\" info=\"%d\"/>\n" % (x, h_width)
    for y in range(0, grid_w):
        txt += 3 * indent + "<y_list index=\"%d\" info=\"%d\"/>\n" % (y, v_width)
    
    txt += indent + "</channels>\n"

    return txt
##########################################################################

##########################################################################
def export_switches(H, V, td_cb = default_cb_delay, td_sb = {}, potential_edges = None):
    """Exports the multiplexer data in the VTR8 RR-graph format.
    Because we do not change the channel width and calculate the area
    ourselves, transistor sizes are all set to zero. Delays are likewise
    lumped into Tdel, similar to COFFE. Because we don't have the SB-taps,
    this is even more realistic than it ever was.

    Parameters
    ----------
    H : List[Tuple[int]]
        The list of track lengths and occurrences per BLE of each,
        in the horizontal channel.
    V : List[Tuple[int]]
        The list of track lengths and occurrences per BLE of each,
        in the vertical channel.
    td_cb : float
        Delay of the connection block multiplexer (and SB-CB wires)
    td_sb : Dict[str, float]
        Delays of all the switch block multiplexers, together with the
        wires they are driving, and the typical loading by other muxes.
    potential_edges : Optional[List[str]], default = None
        List of potential edges obtained with the appropriate >>make_clique<< switches turned on.
        Optionally also a dictionary of their delays can be passed.

    Returns
    -------
    str
        Text of the tags
    Dict[str, int]
        A dictionary of switch ids
    """

    template = 2 * indent + "<switch id=\"%d\" type=\"mux\" name=\"%s\">\n"\
             + 3 * indent + "<timing R=\"0\" Cin=\"0\" Cout=\"0\" Tdel=\"%g\"/>\n"\
             + 3 * indent + "<sizing mux_trans_size=\"0\" buf_size=\"0\"/>\n"\
             + 2 * indent + "</switch>\n"

    id_dict = {}

    switch_id = 0
    switch_name = "__vpr_delayless_switch__"
    id_dict.update({switch_name : switch_id})
    switches = [template % (switch_id, switch_name, 0)]

    switch_id += 1
    switch_name = "cb"
    id_dict.update({switch_name : switch_id})
    switches.append(template % (switch_id, switch_name, td_cb))

    for h_track in H:
        switch_id += 1
        switch_name = "H%d" % h_track[0]
        td = td_sb.get(switch_name, td_sb.get(switch_name, 0.5 * default_delay * (h_track[0] + 1)))
        id_dict.update({switch_name : switch_id})
        switches.append(template % (switch_id, switch_name, td))
    for v_track in V:
        for tap in range(0, tap_M if SEPARATE_TAPS else 1):
            switch_id += 1
            switch_name = "V%d_tap_%d" % (v_track[0], tap)
            td = td_sb.get(switch_name, td_sb.get(switch_name, 0.5 * default_delay * (v_track[0] + 1)))
            id_dict.update({switch_name : switch_id})
            switches.append(template % (switch_id, switch_name, td))
 
    if potential_edges is not None:
        potential_edge_delay = 1e-12
        #NOTE:VPR lookahead breaks with zero delays.
        for u in sorted(potential_edges):
            if isinstance(potential_edges, dict):
                potential_edge_delay = potential_edges[u]
            switch_id += 1
            id_dict.update({u : switch_id})
            switches.append(template % (switch_id, u, potential_edge_delay))

    txt = indent + "<switches>\n"
    txt += ''.join(switches)
    txt += indent + "</switches>\n"
        
    return txt, id_dict
##########################################################################

##########################################################################
def export_segments(H, V, potential_edges = None):
    """Exports the segment information in the VTR8 RR-graph format.
    All RC parameters are set to 0, since the delay is already lumped into
    the driver delay.

    Parameters
    ----------
    H : List[Tuple[int]]
        The list of track lengths and occurrences per BLE of each,
        in the horizontal channel.
    V : List[Tuple[int]]
        The list of track lengths and occurrences per BLE of each,
        in the vertical channel.
    potential_edges : Optional[List[str]], default = None
        List of potential edges obtained with the appropriate >>make_clique<< switches turned on.

    Returns
    -------
    str
        Text of the tags
    Dict[str, int]
        A dictionary of segment ids
    """

    txt = indent + "<segments>\n"

    template = 2 * indent + "<segment id=\"%d\" name=\"%s\">\n"\
             + 3 * indent + "<timing R_per_meter=\"0\" C_per_meter=\"0\"/>\n"\
             + 2 * indent + "</segment>\n"
    
    id_dict = {}
    seg_id = -1 
    segs = []

    for h_track in H:
        seg_id += 1
        seg_name = "H%d" % h_track[0]
        id_dict.update({seg_name : seg_id})
        segs.append(template % (seg_id, seg_name))
    for v_track in V:
        for tap in range(0, tap_M if SEPARATE_TAPS else 1):
            seg_id += 1
            seg_name = "V%d_tap_%d" % (v_track[0], tap)
            id_dict.update({seg_name : seg_id})
            segs.append(template % (seg_id, seg_name))

    if potential_edges is not None:
        for u in sorted(potential_edges):
            seg_id += 1
            id_dict.update({u : seg_id})
            segs.append(template % (seg_id, u))

    txt += ''.join(segs)
    txt += indent + "</segments>\n"
   
    return txt, id_dict
##########################################################################

##########################################################################
def export_clb_block(block_id, export_empty = False):
    """Exports the CLB block type in the VTR8 RR-graph format.
    We assume that the CLB width and height still determine the grid feature size
    and that they are logical equal (unit).

    Parameters
    ----------
    block_id : int
        The id code of the CLB block type.
    export_empty : Optional[bool], default = False
        Specifies that this is an empty clb, part of the empty ring.
    
    Returns
    -------
    str
        Text of the tags.
    """

    txt = indent + "<block_type id=\"%d\" name=\"%sclb\" width=\"1\" height=\"1\">\n" % (block_id, "empty_" if export_empty else "")
    txt += 2 * indent + "<pin_class type=\"INPUT\">\n"
    ptc = -1
    for i in range(0, cluster_inputs):
        ptc += 1
        txt += 3 * indent + "<pin ptc=\"%d\">%sclb.I[%d]</pin>\n" % (ptc, "empty_" if export_empty else "", i)
    txt += 2 * indent + "</pin_class>\n" + 2 * indent + "<pin_class type=\"OUTPUT\">\n"
    for o in range(0, N * O):
        ptc += 1
        txt += 3 * indent + "<pin ptc=\"%d\">%sclb.O[%d]</pin>\n" % (ptc, "empty_" if export_empty else "", o)
    txt += 2 * indent + "</pin_class>\n" + 2 * indent + "<pin_class type=\"INPUT\">\n"
    ptc += 1
    txt += 3 * indent + ("<pin ptc=\"%d\">%sclb.clk[0]</pin>\n" % (ptc, "empty_" if export_empty else ""))\
         + 2 * indent + "</pin_class>\n"
    txt += indent + "</block_type>\n"

    return txt
##########################################################################

##########################################################################
def export_blocks():
    """Exports all block types in the VTR8 RR-graph format.
    All but CLB are taken from the VTR's own output.

    Parameters
    ----------
    None

    Returns
    -------
    str
        Text of the tags.
    """

    empty = 2 * indent + "<block_type id=\"%d\" name=\"EMPTY\" width=\"1\" height=\"1\"/>\n"\
                       % block_ids["empty"]

    io = 2 * indent + "<block_type id=\"%d\" name=\"io\" width=\"1\" height=\"1\">\n"\
                     % block_ids["io"]
    ptc = 0
    for i in range(0, IO_CAPACITY):
        io += 3 * indent + "<pin_class type=\"INPUT\">\n"
        io += 4 * indent + "<pin ptc=\"%d\">io[%d].outpad[0]</pin>\n"\
                         % (ptc, i)
        ptc += 1
        io += 3 * indent + "</pin_class>\n"
        io += 3 * indent + "<pin_class type=\"OUTPUT\">\n"
        io += 4 * indent + "<pin ptc=\"%d\">io[%d].inpad[0]</pin>\n"\
                         % (ptc, i)
        ptc += 1
        io += 3 * indent + "</pin_class>\n"
        io += 3 * indent + "<pin_class type=\"INPUT\">\n"
        io += 4 * indent + "<pin ptc=\"%d\">io[%d].clock[0]</pin>\n"\
                         % (ptc, i)
        ptc += 1
        io += 3 * indent + "</pin_class>\n"
    io += 2 * indent + "</block_type>\n"
   
    txt = indent + "<block_types>\n"
    txt += empty
    txt += io
    txt += export_clb_block(block_ids["clb"])

    if INSERT_EMPTY_RING:
        txt += export_clb_block(block_ids["empty_clb"], export_empty = True)

    if mult_col["freq"] or mem_col["freq"]:
        with open("block_export.xml", "r") as inf:
            txt += inf.read()

    txt += indent + "</block_types>\n"

    return txt
##########################################################################

##########################################################################
def export_grid(grid):
    """Exports the FPGA grid in VTR8 RR-graph format.

    Parameters
    ----------
    grid : Dict[Tuple[int], str]
        A dictionary of block types, indexed by the grid coordinates.

    Returns
    -------
    str
        Text of the tags.
    """

    txt = indent + "<grid>\n"

    empty_clb_locs = []
    if INSERT_EMPTY_RING:
        max_h = max([h[0] for h in H])
        max_v = max([v[0] for v in V])
        ring_xs = [x for x in range(1, max_h + 1)] + [x for x in range(grid_w - 1 - max_h, grid_w - 1)]
        ring_ys = [y for y in range(1, max_v + 1)] + [y for y in range(grid_h - 1 - max_v, grid_h - 1)]
        for x in ring_xs:
            for y in range(1, grid_h - 1):
                empty_clb_locs.append((x, y))
        for y in ring_ys:
            for x in range(1, grid_w -1):
                if x in ring_xs:
                    continue
                empty_clb_locs.append((x, y))

    for coords in sorted(grid):
        width_offset = 0
        height_offset = 0
        x, y = coords
        block_type = grid[coords]
        if (x, y) in empty_clb_locs:
            block_type = "empty_clb"
        if block_type == "mult":
            height_offset = (y - mult_col["start"][1]) % mult_col["height"]
        elif block_type == "mem":
            height_offset = (y - mem_col["start"][1]) % mem_col["height"]
        txt += 2 * indent + "<grid_loc x=\"%d\" y=\"%d\" block_type_id=\"%d\" width_offset=\"%d\" height_offset=\"%d\"/>\n"\
                          % (x, y, block_ids[block_type], width_offset, height_offset)
    txt += indent + "</grid>\n"
       
    return txt
##########################################################################

##########################################################################
def export_rr_nodes(G, grid, init = 0):
    """Exports all nodes of the RR-graph, in the VTR8 RR-graph format.

    Parameters
    ----------
    G : nx.MultiDiGraph
        The routing-resource graph.
    grid : Dict[Tuple[int], str]
        A dictionary of block types, indexed by the grid coordinates.
    init : Optional[int], default = 0
        The initial value of the ID counter.

    Returns
    -------
    str
        Text of the tags.
    Dict[str, Dict[Tuple[int], int]]
        Mapping between the RR-graph nodes in the static form (G)
        and the tile coordinates and node ids.
    """

    node_id = init

    template = 2 * indent + "<node id=\"%d\" type=\"%s\" capacity=\"%d\">\n"\
             + 3 * indent + "<loc xlow=\"%d\" ylow=\"%d\" xhigh=\"%d\" yhigh=\"%d\" %sptc=\"%d\"/>\n"\
             + 3 * indent + "<timing R=\"0\" C=\"0\"/>\n%s"\
             + 2 * indent + "</node>\n"

    txt = indent + "<rr_nodes>\n"

    u_counts = {G.node[u]['p'] : u for u in G if G.node[u]["node_type"] in ("cb_out", "clb_out", "clb_clk")}
    export_u_counts = {}

    io_fanin_dict = {}
    io_fanout_dict = {}
    #IO pins
    for coords in sorted(grid):
        if grid[coords] != "io":
            continue
        x, y = coords
        ptc = 0
        for i in range(0, IO_CAPACITY):
            txt += template % (node_id, "SINK", 1, x, y, x, y, '', ptc, '')
            u = "IO_%d_OPAD_SINK" % i
            try:
                export_u_counts[u].update({coords : node_id})
            except:
                export_u_counts.update({u : {coords : node_id}})
            node_id += 1
            ptc += 1
            txt += template % (node_id, "SOURCE", 1, x, y, x, y, '', ptc, '')
            u = "IO_%d_IPAD_SOURCE" % i
            try:
                export_u_counts[u].update({coords : node_id})
            except:
                export_u_counts.update({u : {coords : node_id}})
            node_id += 1
            ptc += 1
            txt += template % (node_id, "SINK", 1, x, y, x, y, '', ptc, '')
            u = "IO_%d_CLK_SINK" % i
            try:
                export_u_counts[u].update({coords : node_id})
            except:
                export_u_counts.update({u : {coords : node_id}})
            node_id += 1
            ptc += 1
        ptc = 0
        for i in range(0, IO_CAPACITY):
            txt += template % (node_id, "IPIN", 1, x, y, x, y, "side=\"LEFT\" ", ptc, '')
            u = "io_%d_opad_in" % i
            try:
                export_u_counts[u].update({coords : node_id})
            except:
                export_u_counts.update({u : {coords : node_id}})
            node_id += 1
            ptc += 1
            txt += template % (node_id, "OPIN", 1, x, y, x, y, "side=\"LEFT\" ", ptc, '')
            u = "io_%d_ipad_out" % i
            try:
                export_u_counts[u].update({coords : node_id})
            except:
                export_u_counts.update({u : {coords : node_id}})
            node_id += 1
            ptc += 1
            txt += template % (node_id, "IPIN", 1, x, y, x, y, "side=\"LEFT\" ", ptc, '')
            u = "io_%d_clk_in" % i
            try:
                export_u_counts[u].update({coords : node_id})
            except:
                export_u_counts.update({u : {coords : node_id}})
            node_id += 1
            ptc += 1

    #CLB pins
    for coords in sorted(grid):
        if grid[coords] != "clb":
            continue
        x, y = coords

        #Export the cluster_inputs and clk sinks and the O source.
        #TODO: Once we switch to multiple equivalence classes, we will need to create multiple nodes.
        
        txt += template % (node_id, "SINK", cluster_inputs, x, y, x, y, '', 0, '')
        u = "I_SINK"
        try:
            export_u_counts[u].update({coords : node_id})
        except:
            export_u_counts.update({u : {coords : node_id}})
        node_id += 1
        txt += template % (node_id, "SOURCE", N * O, x, y, x, y, '', 1, '')
        u = "O_SOURCE"
        try:
            export_u_counts[u].update({coords : node_id})
        except:
            export_u_counts.update({u : {coords : node_id}})
        node_id += 1
        txt += template % (node_id, "SINK", 1, x, y, x, y, '', 2, '')
        u = "CLK_SINK"
        try:
            export_u_counts[u].update({coords : node_id})
        except:
            export_u_counts.update({u : {coords : node_id}})
        node_id += 1

        ptc = -1
        for p in range(0, cluster_inputs + N * O + 1):
            ptc += 1
            txt += template % (node_id, ('I' if p < cluster_inputs  or p >= cluster_inputs + N * O  else 'O') + "PIN", 1,\
                               x, y, x, y, "side=\"LEFT\" ", ptc, '')
            u = u_counts.get(p)
            try:
                export_u_counts[u].update({coords : node_id})
            except:
                export_u_counts.update({u : {coords : node_id}})
            node_id += 1

    #Tracks
    h_tracks = sorted([u for u in G if G.node[u]["node_type"] == "h_track"], key = lambda t : (t.split('H', 1)[1], t))
    v_tracks = sorted([u for u in G if G.node[u]["node_type"] == "v_track"], key = lambda t : (t.split('V', 1)[1], t))

    #NOTE: First group the >>potential_sb_edges<< by their target node, as this is where they should be inserted.
    chanx_width = get_chan_width(H)
    chany_width = get_chan_width(V)
    potential_edges = {}
    potential_edge_indices = {}
    for u in G:
        if u.startswith("potential_edge"):
            potential_edge_indices.update({u : -1})
            t = u.split("__")[-1]
            try:
                potential_edges[t].append(u)
            except:
                potential_edges.update({t : [u]})

    for ind, u in enumerate(sorted(potential_edge_indices)):
        potential_edge_indices[u] = ind

    min_x = io_crop
    min_y = io_crop
    max_x = grid_w - 1 - io_crop
    max_y = grid_h - 1 - io_crop

    trim_x = lambda x : min(max(min_x, x), max_x)
    trim_y = lambda y : min(max(min_y, y), max_y)

    extended_grid = list(grid.keys())
    for x in range(-1 * max([h[0] for h in H]), 0):
        for y in set([coord[1] for coord in grid]):
            extended_grid.append((x, y))
    for x in range(1, max([h[0] for h in H])):
        for y in set([coord[1] for coord in grid]):
            extended_grid.append((grid_w - 1 + x, y))
    for y in range(-1 * max([v[0] for v in V]), 0):
        for x in set([coord[0] for coord in grid]):
            extended_grid.append((x, y))
    for y in range(1, max([v[0] for v in V])):
        for x in set([coord[0] for coord in grid]):
            extended_grid.append((x, grid_h - 1 + y))

    valid_x = set([coord[0] for coord in grid])
    valid_y = set([coord[1] for coord in grid])
    for coords in sorted(extended_grid):
        x, y = coords
        if y in valid_y:
            ptc = -1
            visited = set()
            for u in h_tracks:
                d = u.split('_')[-2]
                L = int(u.split('_')[2][1:])
                if L in visited:
                    ptc += 1
                else:
                    ptc = 0
                    for L_visited in visited: 
                        ptc += len([t for t in h_tracks if "H%d" % L_visited in t]) * L_visited
                    ptc += (x % L) * len([t for t in h_tracks if "H%d" % L in t])
                    visited.add(L)
                seg_id = seg_ids[u.split('_')[2]]
                seg_decl = 3 * indent + "<segment segment_id=\"%d\"/>\n" % seg_id
                if d == 'L':
                    track_type = "CHANX\" direction=\"DEC_DIR"
                    xlow = x - L + 1
                    xhigh = x
                    ylow = yhigh = y
                    if xlow > max_x or xhigh < min_x:
                        continue
                    #if trim_x(xhigh) != xhigh:
                    #    continue
                    xlow = trim_x(xlow)
                    xhigh = trim_x(xhigh)
                    if y > max_y:
                        continue
                    mux_type = u.split('_')[2]
                    if y == 0:
                        try:
                            io_fanin_dict[(xlow, y)].append((node_id, "cb"))
                        except:
                            io_fanin_dict.update({(xlow, y) : [(node_id, "cb")]})
                        try:
                            io_fanout_dict[(xhigh, y)].append((node_id, mux_type))
                        except:
                            io_fanout_dict.update({(xhigh, y) : [(node_id, mux_type)]})
                    elif y == max_y:
                        try:
                            io_fanin_dict[(xlow, y + 1)].append((node_id, "cb"))
                        except:
                            io_fanin_dict.update({(xlow, y + 1) : [(node_id, "cb")]})
                        try:
                            io_fanout_dict[(xhigh, y + 1)].append((node_id, mux_type))
                        except:
                            io_fanout_dict.update({(xhigh, y + 1) : [(node_id, mux_type)]})
                elif d == 'R':
                    track_type = "CHANX\" direction=\"INC_DIR"
                    xlow = x
                    xhigh = x + L - 1
                    ylow = yhigh = y
                    if xlow > max_x or xhigh < min_x:
                        continue
                    #if trim_x(xlow) != xlow:
                    #    continue
                    xlow = trim_x(xlow)
                    xhigh = trim_x(xhigh)
                    if y > max_y:
                        continue
                    mux_type = u.split('_')[2]
                    if y == 0:
                        try:
                            io_fanin_dict[(xhigh, y)].append((node_id, "cb"))
                        except:
                            io_fanin_dict.update({(xhigh, y) : [(node_id, "cb")]})
                        try:
                            io_fanout_dict[(xlow, y)].append((node_id, mux_type))
                        except:
                            io_fanout_dict.update({(xlow, y) : [(node_id, mux_type)]})
                    elif y == max_y:
                        try:
                            io_fanin_dict[(xhigh, y + 1)].append((node_id, "cb"))
                        except:
                            io_fanin_dict.update({(xhigh, y + 1) : [(node_id, "cb")]})
                        try:
                            io_fanout_dict[(xlow, y + 1)].append((node_id, mux_type))
                        except:
                            io_fanout_dict.update({(xlow, y + 1) : [(node_id, mux_type)]})
                try:
                    export_u_counts[u].update({coords : node_id})
                except:
                    export_u_counts.update({u : {coords : node_id}})
    
                txt += template % (node_id, track_type, 1, xlow, ylow, xhigh, yhigh, '', ptc, seg_decl)
                node_id += 1
                for potential_edge in potential_edges.get(u, []):
                    pptc = potential_edge_indices[potential_edge]
                    try:
                        export_u_counts[potential_edge].update({coords : node_id})
                    except:
                        export_u_counts.update({potential_edge : {coords : node_id}})
                    seg_id = seg_ids[lut_canonical_potential_edge(potential_edge)]
                    seg_decl = 3 * indent + "<segment segment_id=\"%d\"/>\n" % seg_id
                    txt += template % (node_id, track_type, 1, xlow if d == 'R' else xhigh, ylow,\
                                       xlow if d == 'R' else xhigh, yhigh, '', pptc + chanx_width, seg_decl)
                    node_id += 1

        if x < 0 or x > max_x:
            continue

        ptc = -1
        visited = set()
        for u in v_tracks:
            d = u.split('_')[-2]
            L = int(u.split('_')[2][1:])
            if "_tap" in u:
                tap = int(u.split('_')[-1])
                d = u.split("_tap")[0].split('_')[-2]
                L = L if not SEPARATE_TAPS else (tap_phi if tap > 0 else L - (tap_M * tap_phi - 1))
                seg_id = seg_ids[u.split('_')[2] + ("_tap_%d" % tap)]
            else:
                seg_id = seg_ids[u.split('_')[2]]
            if L in visited:
                ptc += 1
            else:
                ptc = 0
                for L_visited in visited: 
                    ptc += len([t for t in v_tracks if "V%d" % L_visited in t]) * L_visited
                ptc += (y % L) * len([t for t in v_tracks if "V%d" % L in t])
                visited.add(L)
            seg_decl = 3 * indent + "<segment segment_id=\"%d\"/>\n" % seg_id
            if d == 'D':
                track_type = "CHANY\" direction=\"DEC_DIR"
                xlow = xhigh = x
                ylow = y - L + 1
                yhigh = y
                if ylow > max_y or yhigh < min_y:
                    continue
                #if trim_y(yhigh) != yhigh:
                #    continue
                ylow = trim_y(ylow)
                yhigh = trim_y(yhigh)
                if x > max_x:
                    continue
                mux_type = u.split('_')[2] + "_tap_0"
                if x == 0:
                    for tap in range(0, tap_M):
                        if ylow + tap >= grid_h - 1:
                            break
                        try:
                            io_fanin_dict[(x, ylow + tap)].append((node_id, "cb"))
                        except:
                            io_fanin_dict.update({(x, ylow + tap) : [(node_id, "cb")]})
                    try:
                        io_fanout_dict[(x, yhigh)].append((node_id, mux_type))
                    except:
                        io_fanout_dict.update({(x, yhigh) : [(node_id, mux_type)]})
                elif x == max_x:
                    for tap in range(0, tap_M):
                        if ylow + tap >= grid_h - 1:
                            break
                        try:
                            io_fanin_dict[(x + 1, ylow + tap)].append((node_id, "cb"))
                        except:
                            io_fanin_dict.update({(x + 1, ylow + tap) : [(node_id, "cb")]})
                    try:
                        io_fanout_dict[(x + 1, yhigh)].append((node_id, mux_type))
                    except:
                        io_fanout_dict.update({(x + 1, yhigh) : [(node_id, mux_type)]})
            elif d == 'U':
                track_type = "CHANY\" direction=\"INC_DIR"
                xlow = xhigh = x
                ylow = y
                yhigh = y + L - 1
                if ylow > max_y or yhigh < min_y:
                    continue
                #if trim_y(ylow) != ylow:
                #    continue
                ylow = trim_y(ylow)
                yhigh = trim_y(yhigh)
                if x > max_x:
                    continue
                mux_type = u.split('_')[2] + "_tap_0"
                if x == 0:
                    for tap in range(0, tap_M):
                        if yhigh - tap <= 0:
                            break
                        try:
                            io_fanin_dict[(x, yhigh - tap)].append((node_id, "cb"))
                        except:
                            io_fanin_dict.update({(x, yhigh - tap) : [(node_id, "cb")]})
                    try:
                        io_fanout_dict[(x, ylow)].append((node_id, mux_type))
                    except:
                        io_fanout_dict.update({(x, ylow) : [(node_id, mux_type)]})
                elif x == max_x:
                    for tap in range(0, tap_M):
                        if yhigh - tap <= 0:
                            break
                        try:
                            io_fanin_dict[(x + 1, yhigh - tap)].append((node_id, "cb"))
                        except:
                            io_fanin_dict.update({(x + 1, yhigh - tap) : [(node_id, "cb")]})
                    try:
                        io_fanout_dict[(x + 1, ylow)].append((node_id, mux_type))
                    except:
                        io_fanout_dict.update({(x + 1, ylow) : [(node_id, mux_type)]})
            try:
                export_u_counts[u].update({coords : node_id})
            except:
                export_u_counts.update({u : {coords : node_id}})

            txt += template % (node_id, track_type, 1, xlow, ylow, xhigh, yhigh, '', ptc, seg_decl)
            node_id += 1
            for potential_edge in potential_edges.get(u, []):
                pptc = potential_edge_indices[potential_edge]
                try:
                    export_u_counts[potential_edge].update({coords : node_id})
                except:
                    export_u_counts.update({potential_edge : {coords : node_id}})
                seg_id = seg_ids[lut_canonical_potential_edge(potential_edge)]
                seg_decl = 3 * indent + "<segment segment_id=\"%d\"/>\n" % seg_id
                txt += template % (node_id, track_type, 1, xlow, ylow if d == 'U' else yhigh, xhigh,\
                                   ylow if d == 'U' else yhigh, '', pptc + chany_width, seg_decl)
                node_id += 1

    txt += "</rr_nodes>\n"

    return txt, export_u_counts, io_fanin_dict, io_fanout_dict
##########################################################################

##########################################################################
def export_rr_edges(G, u_counts, io_fanin_dict, io_fanout_dict):
    """Exports the RR graph edges in the VTR8 format. 

    Parameters
    ----------
    G : nx.MultiDiGraph
        The routing-resource graph.
    u_counts : Dict[str, Dict[Tuple[int], int]]
        Mapping between the RR-graph nodes in the static form (G)
        and the tile coordinates and node ids.

    Returns
    -------
    str
        Text of the tags.
    """
   
    txt = indent + "<rr_edges>\n"
    beg = 2 * indent + "<edge src_node=\""
    lines = set()

    #------------------------------------------------------------------------#
    def terminates_early(u, coords, u_counts):
        """Determines if a wire terminates early. If so, it finds its reflection.
    
        Parameters
        ----------
        u : str
            Node name in the static graph.
        coords : Tuple[int]
            Coordinates of the replica in the periodic graph.
        u_counts : Dict[str, Dict[Tuple[int], int]]
            Mapping between the RR-graph nodes in the static form (G)
            and the tile coordinates and node ids.
       
        Returns
        -------
        int
            The reflected child. If it does not exist, None.
        """
   
        #NOTE: Remember that channels are created above and to the right of blocks,
        #which determines the hanging conditions.

        #TODO: Check if it works for separated taps as well!
    
        x, y = coords
    
        L = int(u.split('_')[2][1:])
        d = u.split("_tap")[0].split('_')[-2]
        
        if d == 'R':
            hang = x + (L - 1) - (grid_w - 1 - io_crop)
            if hang <= 0:
                return None
            return u_counts[u.replace("_R_", "_L_")][(x + (L - 1), y)]
        elif d == 'L':
            hang = L - x
            if hang <= 0:
                return None
            return u_counts[u.replace("_L_", "_R_")][(x - (L - 1), y)]
        elif d == 'U':
            hang = y + (L - 1) - (grid_h - 1 - io_crop)
            if hang <= 0:
                return None
            return u_counts[u.replace("_U_", "_D_")][(x, y + (L - 1))]
        elif d == 'D':
            hang = L - y
            if hang <= 0:
                return None
            return u_counts[u.replace("_D_", "_U_")][(x, y - (L - 1))]
        else:
            print "Unknown direction", d, u
            raise ValueError
    #------------------------------------------------------------------------#

    reflections = set()
    reflected = set()
    for u in G:
        for coords in u_counts[u]:
            x, y = coords
            uc = u_counts[u][coords]
            vs = list(G[u])
            u_proxy = u
            if G.node[u]["node_type"] in ("h_track", "v_track") and not "potential_edge" in u: 
                reflection = terminates_early(u, coords, u_counts)
                if reflection is not None:
                    mux_id = " switch_id=\"%d\"/>\n" % mux_ids["__vpr_delayless_switch__"]
                    lines.add("%s%d\" sink_node=\"%d\"%s" % (beg, uc, reflection, mux_id))
                    reflections.add(reflection)
                    reflected.add(uc)
                    continue
   
            if SEPARATE_TAPS:
                L = None
                try:
                    L = int(u.split('_')[2][1:])
                    if "_tap" in u:
                        tap = int(u.split('_')[-1])
                        d = u.split("_tap")[0].split('_')[-2]
                        L = tap_phi if tap > 0 else L - (tap_M * tap_phi - 1)
                        if (d == 'U' and coords[1] + L >= grid_h - 1 and tap != tap_M - 1)\
                        or (d == 'D' and coords[1] - L <= 1 and tap != tap_M - 1):
                            #NOTE: This is a prematurely terminating wire, so
                            #we should reconnect it.
                            u_proxy = u.replace("_tap_%d" % tap, "_tap_%d" % (tap_M - 1))
                            vs = list(G[u_proxy])
                except:
                    pass
            for v in vs:
                for e in G[u_proxy][v]:
                    attrs = G[u_proxy][v][e]
                    offset = list(attrs.get("offset", (0, 0)))
                    if G.node[u]["node_type"] == "h_track":
                        if x + offset[0] > grid_w - 1 - io_crop and G.node[v]["node_type"].endswith("v_track"):
                            offset[0] = grid_w - 1 - io_crop - x
                        elif x + offset[0] < 0 and G.node[v]["node_type"].endswith("v_track"):
                            offset[0] = - x
                            #NOTE: There is a (0, y) vertical channel, but not (grid_w - 1, y).
                    elif G.node[u]["node_type"] == "v_track":
                        if y + offset[1] > grid_h - 1 - io_crop and G.node[v]["node_type"].endswith("h_track"):
                            offset[1] = grid_h - 1 - io_crop - y
                        elif y + offset[1] < 0 and G.node[v]["node_type"].endswith("h_track"):
                            offset[1] = - y
                            #NOTE: There is a (x, 0) horizontal channel, but not (x, grid_h - 1).

                    mux_id = " switch_id=\"%d\"/>\n" % mux_ids[attrs["mux_type"]]
                    try:
                        vc = u_counts[v][(x + offset[0], y + offset[1])]
                    except:
                        continue
                    if uc == vc:
                        #NOTE: This could happen at the grid boundary.
                        continue
                    if HUMAN_READABLE:
                        u_str = u + '_' + str(coords)
                        v_str = v + '_' + str((coords[0] + attrs.get("offset", (0, 0))[0],\
                                               coords[1] + attrs.get("offset", (0, 0))[1]))
                        lines.add("%s%s\" sink_node=\"%s\"%s" % (beg, u_str, v_str, mux_id))
                    else:
                        lines.add("%s%d\" sink_node=\"%d\"%s" % (beg, uc, vc, mux_id))

    for coords in sorted(io_fanin_dict):
        for i in range(0, IO_CAPACITY):
            v = "io_%d_opad_in" % i
            vc = u_counts[v][coords]
            for u in io_fanin_dict[coords]:
                uc, mux = u
                if uc in reflected:
                    continue
                mux_id = " switch_id=\"%d\"/>\n" % mux_ids[mux]
                lines.add("%s%d\" sink_node=\"%d\"%s" % (beg, uc, vc, mux_id))
            v = "io_%d_clk_in" % i
            vc = u_counts[v][coords]
            for u in io_fanin_dict[coords]:
                uc, mux = u
                mux_id = " switch_id=\"%d\"/>\n" % mux_ids[mux]
                lines.add("%s%d\" sink_node=\"%d\"%s" % (beg, uc, vc, mux_id))

    for coords in sorted(io_fanout_dict):
        for i in range(0, IO_CAPACITY):
            for v in io_fanout_dict[coords]:
                vc, mux = v
                if vc in reflections:
                    continue
                mux_id = " switch_id=\"%d\"/>\n" % mux_ids[mux]
                u = "io_%d_ipad_out" % i
                uc = u_counts[u][coords]
                lines.add("%s%d\" sink_node=\"%d\"%s" % (beg, uc, vc, mux_id))

    #Now remove any remaining reflections on edge-splitters from the edge list.
    get_attr = lambda line, attr : line.split("%s=\"" % attr)[1].split('"')[0]
    rm_set = set()
    for line in lines:
        if int(get_attr(line, "sink_node")) in reflections or int(get_attr(line, "src_node")) in reflected:
            if int(get_attr(line, "switch_id")) != mux_ids["__vpr_delayless_switch__"]:
                rm_set.add(line)
    for line in rm_set:
        lines.remove(line)
            

    sort_key = lambda l : (int(l.split("=\"")[1].split('"')[0]), int(l.split("=\"")[2].split('"')[0]))\
                          if not HUMAN_READABLE else l

    txt += ''.join(sorted(lines, key = sort_key))
    txt += indent + "</rr_edges>\n"

    return txt
##########################################################################            

##########################################################################
def generate_rr_graph(make_sb_clique = False):
    """Top-level function generating the entire RR-graph.

    Parameters
    ----------
    make_sb_clique : Optional[bool], default = False
        Local equivalent to args.make_sb_clique.
        It does not make sense to construct the clique
        while padding is in progress.    

    Returns
    -------
    nx.MultiDiGraph
        The routing-resource graph.
    Dict[Tuple[int], str]
        A dictionary of block types, indexed by the grid coordinates.
    """

    G = nx.MultiDiGraph()
    grid = generate_grid(grid_w, grid_h)
    add_io_pins(G)
    add_cluster_pins(G)
    add_io_sinks_and_sources(G)
    add_cluster_sinks_and_sources(G)
    compose_channels(G, H, V)
    add_taps(G)
    add_clb_to_sb(G)
    add_sb_to_sb(G, make_clique = make_sb_clique)

    return G, grid
##########################################################################

##########################################################################
def export_rr_graph(G, grid, filename, potential_edge_delays = None):
    """Exports the RR-graph in the VTR8 RR-graph format.

    Parameters
    ----------
    G :  nx.MultiDiGraph
        The routing-resource graph.
    grid : Dict[Tuple[int], str]
        A dictionary of block types, indexed by the grid coordinates.
    potential_edge_switches : Optional[Dict[str, float]], default = None
        A dictionary of delays assigned to the potential edges.

    Returns
    -------
    None
    """

    header = "<rr_graph tool_name=\"vpr\" tool_version=\"8.0.0+unkown\""\
           +  " tool_comment=\"Generated from arch file %s\">\n" % args.arc_name
    footer = "</rr_graph>"

    td_dict = {}
    cb_delay = default_cb_delay
    inherit = None
    try:
        inherit = args.change_grid_dimensions
    except:
        pass
    if inherit is None:
        td_dict = spice_all_wires(G)
        cb_delay = 0.5 * (td_dict["cb_h"] + td_dict["cb_v"])
    else:
        td_dict = read_delays_from_arc(inherit)
        cb_delay = td_dict["cb"]

    #Now remove the delay-measurement-only edges.
    if MAKE_SB_CLIQUE:
        rem_list = []
        for u,v, attrs in G.edges(data = True):
            if attrs.get("delay_operating_point_only", False):
                rem_list.append((u, v))
        for u, v in rem_list:
            G.remove_edge(u, v)

    with open(filename, "w") as outf:
        potential_edges = sorted(set([attrs["seg"] for u, attrs in G.nodes(data = True)\
                                     if u.startswith("potential_edge")])) if MAKE_SB_CLIQUE else None
        txt = header + export_chan_tags(H, V, grid_w, grid_h,\
                                        potential_edge_no = N * len(potential_edges) if MAKE_SB_CLIQUE else 0)

        #Determine which potential delays will be used (default, measured, or inherited).
        potential_edge_fwd = potential_edges
        if potential_edge_delays is not None:
            potential_edge_fwd = potential_edge_delays
        if inherit is not None:
            if not potential_edges:
                potential_edge_fwd = {}
            else:
                potential_edge_fwd = {e : td_dict[e] for e in potential_edges}

        outf.write(txt)
        global mux_ids
        txt, mux_ids = export_switches(H, V, cb_delay, td_dict, potential_edges = potential_edge_fwd)
        outf.write(txt)
        global seg_ids
        txt, seg_ids = export_segments(H, V, potential_edges = potential_edge_fwd)
        outf.write(txt)
        txt = export_blocks()
        outf.write(txt)
        txt = export_grid(grid)
        outf.write(txt)
        txt, counts, io_fanin_dict, io_fanout_dict = export_rr_nodes(G, grid)
        outf.write(txt)
        txt = export_rr_edges(G, counts, io_fanin_dict, io_fanout_dict) + footer
        outf.write(txt)
    if COMPRESS_RR:
        os.system("lz4 --rm %s %s.lz4" % (filename, filename)) 

    if True:#inherit is None
        fill_in_template(G, cb_delay, td_dict, potential_edges = potential_edge_fwd)
    else:
        with open(inherit, "r") as inf:
            lines = inf.readlines()
        txt = ""
        for line in lines:
            if "<fixed_layout" in line:
                words = line.split()
                words[2] = "width=\"%d\"" % grid_w
                words[3] = "height=\"%d\">" % grid_h 
                txt += 2 * indent + ' '.join(words) + "\n"
            else:
                txt += line
        with open(args.arc_name, "w") as outf:
            outf.write(txt)

    placement_delay_matrix_filename = args.arc_name.rsplit(".xml", 1)[0] + "_placement_delay.matrix" 
    generate_optimal_placement_delay_matrix(cb_delay, td_dict, placement_delay_matrix_filename)
##########################################################################

##########################################################################
def export_fpga_layout():
    """Exports the grid layout information in VTR8 xml format.

    Parameters
    ----------
    None

    Returns
    -------
    str
        Text of the tags.
    """

    io_priority = 100
    clb_priority = 10
    mult_priority = mem_priority = 20

    txt = indent + "<layout>\n"
    txt += 2 * indent + "<fixed_layout name=\"fix\" width=\"%d\" height=\"%d\">\n"\
                     % (grid_w, grid_h)
    txt += 3 * indent + "<perimeter type=\"io\" priority=\"%d\"/>\n" % io_priority

    if cut_corners:
        txt += 3 * indent + "<corners type=\"EMPTY\" priority=\"%d\"/>\n" % (io_priority + 1)
    if TOP_BOTTOM_IO:
        txt += 3 * indent + "<col type=\"EMPTY\" startx=\"0\" starty=\"1\" priority=\"%d\"/>\n" % (io_priority + 2)
        txt += 3 * indent + "<col type=\"EMPTY\" startx=\"%d\" starty=\"1\" priority=\"%d\"/>\n" % (grid_w - 1, io_priority + 2)

    if INSERT_EMPTY_RING:
        max_h = max([h[0] for h in H])
        max_v = max([v[0] for v in V])
        ring_xs = [x for x in range(1, max_h + 1)] + [x for x in range(grid_w - 1 - max_h, grid_w - 1)]
        ring_ys = [y for y in range(1, max_v + 1)] + [y for y in range(grid_h - 1 - max_v, grid_h - 1)]
        for x in ring_xs:
            for y in range(1, grid_h - 1):
                txt += 3 * indent + "<single type=\"empty_clb\" x=\"%d\" y=\"%d\" priority=\"%d\"/>\n"\
                                  % (x, y, clb_priority + 1)
        for y in ring_ys:
            for x in range(1, grid_w -1):
                if x in ring_xs:
                    continue
                txt += 3 * indent + "<single type=\"empty_clb\" x=\"%d\" y=\"%d\" priority=\"%d\"/>\n"\
                                  % (x, y, clb_priority + 1)

    txt += 3 * indent + "<fill type=\"clb\" priority=\"%d\"/>\n" % clb_priority

    if mult_col["freq"]:
        txt += 3 * indent + "<col type=\"mult_36\" startx=\"%d\" starty=\"%d\" repeatex=\"%d\" priority=\"%d\"/>\n"\
                          % (mult_col["start"][0], mult_col["start"][1], mult_col["freq"], mult_priority)
        txt += 3 * indent + "<col type=\"EMPTY\" startx=\"%d\" starty=\"%d\" repeatex=\"%d\" priority=\"%d\"/>\n"\
                          % (mult_col["start"][0], mult_col["start"][1], mult_col["freq"], mult_priority - 1)
    if mem_col["freq"]:
        txt += 3 * indent + "<col type=\"memory\" startx=\"%d\" starty=\"%d\" repeatex=\"%d\" priority=\"%d\"/>\n"\
                          % (mem_col["start"][0], mem_col["start"][1], mem_col["freq"], mem_priority)
        txt += 3 * indent + "<col type=\"EMPTY\" startx=\"%d\" starty=\"%d\" repeatex=\"%d\" priority=\"%d\"/>\n"\
                          % (mem_col["start"][0], mem_col["start"][1], mem_col["freq"], mem_priority - 1)

    txt += 2 * indent + "</fixed_layout>\n"
    txt += indent + "</layout>\n"

    return txt
##########################################################################

##########################################################################
def export_fpga_switches(H, V, td_cb = default_delay, td_sb = {}, potential_edges = None):
    """Exports the switch information in the VTR8 xml format.

    Parameters
    ----------
    H : List[Tuple[int]]
        The list of track lengths and occurrences per BLE of each,
        in the horizontal channel.
    V : List[Tuple[int]]
        The list of track lengths and occurrences per BLE of each,
        in the vertical channel.
    td_cb : float
        Delay of the connection block multiplexer (and SB-CB wires)
    td_sb : Dict[str, float]
        Delays of all the switch block multiplexers, together with the
        wires they are driving, and the typical loading by other muxes.
    potential_edges : Optional[List[str]], default = None
        List of potential edges obtained with the appropriate >>make_clique<< switches turned on.
        Optionally also a dictionary of their delays can be passed.

    Returns
    -------
    str
        Text of the tags
    """

    template = 2 * indent + "<switch type=\"mux\" name=\"%s\" R=\"0\" Cin=\"0\" Cout=\"0\" Tdel=\"%g\""\
                          + " mux_trans_size=\"0\" buf_size=\"0\"/>\n"

    txt = indent + "<switchlist>\n"
    for h in H:
        name = "H%d" % h[0]
        txt += template % (name, td_sb.get(name, 0.5 * default_delay * (h[0] + 1)))
    for v in V:
        for tap in range(0, tap_M if SEPARATE_TAPS else 1):
            name = "V%d_tap_%d" % (v[0], tap)
            txt += template % (name, td_sb.get(name, 0.5 * default_delay * (v[0] + 1)))
    txt += template % ("cb", td_cb)

    if potential_edges is not None:
        potential_edge_delay = 1e-12
        #NOTE:VPR lookahead may break with zero delays. Be cautious.
        for u in potential_edges:
            if isinstance(potential_edges, dict):
                potential_edge_delay = potential_edges[u]
            txt += template % (u, potential_edge_delay)

    txt += indent + "</switchlist>\n"

    return txt
##########################################################################

##########################################################################
def export_fpga_segments(H, V, potential_edges = None):
    """Exports the segment information in the VTR8 xml format.

    Parameters
    ----------
    H : List[Tuple[int]]
        The list of track lengths and occurrences per BLE of each,
        in the horizontal channel.
    V : List[Tuple[int]]
        The list of track lengths and occurrences per BLE of each,
        in the vertical channel.
    potential_edges : Optional[List[str]], default = None
        List of potential edges obtained with the appropriate >>make_clique<< switches turned on.

    Returns
    -------
    str
        Text of the tags
    """

    template = 2 * indent + "<segment freq=\"1.0\" length=\"%d\" type=\"unidir\""\
             + " Rmetal=\"0\" Cmetal=\"0\" name=\"%s\">\n"
    template += 3 * indent + "<mux name=\"%s\"/>\n"
    template += 3 * indent + "<sb type=\"pattern\">1%s1</sb>\n"
    template += 3 * indent + "<cb type=\"pattern\">1%s</cb>\n"
    template += 2 * indent + "</segment>\n"

    txt = indent + "<segmentlist>\n"

    for h in H:
        L = h[0]
        name = "H%d" % L
        txt += template % (L, name, name, ' ' + (L - 1) * "0 ", ' ' + (L - 2) * "0 "  + '1' if L > 1 else '')
    for v in V:
        if SEPARATE_TAPS:
            for tap in range(0, tap_M):
                L = v[0]
                name = "V%d_tap_%d" % (L, tap)
                L = tap_phi if tap > 0 else L - (tap_M * tap_phi - 1)
                txt += template % (L, name, name, ' ' + (L - 1) * "0 ", ' ' + (L - 2) * "0 "  + '1' if L > 1 else '')
        else:
            L = v[0]
            name = "V%d_tap_0" % L
            txt += template % (L, name, name, ' ' + (L - 1) * "0 ",\
                              ' ' + (L - 1 - tap_M) * "0 "  + (min(tap_M, L - 1) * "1 ") if L > 1 else '')

    if potential_edges is not None:
        potential_edge_logical_length = 0
        #NOTE: VPR lookahead may break with length-0 wires. Be cautious.
        for u in sorted(potential_edges):
            u_txt = template % (potential_edge_logical_length, u, u, ' ', '')
            u_lines = u_txt.splitlines()
            txt += "\n".join(u_lines[:2] + [u_lines[-1]]) + "\n"

    txt += indent + "</segmentlist>\n"

    return txt
##########################################################################

##########################################################################
def export_mux_sizes(G):
    """Exports the multiplexer sizes.

    Parameters
    ----------
    G : nx.MultiDiGraph
        The routing-resource graph.

    Returns
    -------
    List[int]
        Connection block multiplexer sizes.
    List[int] 
        Switch block multiplexer sizes.
    """

    nodes = [u for u in G.nodes() if u.startswith("ble_%d_" % NEUTRAL_BLE) and G.in_degree(u) and not "_o_" in u]

    cb_nodes = [u for u in nodes if "cb_out" in u]
    sb_nodes = [u for u in nodes if not "cb_out" in u and (not "_tap" in u or "_tap_0" in u)]

    cb_sizes = {u : sum([len(G[p][u]) for p in G.pred[u] if not "io" in p and not "potential_edge" in p]) for u in cb_nodes}
    sb_sizes = {u : sum([len(G[p][u]) for p in G.pred[u] if not "io" in p and not "potential_edge" in p]) for u in sb_nodes}

    return cb_sizes, sb_sizes
##########################################################################

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

##########################################################################
def stack_muxes(G, get_pins = False):
    """Stacks the routing multiplexers.

    Parameters
    ----------
    G : nx.MultiDiGraph
        The routing-resource graph.
    get_pins : Optional[bool], default = False
        Specifies if mux pins should be returned too.

    Returns
    -------
    int
        Total routing mux width in fin pitches.
    Dict[str, Tiple[int]]
        A dictionary of 
    """

    crossbar_muxes = {"crossbar%d" % i : crossbar_mux_size for i in range(0, K)}
    cb_muxes, sb_muxes = export_mux_sizes(G)

    all_sizes = {}
    all_sizes.update(crossbar_muxes)
    all_sizes.update(cb_muxes)
    all_sizes.update(sb_muxes)

    muxes = list(sorted(crossbar_muxes, key = lambda m : crossbar_muxes[m], reverse = True))
    muxes += list(sorted(cb_muxes, key = lambda m : cb_muxes[m], reverse = True))
    if SB_MUX_ORDER is None:
        muxes += list(sorted(sb_muxes, key = lambda m : sb_muxes[m], reverse = True))
    else:
        muxes += copy.copy(SB_MUX_ORDER)

    cols = [[0, []]]

    get_w = lambda row_cnt : 21 + row_cnt
    #Width in fin pitches

    get_h = lambda row_cnt, col_cnt : 2 * max(row_cnt, col_cnt)
    #Height in gate pitches

    #-------------------------------------------------------------------------#
    def get_driver(mux, w):
        """Returns the gate pitches required by the driver.

        Parameters
        ----------
        mux : str
            Multiplexer identifier.
        w : int
            Multiplexer width in fin pitches.
        
        Returns
        -------
        int
            Number of gate pitches to add to mux height.
        """

        get_height = lambda driver : int(math.ceil(((driver[0] * 2 + 1 + driver[1] * 2 + 1 + 1) / float(w))))

        if "crossbar" in mux:
            return 0
        if "cb_out" in mux:
            return get_height(local_driver)
        L = int(mux.split('_')[2][1:])
        if 'H' in mux:
            return get_height(H_drivers.get(L, (0, 0)))
        if 'V' in mux:
            return get_height(V_drivers.get(L, (0, 0)))

        print "Illegal mux name."
        raise ValueError
    #-------------------------------------------------------------------------#

    pins = {}
    cur_x = 0
    cur_y = 0
    for mux in muxes:
        size = all_sizes[mux]
        rcnt, ccnt = get_mux_dimensions(size, lut_height)
        w = get_w(rcnt)
        h = get_h(rcnt, ccnt) + get_driver(mux, w)
        if cols[-1][0] + h > lut_height:
            cur_x -= max([added_mux[0] for added_mux in cols[-1][1]])
            cols.append([0, []])
            cur_y = 0
        pins.update({mux : {'i' : (cur_x - 0.5 * w, cur_y), 'o' : (cur_x, cur_y + 0.5 * h)}})
        cur_y += h
        cols[-1][0] += h
        cols[-1][1].append((w, h))

    if get_pins:
        return pins, all_sizes

    return sum([max(mux[0] for mux in col[1]) for col in cols])
##########################################################################

##########################################################################
def get_tile_dimensions(G):
    """Returns physical dimensions of a tile in nanometers.

    Parameters
    ----------
    G : nx.MultiDiGraph
        The routing-resource graph.

    Returns
    -------
    int
        Tile width.
    int
        Tile height.
    """

    return (lut_width + stack_muxes(G)) * FP, lut_height * N * GP
##########################################################################

##########################################################################
def get_metal_dimensions():
    """Returns the width and the height needed to trace all metal.

    Parameters
    ----------
    None

    Returns
    -------
    Tuple(int)
        Width and height needed by metal.
    """

    width = 0
    for v in V:
        width += 2 * v[0] * v[1] * N

    height = 0
    for h in H:
        height += 2 * h[0] * h[1] * N

    return width * MyP, height * MyP
##########################################################################

##########################################################################
def pad_LEN1(G, max_mux_width = 16):
    """Determines the number of LEN-1 wires that can be added until
    either the tile dimensions no longer fit any tracks, or mux size
    increases beyond an allowed amount.

    Parameters
    ----------
    G : nx.MultiDiGraph
        The routing-resource graph.

    max_mux_width : Optional[int], default = 16
        Maximum allowed mux size.

    Returns
    -------
    Tuple(int)
        The number of H1 and V1 wires to be added.
    nx.MultiDiGraph
        The updated routing graph.

    Notes
    -----
    We assume one layer per direction.
    """

    tile_w, tile_h = get_tile_dimensions(G)
    metal_w, metal_h = get_metal_dimensions()

    global H
    global V
    #We assume that LEN-1 wires don't enter the original specification.

    K6N8_LUT4 = 8 * 4
    KN_LUT4 = N * 2 ** (K - 4)
    VL = max(1, K6N8_LUT4 / KN_LUT4)

    h1_exists = False
    for i, h in enumerate(H):
        if h[0] == 1:
            h1_exists = True
            break
    h = 0
    if h1_exists:
        h = H.pop(i)[1]
            
    H.append((1, h))

    v1_exists = False
    for i, v in enumerate(V):
        if v[0] == VL:
            v1_exists = True
            break
    v = 0
    if v1_exists:
        v = V.pop(i)[1]

    V.append((VL, v))

    #-------------------------------------------------------------------------#
    def get_largest_mux(G):
        """Returns the size of the largest mux, modulo the twists
        
        Parameters
        ----------
        None

        Returns
        -------
        int
            Largest mux size.
        """

        return 0

        cb_sizes, sb_sizes = export_mux_sizes(G)
        return max(max([cb_sizes[mux] for mux in cb_sizes]), max([sb_sizes[mux] for mux in sb_sizes])) 
    #-------------------------------------------------------------------------#

    if args.import_padding is not None:
        with open(args.import_padding, "r") as inf:
            lines = inf.readlines()
        for line in lines:
            if line.startswith("H1"):
                H[-1] = (1, H[-1][1] + int(line.split()[-1]))
            elif line.startswith("V1"):
                V[-1] = (VL, V[-1][1] + int(line.split()[-1]))
        G, grid = generate_rr_graph()
        metal_w, metal_h = get_metal_dimensions()
        tile_w, tile_h = get_tile_dimensions(G)
        print("Padding imported.")
    else:
        h_tracks = 1 #int(math.floor(float(tile_h - metal_h) / (2 * N * MyP)))
        while metal_h <= tile_h:
            H[-1] = (1, H[-1][1] + max(1, h_tracks))
            G, grid = generate_rr_graph()
            metal_w, metal_h = get_metal_dimensions()
            tile_w, tile_h = get_tile_dimensions(G)
            if get_largest_mux(G) > max_mux_width:
                break
        print("H added.")
    
        v_tracks = 1 #int(math.floor(float(tile_w - metal_w) / (2 * VL * N * MyP)))
        while metal_w <= tile_w:
            V[-1] = (VL, V[-1][1] + max(1, v_tracks))
            G, grid = generate_rr_graph()
            metal_w, metal_h = get_metal_dimensions()
            tile_w, tile_h = get_tile_dimensions(G)
            #v_tracks = int(math.floor(float(tile_w - metal_w) / (2 * VL * N * MyP)))
            if get_largest_mux(G) > max_mux_width:
                break
        print("V added.")

    #------------------------------------------------------------------------#
    def log_padding_results():
        """Logs the results of padding.

        Parameters
        ----------
        None
        
        Returns
        -------
        None
        """

        filename = args.arc_name.rsplit('.', 1)[0] + "_padding.log"

        txt = "Channel composition after padding:\n\n"
        for h in sorted(H):
            txt += "H%d %d\n" % (h[0], h[1])
        for v in sorted(V):
            txt += "V%d %d\n" % (v[0], v[1])

        active_w, active_h =  get_tile_dimensions(G)
        txt += "\nActive dimensions: %d X %d nm\n" % (active_w, active_h)
        metal_w, metal_h = get_metal_dimensions()
        txt += "Metal dimensions: %d X %d nm\n\n" % (metal_w, metal_h)

        cb_sizes, sb_sizes = export_mux_sizes(G)
        mux_sizes = cb_sizes
        mux_sizes.update(sb_sizes)
        size_indexed = {}
        for mux in mux_sizes:
            size = mux_sizes[mux]
            mux_id = mux.split("ble_%d_" % NEUTRAL_BLE)[1]
            try:
                size_indexed[size].append(mux_id)
            except:
                size_indexed.update({size : [mux_id]})

        txt += "Multiplexer sizes:\n\n"
        for size in sorted(size_indexed):
            txt += "%d: " % size
            for mux in sorted(size_indexed[size]):
                txt += "%s " % mux
            txt = txt[:-1] + "\n"

        with open(filename, "w") as outf:
            outf.write(txt[:-1])
    
        print("\n" + txt)
    #------------------------------------------------------------------------#

    log_padding_results()
        

    return H[-1][-1], V[-1][-1], G, grid
##########################################################################

##########################################################################
def meas_lut_access_delay(G):
    """Measures the LUT to SB connection access delay.

    Parameters
    ----------
    G : nx.MultiDiGraph
        The routing-resource graph.

    Returns
    -------
    nx.Graph
        Netlist for meaurement.    
    """

    lut_x = 2 ** 4 * 10 / 2.0
    lut_y = lut_height / (2 ** (K - 4) / 2.0)

    pins, all_sizes = stack_muxes(G, get_pins = True)
    for u, attrs in G.nodes(data = True):
        if attrs["node_type"] == "clb_out" and "ble_%d_" % NEUTRAL_BLE in u:
            break

    fanout = [v for v in pins if v in G[u]]
    fanout.sort(key = lambda u : (pins[u]['i'][1], pins[u]['i'][0]))
    fanout.sort(key = lambda m : abs(pins[m]['i'][0] - lut_x) * FP + abs(pins[m]['i'][1] - lut_y) * GP)

    all_xs = [pins[f]['i'][0] for f in fanout]
    avg_x = float(sum(all_xs)) / len(all_xs)

    target_mux = "t%d" % (len(fanout) / 2)

    #------------------------------------------------------------------------#
    def connect_targets(fanouts, pins, all_sizes, lut_x, lut_y):
        """Connects the targets into a lattice to approximate what a real router would do.

        Parameters
        ----------
        fanouts : List[str]
            The list of targets.
        pins : Dict[str, various]
            A dictionary of pin coordinates for all multiplexers.
        all_sizes : Dict[str, int]
            A dictionary of multiplexer sizes.
        lut_x : float
            Horizontal offset of the LUT output.
        lut_y : float
            Height of the LUT output.

        Returns
        -------
        nx.Graph()
            The lattice graph.
        """

        lattice = nx.Graph()
        cols = {}
        col_thr = 10
        #The number of fin pitches horizontal offset to be considered part of the same column.
        for fanout in fanouts:
            x, y = pins[fanout]['i']
            u = "t%d" % fanouts.index(fanout)
            if u == target_mux:
                u = 't'
            lattice.add_node(u, coords = (x, y), mux = True, size = all_sizes[fanout], potential_target = True)
            found = False
            for col in cols:
                if abs(x - col) <= col_thr:
                    cols[col].append((y, u))
                    found = True
                    break
            if not found:
                cols.update({x : [(y, u)]})

        for i, col in enumerate(sorted(cols, key = lambda c : abs(lut_x - c))):
            lattice.add_node("junction_%d" % i, coords = (col, lut_y), mux = False)
            cols[col].append((lut_y, "junction_%d" % i))
            if i == 0:
                lattice.add_node("wire_s", coords = (lut_x, lut_y), mux = False)
                lattice.add_edge("wire_s", "junction_0", layer = "Mx")
            else:
                lattice.add_edge("junction_%d" % i, "junction_%d" % (i - 1), layer = "Mx")
            cols[col].sort()
            for j, node in enumerate(cols[col][1:], 1):
                #By construction there must be at least 2 nodes in the column, one being the junction.
                u = node[1]
                v = cols[col][j - 1][1]
                lattice.add_edge(u, v, layer = "Mx")

        return lattice
    #------------------------------------------------------------------------#

    net = connect_targets(fanout, pins, all_sizes, lut_x, lut_y)
    node_mapping = {}

    for t, target in enumerate(fanout):
        u = "t%d" % t
        if u == target_mux:
            u = 't'
        node_mapping.update({target : u})

    #Determine the states of muxes.
    total_cols = 0
    for u in fanout:
        t = node_mapping[u]
        if t == target_mux:
            continue
        rcnt, ccnt = get_mux_dimensions(net.node[t]["size"], lut_height)
        total_cols += ccnt
    partial_prob = 1.0 / (float(total_cols) / len(fanout))
    partial_fanout = int(math.ceil(partial_prob * len(fanout)))

    partial_space = len(fanout) / partial_fanout
    for i, f in enumerate(fanout):
        m = node_mapping[f]
        if i % partial_space == 0:
            net.node[m]["state"] = "partial"
        else:
            net.node[m]["state"] = "off"

    net.node['t']["state"] = "on"

    return net 
##########################################################################

##########################################################################
def get_netlist(G, wire, source, get_cb_delay = False):
    """Assembles a netlist for Spice measurement of the given wire type.

    Parameters
    ----------
    G : nx.MultiDiGraph
        The routing-resource graph.
    wire : str
        Wire type.
    source : str
        Source node.
    get_cb_delay : Optional[bool], default = False
        Determines the position of the wire and the connection block and then calls
        >>local_wires.py<< to obtain the delay from the wire to a LUT input pin.
   
    Returns
    -------
    nx.Graph
        The netlist graph.
    """

    #------------------------------------------------------------------------#
    def connect_targets(fanouts, pins, all_sizes, wire_x, wire_y, tap = None):
        #FIXME: Almost a duplicate from meas_lut_access().
        """Connects the targets into a lattice to approximate what a real router would do.

        Parameters
        ----------
        fanouts : List[str]
            The list of targets.
        pins : Dict[str, various]
            A dictionary of pin coordinates for all multiplexers.
        all_sizes : Dict[str, int]
            A dictionary of multiplexer sizes.
        wire_x : float
            Horizontal offset of the global wire junction point.
        wire_y : float
            Height of the global wire junction point.
        tap : Optional[int], default = None
            Tap, if the wire is vertical.

        Returns
        -------
        nx.Graph()
            The lattice graph.
        """

        lattice = nx.Graph()
        cols = {}
        col_thr = 10
        #The number of fin pitches horizontal offset to be considered part of the same column.
        for fanout in fanouts:
            x, y = pins[fanout]['i']
            u = ("t%d" % fanouts.index(fanout)) if tap is None else "t_%d_%d" % (tap, fanouts.index(fanout))
            if u == target_mux:
                u = 't'
            lattice.add_node(u, coords = [x, y], mux = True, size = all_sizes[fanout], potential_target = True)
            found = False
            for col in cols:
                if abs(x - col) <= col_thr:
                    cols[col].append((y, u))
                    found = True
                    break
            if not found:
                cols.update({x : [(y, u)]})

        sorted_cols = sorted(cols)
        fuse, fuse_index = sorted((val, ind) for (ind, val) in enumerate(sorted(cols)))[len(cols) / 2]
        for i, col in enumerate(sorted(cols)):
            junction = "%sjunction_%d" % (("tp_%d_" % tap) if tap is not None else '', i)
            lattice.add_node(junction, coords = [col, wire_y], mux = False)
            cols[col].append((wire_y, junction))
            if (i == fuse_index):
                lattice.add_node("wire_t" if tap is None else "tap_%d_t" % tap, coords = [wire_x, wire_y], mux = False)
                lattice.add_edge("wire_t" if tap is None else "tap_%d_t" % tap, junction, layer = "Mx")
            if i != 0:
                lattice.add_edge(junction, junction.rsplit('_', 1)[0] + "_%d" % (i - 1), layer = "Mx")
            cols[col].sort()
            for j, node in enumerate(cols[col][1:], 1):
                #By construction there must be at least 2 nodes in the column, one being the junction.
                u = node[1]
                v = cols[col][j - 1][1]
                lattice.add_edge(u, v, layer = "Mx")

        return lattice
    #------------------------------------------------------------------------#

    tile_width = max(get_metal_dimensions()[0], get_tile_dimensions(G)[0]) / FP
    tile_height = max(get_metal_dimensions()[1], get_tile_dimensions(G)[1]) / GP

    pins, all_sizes = stack_muxes(G, get_pins = True)

    L = int(wire[1:])

    net = nx.Graph()
    node_mapping = {}

    if wire[0] == 'V':
        source = source.split("_tap")[0]
        y_cur = pins[source + "_tap_0"]['o'][1]
        net.add_node("wire_s", coords = [-1, y_cur])
        net.add_edge("wire_s", "tap_0_s", layer = "My")
        min_x = float("inf")
        max_x = -1 * float("inf")
        for tap in range(0, tap_M):
            L_local = tap_phi if tap > 0 else L - (tap_M * tap_phi - 1)
            y_cur += L_local * tile_height
            net.add_node("tap_%d_s" % tap, coords = [0, y_cur])
            net.add_edge("tap_%s_s" % tap, "tap_%d_t" % tap, layer = "Mx" if tap != tap_M - 1 else "My")
            fanout = []
            if not SEPARATE_TAPS:
                all_fanout = G[source + "_tap_0"]
                all_fanout = [child for child in G[source + "_tap_0"] if not "potential_edge" in child]
                if not FIRST_CLIQUE_ITER:
                    real_fanout = [c for c in all_fanout\
                                   if (not G[source + "_tap_0"][c].get("delay_operating_point_only", False))]# and\
                                       #(G.node[c]["node_type"] in ("h_track", "v_track"))]
                    if real_fanout:
                        all_fanout = real_fanout
                local_fanout = set()
                for child in all_fanout:
                    if "potential_edge" in child:
                        continue
                    if not child in pins:
                        if G.node[child]["node_type"] == "h_track":
                            child_ble = int(child.split("_H")[0].split("ble_")[1])
                        elif G.node[child]["node_type"] == "v_track":
                            child_ble = int(child.split("_V")[0].split("ble_")[1])
                        else:
                            child_ble = NEUTRAL_BLE
                        child_ble_offset = child_ble - NEUTRAL_BLE
                        child_pins = pins[child.replace("ble_%d" % child_ble, "ble_%d" % NEUTRAL_BLE)]
                        child_i = (child_pins['i'][0], child_pins['i'][1] + child_ble_offset * lut_height)
                        child_o = (child_pins['o'][0], child_pins['o'][1] + child_ble_offset * lut_height)
                        pins.update({child : {'i' : child_i, 'o' : child_o}})
                        all_sizes.update({child : all_sizes[child.replace("ble_%d" % child_ble, "ble_%d" % NEUTRAL_BLE)]})

                    for e in G[source + "_tap_0"][child]:
                        e_tap = G[source + "_tap_0"][child][e]["tap"]
                        if e_tap == -1:
                            e_tap = tap_M - 1
                        if e_tap == tap:
                            local_fanout.add(child)

                #print tap, source + "_tap_0", sorted(local_fanout)
                #raw_input()

                fanout = list(local_fanout)
                fanout.sort(key = lambda u : pins[u]['i'][0] * FP - abs(pins[u]['i'][1]) * GP)
                #NOTE: x-coordinates of all multiplexers are negative.
            else:
                try:
                    fanout = [u for u in pins if u in G[source + ("_tap_%d" % tap)]]
                    fanout.sort(key = lambda u : (pins[u]['i'][1], pins[u]['i'][0]))
                except:
                    pass
            target_mux = None
            if tap == tap_M - 1:
                sb_muxes = [f for f in fanout if not "cb_out" in f]
                target_mux = sb_muxes[len(sb_muxes) / 2]
                target_mux = "t_%d_%d" % (tap, fanout.index(target_mux))
            avg_load_y = sum([pins[u]['i'][1] for u in fanout]) / len(fanout)
            lattice = connect_targets(fanout, pins, all_sizes, 0, avg_load_y, tap = tap)
            for node in lattice:
                if node == "wire_t":
                    continue
                net.add_node(node)
                net.node[node].update(lattice.node[node])
                net.node[node]["coords"][1] += y_cur
            for u, v, attrs in lattice.edges(data = True):
                if u == "wire_t":
                    u = "tap_%d_t" % tap
                if v == "wire_t":
                    v = "tap_%d_t" % tap
                net.add_edge(u, v)
                net[u][v].update(attrs)

            local_fanout_x = 0
            for t, target in enumerate(fanout):
                t_x, t_y = pins[target]['i']
                local_fanout_x += t_x
                if t_x < min_x:
                    min_x = t_x
                if t_x > max_x:
                    max_x = t_x
                u = "t_%d_%d" % (tap, t)
                if u == target_mux:
                    u = 't'
                node_mapping.update({target : u})
                if tap == tap_M - 1:
                    net.node[u]["potential_target"] = True
            if fanout:
                local_fanout_x /= float(len(fanout))
                net.node["tap_%d_t" % tap]["coords"][0] = local_fanout_x
            if tap != tap_M - 1:
                net.add_edge("tap_%d_s" % tap, "tap_%d_s" % (tap + 1), layer = "My")

            #Determine the states of muxes.
            total_cols = 0
            for u in fanout:
                t = node_mapping[u]
                rcnt, ccnt = get_mux_dimensions(net.node[t]["size"], lut_height)
                total_cols += ccnt
            partial_prob = 1.0 / (float(total_cols) / len(fanout))
            partial_fanout = int(math.ceil(partial_prob * len(fanout)))
    
            partial_space = len(fanout) / partial_fanout
            for i, f in enumerate(fanout):
                m = node_mapping[f]
                if i % partial_space == 0:
                    net.node[m]["state"] = "partial"
                else:
                    net.node[m]["state"] = "off"

            if tap == tap_M - 1:
                cb_muxes = [u for u in fanout if "cb_out" in u]
                cb_mux_on = cb_muxes[len(cb_muxes) / 2]
                net.node[node_mapping[cb_mux_on]]["state"] = "on"

        net.add_node('s', coords = pins[source + "_tap_0"]['o'], mux = True, size = all_sizes[source + "_tap_0"])
        net.add_edge('s', "wire_s", layer = "My")

        net.node['s']["state"] = "on"
        net.node['t']["state"] = "on"

        #Offset the spine to simulate an average case of horizontal distance between the mux and the wire.
        #TODO: We can always move the wire that we know is intended to replace a feedback closer to the optimal position.
        lut_x = 2 ** 4 * 10 / 2.0
        lut_y = lut_height / (2 ** (K - 4) / 2.0)
        avg_tile_x = 0.5 * (min_x + lut_width)
        avg_load_x = 0.5 * (min_x + max_x)
        avg_cb_x = 0.5 * (lut_x + min([pins[b]['i'][0] for b in pins if "cb_out" in b]))
        wire_s_x = avg_tile_x

        for u in net:
            if "tap" in u and u.endswith("_s"):
                net.node[u]["coords"][0] += wire_s_x
        net.node["wire_s"]["coords"][0] = wire_s_x

        #for node in sorted(node_mapping, key = lambda n : node_mapping[n]):
        #    print node, node_mapping[node]
        #print "..............................................................."
        #for u, v in sorted(net.edges()):
        #    print u, net.node[u]["coords"], "->", v, net.node[v]["coords"], net[u][v]["layer"]
        #exit(0)
    else:
        wire_s_x = pins[source]['o'][0]
        wire_s_y = 0.5 * lut_height
        all_fanout = [c for c in G[source] if not "potential_edge" in c]
        if not FIRST_CLIQUE_ITER:
            real_fanout = [c for c in all_fanout\
                           if (not G[source][c].get("delay_operating_point_only", False))]# and\
                               #(G.node[c]["node_type"] in ("h_track", "v_track"))]
            if real_fanout:
                all_fanout = real_fanout
        fanout = set()
        for child in all_fanout:
            if not child in pins:
                if G.node[child]["node_type"] == "h_track":
                    child_ble = int(child.split("_H")[0].split("ble_")[1])
                elif G.node[child]["node_type"] == "v_track":
                    child_ble = int(child.split("_V")[0].split("ble_")[1])
                else:
                    child_ble = NEUTRAL_BLE
                child_ble_offset = child_ble - NEUTRAL_BLE
                child_pins = pins[child.replace("ble_%d" % child_ble, "ble_%d" % NEUTRAL_BLE)]
                child_i = (child_pins['i'][0], child_pins['i'][1] + child_ble_offset * lut_height)
                child_o = (child_pins['o'][0], child_pins['o'][1] + child_ble_offset * lut_height)
                pins.update({child : {'i' : child_i, 'o' : child_o}})
                all_sizes.update({child : all_sizes[child.replace("ble_%d" % child_ble, "ble_%d" % NEUTRAL_BLE)]})

            fanout.add(child)
        fanout = list(fanout)
        fanout.sort(key = lambda u : abs(pins[u]['i'][0] - wire_s_x) * FP + abs(pins[u]['i'][1] - wire_s_y) * GP)
        avg_load_x = sum([pins[u]['i'][0] for u in fanout]) / len(fanout)
        avg_load_y = sum([pins[u]['i'][1] for u in fanout]) / len(fanout)
        
        #print source, '|', sorted(fanout)
        #raw_input()

        #print all_fanout
        #for c in all_fanout:
        #    print G[source][c]
        #print fanout
        #exit(0)

        sb_muxes = [f for f in fanout if not "cb_out" in f]
        target_mux = sb_muxes[len(sb_muxes) / 2]
        target_mux = "t%d" % fanout.index(target_mux)

        net = connect_targets(fanout, pins, all_sizes, avg_load_x, avg_load_y)
        for u in net:
            net.node[u]["coords"] = (net.node[u]["coords"][0] + L * tile_width, net.node[u]["coords"][1])
        net.add_node('s', coords = pins[source]['o'], mux = True, size = all_sizes[source])
        net.add_node("wire_s", coords = (wire_s_x, wire_s_y))
        for t, target in enumerate(fanout):
            t_x, t_y = pins[target]['i']
            u = "t%d" % t
            if u == target_mux:
                u = 't'
            node_mapping.update({target : u})

        net.add_edge("s", "wire_s", layer = "My")
        net.add_edge("wire_s", "wire_t", layer = "My")

        #for node in sorted(node_mapping, key = lambda n : node_mapping[n]):
        #    print node, node_mapping[node]
        #print "..............................................................."
        #for u, v in sorted(net.edges()):
        #    print u, net.node[u]["coords"], "->", v, net.node[v]["coords"], net[u][v]["layer"]
        #exit(0)

        #Determine the states of muxes.
        total_cols = 0
        for u in fanout:
            t = node_mapping[u]
            rcnt, ccnt = get_mux_dimensions(net.node[t]["size"], lut_height)
            total_cols += ccnt
        partial_prob = 1.0 / (float(total_cols) / len(fanout))
        partial_fanout = int(math.ceil(partial_prob * len(fanout)))

        partial_space = len(fanout) / partial_fanout
        for i, f in enumerate(fanout):
            m = node_mapping[f]
            if i % partial_space == 0:
                net.node[m]["state"] = "partial"
            else:
                net.node[m]["state"] = "off"

        net.node['s']["state"] = "on"
        net.node['t']["state"] = "on"
        cb_muxes = [u for u in fanout if "cb_out" in u]
        cb_mux_on = cb_muxes[len(cb_muxes) / 2]
        net.node[node_mapping[cb_mux_on]]["state"] = "on"
   
    #coords = []
    #for f in fanout:
    #    coords.append((pins[f]['i']))
    #for c in sorted(coords):
    #    print c
    #print "center", (wire_s_x, wire_s_y)

    if get_cb_delay:
        cb_x = 0.5 * (pins[cb_mux_on]['i'][0] + pins[cb_mux_on]['o'][0])
        wire_x = net.node['t']["coords"][0] - (0 if wire[0] == 'V' else L * tile_width) 
        #NOTE: We already counted the return delay from the wire offset to the target SB
        #in the delay of the wire itself so there is no need to count it again in the CB delay.
        cb_size = all_sizes[cb_mux_on]

        wd = os.getcwd()
        dump_filename = "cb_meas_%s_%s_%s.dump" % (args.arc_name, wd.rsplit('/', 1)[1], str(time.time()))
        os.chdir("../wire_delays")
        sim_dir = "cb_sim_%s_%s_%s/" % (args.arc_name, wd.rsplit('/', 1)[1], str(time.time()))
        os.system("mkdir %s" % sim_dir)
        os.system("python -u local_wires.py --K %d --N %d --tech % s --density %f --meas_cb %s --wd %s > %s"\
                  % (K, N, args.tech, density, "\"%f %f %d\"" % (wire_x, cb_x, cb_size), sim_dir, dump_filename))
        with open(dump_filename, "r") as inf:
            lines = inf.readlines()
            td = float(lines[-1].strip())
        os.system("rm -f %s" % dump_filename)
        os.system("rm -rf %s" % sim_dir)
        os.chdir(wd)
        return td

    return net
##########################################################################

##########################################################################
def look_up_load_delay(G, mux):
    """Extracts the delay change from the precomputed load-delay model.

    Parameters
    ----------
    G : nx.MultiDiGraph
        The routing-resource graph.
    mux : str
        Multiplexer whose load is being assessed.

    Returns
    -------
    float
        The predicted delay.
    float
        Loading length.
    """

    wire = mux.split('_')[2]
    net = get_netlist(G, wire, mux)

    load_length = 0.0
    for u, v in sorted(net.edges()):
        if "wire" in u and "wire" in v:
            continue
        u_coords = net.node[u]["coords"]
        v_coords = net.node[v]["coords"]

        L = abs(u_coords[1] - v_coords[1]) * GP
        L += abs(u_coords[0] - v_coords[0]) * FP
        
        #NOTE: We sum up all Mx wires, as this will approximate also the effect
        #of moving the multiplexer away from the LUT driver, which is important,
        #although not modeled separately at the moment. In the future, we may want
        #to be more precise.
        if net[u][v]["layer"] == "Mx":
            load_length += L
    
    #print load_length

    #NOTE: The polynomial approximation might be off and it may return negative values,
    #which we of course want to cap.

    return max(0.0, load_model[wire].evaluate(load_length)), load_length
##########################################################################

##########################################################################
def optimize_pattern_layout(G, criticalities = None):
    """Optimizes the switch-pattern using simulated annealing.

    Parameters
    ----------
    G : nx.MultiDiGraph
        The routing-resource graph.
    criticalities : Optional[Dict[str, float]], default = None
        A dictionary of wire-type criticalities, extracted from some VPR run.
        If not specified, all wire types will be assigned a criticality of 1.

    Returns
    -------
    List[str]
        The final switch-block multiplexer order.
    """

    #------------------------------------------------------------------------#
    def evaluate_move(prev_td, prev_wl, new_td, new_wl, T):
        """Evaluates the move

        Parameters
        ----------
        prev_td : float
            Previous timing cost.
        prev_wl : float
            Previous wirelength cost.
        new_td : float
            New timing cost.
        new_wl : float
            New wirelength cost.
        T : float
            Temperature.
        
        Returns
        -------
        bool
            Indicator of acceptance.
        """

        td_vs_wl = 0.5

        delta_td = (new_td - prev_td) / float(prev_td)
        delta_wl = (new_wl - prev_wl) / float(prev_wl)
        delta_c = td_vs_wl * delta_td + (1.0 - td_vs_wl) * delta_wl

        if delta_c < 0:
            return True

        energy_threshold = random.random()
        energy = math.exp(-1 * delta_c / T)

        return energy > energy_threshold
    #------------------------------------------------------------------------#

    #------------------------------------------------------------------------#
    def swap_pair(a, b, l):
        """Swaps two pairs in a list.

        Parameters
        ----------
        a : int
            First index.
        b : int
            Second index.
        l : List[str]
            List.
    
        Returns
        -------
        List[str]
            List with the swap committed.
        """

        cp = l[a]
        l[a] = l[b]
        l[b] = cp

        return l
    #------------------------------------------------------------------------#

    global SB_MUX_ORDER

    #Make a copy of the routing graph without the potential edges, so as to
    #make the optimization faster.
    G_real = copy.deepcopy(G)
    rm_list = [u for u in G if "potential_edge" in u]
    for u in rm_list:
        G_real.remove_node(u)

    cb_muxes, sb_muxes = export_mux_sizes(G_real)
    sb_muxes = list(sb_muxes.keys())

    random.seed(19225)
    random.shuffle(sb_muxes)
    SB_MUX_ORDER = sb_muxes
    sb_mux_indices = [i for i in range(0, len(sb_muxes))]
    prev_td, prev_wl = evaluate_pattern_layout_cost(G_real, criticalities = criticalities)

    #SA parameters:
    init_t = 1.0
    t_red = 0.8
    outer_iters = 30
    inner_iters = 100

    T = init_t
    move = -1
    for tcnt in range(0, outer_iters):
        for mcnt in range(0, inner_iters):
            move += 1
            swap_a, swap_b = random.sample(sb_mux_indices, 2)
            print move, T, "swapping", sb_muxes[swap_a], sb_muxes[swap_b]
            sb_muxes = swap_pair(swap_a, swap_b, sb_muxes)
            SB_MUX_ORDER = sb_muxes
            new_td, new_wl = evaluate_pattern_layout_cost(G_real, criticalities = criticalities)
            accept = evaluate_move(prev_td, prev_wl, new_td, new_wl, T)
            if accept:
                prev_td = new_td
                prev_wl = new_wl
                print "accepted", new_td, new_wl
            else:
                sb_muxes = swap_pair(swap_a, swap_b, sb_muxes)
                SB_MUX_ORDER = sb_muxes
                print "rejected"
        T *= t_red

    with open("sb_mux.order", "w") as outf:
        outf.write(str(SB_MUX_ORDER))
##########################################################################

##########################################################################
def evaluate_pattern_layout_cost(G, criticalities = None, do_spice = False):
    """Returns the delay and wirelength cost of the pattern's current layout.

    Parameters
    ----------
    G : nx.MultiDiGraph
        The routing-resource graph.
    criticalities : Optional[Dict[str, float]], default = None
        A dictionary of wire-type criticalities, extracted from some VPR run.
        If not specified, all wire types will be assigned a criticality of 1.
    do_spice : Optional[bool], default = False
        Specifies that the delays should be measured again.

    Returns
    -------
    float
        Timing cost.
    float
        Wirelength cost.
    """

    #TODO: See if we want to do full measurements at all, or operate merely
    #on the variable parts of the delay, which may bring more fidelity.
    #If we do want, then we should save the delays of the first architecture
    #in the clique search, which has the absolute minimum loading of wires,
    #then add the variable part that is being looked up.


    #Cost function parameters:
    crit_exp = 8.0

    if criticalities is None:
        criticalities = {}

    #TODO: Separate the wire types into individual wires per BLE as their delays may differ!

    total_wl = 0.0
    total_td = 0.0

    cb_muxes, sb_muxes = export_mux_sizes(G)
    for mux in sb_muxes:
        wire = mux.split('_')[2]
        crit = criticalities.get(wire, 1.0)
        delay, wirelength = look_up_load_delay(G, mux)
        delay *= (crit ** crit_exp)
        total_wl += wirelength
        total_td += delay

    return total_td, total_wl
##########################################################################

##########################################################################
def conv_nx_to_spice(net, meas_lut_access = False):
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

    txt = ".TITLE GLOBAL_WIRE_MEAS\n\n"
    txt += ".LIB %s\n" % (spice_model_path % (int(tech_node), int(tech_node)))
    txt += ".TRAN 0.1p 16n\n.OPTIONS BRIEF=1\n\n"

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
    txt += add_wire_subckt(N * 2 ** (K - 4), my_wire = True)

    txt += ".SUBCKT buf n_in n_out vdd STRENGTH0=D0 STRENGTH1=D1\n"
    txt += "MN1 n_mid n_in gnd gnd nmos L=GL nfin=STRENGTH0\n"
    txt += "MP1 n_mid n_in vdd vdd pmos L=GL nfin=STRENGTH0\n"
    txt += "MN2 n_out_pre_via n_mid gnd gnd nmos L=GL nfin=STRENGTH1\n"
    txt += "MP2 n_out_pre_via n_mid vdd vdd pmos L=GL nfin=STRENGTH1\n"
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

    mux_models = set()
    for u, attrs in net.nodes(data = True):
        if attrs.get("mux", False):
            size = attrs["size"]
            rcnt, ccnt = get_mux_dimensions(size, lut_height)
            model = add_mux_subckt(rcnt, ccnt, "on")
            model += add_mux_subckt(rcnt, ccnt, "off")
            model += add_mux_subckt(rcnt, ccnt, "partial")
            mux_models.add(model)
            net.node[u]["model"] = "mux_%d_%d_%s" % (rcnt, ccnt, attrs["state"])
    for model in sorted(mux_models):
        txt += model

    txt += "Vps vdd gnd supply_v\n"
    txt += "Vin n_in gnd PULSE (0 supply_v 0 0 0 2n 4n)\n\n"
    if meas_lut_access:
        txt += "Xtg_drv_mux_2_1 n_in n_in_mux vdd tg_on\n"
        txt += "Xbuf_drv n_in_mux n_wire_s vdd buf STRENGTH0=%d STRENGTH1=%d\n\n" % (local_driver[0], local_driver[1])
        txt += "Xbuf_local n_in_mux n_local_drv vdd buf STRENGTH0=%d STRENGTH1=%d\n\n" % (local_driver[0], local_driver[1])
    else:
        txt += "Xdriver_mux n_in n_in_drv vdd %s\n" % net.node['s']["model"]
        txt += "Xbuf_drv n_in_drv n_s vdd buf\n\n"

        try:
            u_coords = net.node["wire_s"]["coords"]
            v_coords = net.node["wire_t"]["coords"]
            L = abs(u_coords[1] - v_coords[1]) * GP
            L += abs(u_coords[0] - v_coords[0]) * FP
            txt += "Xwire_spine n_wire_s n_wire_t_pre_via my_wire L=%.2f\n" % L
            txt += "Rvia_spine_out n_wire_t_pre_via n_wire_t R=Rvia\n\n"
        except:
            pass

    ecnt = -1
    mux_cnt = -1
    for u, v in sorted(net.edges()):
        if "wire" in u and "wire" in v and not meas_lut_access:
            continue
        ecnt += 1
        u_coords = net.node[u]["coords"]
        v_coords = net.node[v]["coords"]

        L = abs(u_coords[1] - v_coords[1]) * GP
        L += abs(u_coords[0] - v_coords[0]) * FP

        model = "my_wire" if net[u][v]["layer"] == "My" else "wire"
        n_in = str(u)
        n_out = str(v)
        txt += "Xwire%d n_%s n_%s %s L=%.2f\n" % (ecnt, n_in, n_out, model, L)

        mux_node = None
        if net.node[v].get("mux", False) and v != 's':
            mux_node = v
        elif net.node[u].get("mux", False) and u != 's':
            mux_node = u
        if mux_node is not None:
            mux_cnt += 1
            n_in = str(mux_node)
            n_out = "n_dummy_mux_out%d" % mux_cnt
            txt += "Xmux%d n_%s %s vdd %s\n\n" % (mux_cnt, n_in, n_out, net.node[mux_node]["model"])

    template = ".MEASURE tfall%s TRIG V(n_%s) VAL='supply_v/2' FALL=2\n"
    template += "+                  TARG V(n_%s) VAL supply_v/2 FALL=2\n\n"
    template += ".MEASURE trise%s TRIG V(n_%s) VAL='supply_v/2' RISE=2\n"
    template += "+                  TARG V(n_%s) VAL supply_v/2 RISE=2\n\n"

    txt += template % ('', "in", 't', '', "in", 't')

    if meas_lut_access:
        txt += template % ("_ble_mux", "in", "in_mux", "_ble_mux", "in", "in_mux")
    else:
        for u in sorted(net):
            if "tap" in u and not u.endswith("_t"):
                tap = int(u.split('_')[1])
                if tap == 0:
                    txt += template % ("_tap_0", "in", "tap_0_s", "_tap_0", "in", "tap_0_s")
                else:
                    suffix = '_' + u.rsplit('_', 1)[0]
                    target = u
                    if tap == tap_M - 1:
                        target = 't'
                    txt += template % (suffix, "tap_%d_s" % (tap - 1), target, suffix, "tap_%d_s" % (tap - 1), target)
    
    txt += ".END"

    return txt
##########################################################################

##########################################################################
def measure(G, wire, get_cb_delay = False, meas_lut_access = False):
    """Calls HSPICE to obtain the delay of the wire.
    
    Parameters
    ----------
    G : nx.MultiDiGraph
        The routing-resource graph.
    wire : str
        Wire type.
    get_cb_delay : Optional[bool], default = False
        Determines the position of the wire and the connection block and then calls
        >>meas_local_wire.py<< to obtain the delay from the wire to a LUT input pin.

    Returns
    -------
    float
        Delay.
    """

    #------------------------------------------------------------------------#
    def run():
        """Runs HSPICE and parses the delay."""
 
        with open(netlist_filename, "w") as outf:
           outf.write(conv_nx_to_spice(net, meas_lut_access = meas_lut_access))
       
        hspice_call = os.environ["HSPICE"] + " %s > %s" % (netlist_filename, hspice_dump)
        os.system(hspice_call)
       
        scale_dict = {'f' : 1e-15, 'p' : 1e-12, 'n' : 1e-9}
       
        with open(hspice_dump, "r") as inf:
            lines = inf.readlines()
    
        #os.system("rm " + hspice_dump)
       
        td_dict = {}
      
        get_td = lambda l :  round(float(l.split()[1][:-1]), 1) * scale_dict[l.split()[1][-1]]
        get_tap = lambda l : wire + '_' + l.split('=', 1)[0].split('_', 1)[1]
        for line in lines:
            if "tfall=" in line:
                tfall = get_td(line) 
            elif "trise=" in line:
                trise = get_td(line)
            elif meas_lut_access:
                if "tfall_ble_mux" in line or "trise_ble_mux" in line:
                    td = get_td(line)
                    if td < 0:
                        print "Negative time!"
                        raise ValueError
                    try:
                        td_dict["ble_mux"] = 0.5 * (td_dict["ble_mux"] + td)
                    except:
                        td_dict.update({"ble_mux" : td})
            elif wire[0] == 'V':
                if "tfall_tap" in line or "trise_tap" in line:
                    tap = get_tap(line)
                    td = get_td(line)
                    if td < 0:
                        print "Negative time!"
                        raise ValueError
                    try:
                        td_dict[tap] = 0.5 * (td_dict[tap] + td)
                    except:
                        td_dict.update({tap : td})
                
        if trise < 0 or tfall < 0:
            print "Negative time!"
            raise ValueError
     
        if wire[0] == 'V':
            td_dict.update({"whole" : 0.5 * (trise + tfall)})

        if meas_lut_access:
            td_dict.update({"lut_access" : 0.5 * (trise + tfall) - td_dict["ble_mux"]})
            return td_dict
        if wire[0] == 'V':
            return td_dict
      
        return 0.5 * (trise + tfall)
    #------------------------------------------------------------------------#

    netlist_filename = "sim_global_%s_%s.sp" % (args.arc_name, wire)
    hspice_dump = "hspice_%s_%s.dump" % (args.arc_name, wire)

    if meas_lut_access:
        net = meas_lut_access_delay(G)
        return run()
    else:
        pins, all_sizes = stack_muxes(G, get_pins = True)
        source_dict = {}
        for mux in pins:
            if wire in mux and mux.startswith("ble_%d_" % NEUTRAL_BLE):

                if ROBUSTNESS_LEVEL == 0:
                    source = mux
                    if get_cb_delay:
                        return get_netlist(G, wire, source, get_cb_delay = True)
                    net = get_netlist(G, wire, source)
                    return run()                   

                key = mux.split("_tap")[0]
                offset = pins[mux]['o'][0 if wire[0] == 'V' else 1]
                deg = 0
                for fanout in G:
                    if fanout.startswith(key):
                        deg += G.in_degree(fanout) + G.out_degree(fanout)
                source_dict.update({key : {"mux" : mux, "deg" : deg, "offset" : offset}})

        sorted_keys = sorted(source_dict, key = lambda s : source_dict[s]["deg"]\
                             * abs(source_dict[s]["offset"]))

        if ROBUSTNESS_LEVEL == 1 or get_cb_delay:
            #NOTE: Connection-block delays are very robust to changing the multiplexer as they usually
            #assume only one or two columns, immediately next to the crossbar. Hence, the x-offset is
            #less varialbe. Also, the load is within the cluster itself. If there is any variation in
            #multiplexer sizes, that is more of an artifact of parametrized architecture generation.
            #Median fanin should be a good representative in this case.
            source = source_dict[sorted_keys[len(source_dict) / 2]]["mux"]

            if get_cb_delay:
                return get_netlist(G, wire, source, get_cb_delay = True)
            net = get_netlist(G, wire, source)
            return run()                   

        td_dicts = []
        for source_key in sorted_keys:
            source = source_dict[source_key]["mux"]
            net = get_netlist(G, wire, source)
            td_dicts.append(run())
           
            if ROBUSTNESS_LEVEL == 3: 
                potential_targets = [u for u, attrs in net.nodes(data = True) if attrs.get("potential_target", False)]
                for i, u in enumerate(potential_targets):
                    relabeling_dict = {}
                    if u == 't':
                        continue
                    relabeling_dict.update({'t' : "prev_t_%d" % i})
                    relabeling_dict.update({u : 't'})
                    net = nx.relabel_nodes(net, relabeling_dict)
                    td_dicts.append(run())
            
        if (wire[0] == 'H' and not meas_lut_access) or get_cb_delay:
            return sum(td_dicts) / len(td_dicts)

        for v in td_dicts[0]:
            for td_dict in td_dicts[1:]:
                td_dicts[0][v] += td_dict[v]
            td_dicts[0][v] /= len(td_dicts)
    
        return td_dicts[0]
##########################################################################

##########################################################################
class TimingModelEntry(object):
    """Defines a timing model entry, used for switch-pattern optimization.

    Parameters
    ----------
    measurements : Dict[int, float]
        Performed spice measurements
    
    Methods
    -------
    float evaluate(l : int)
        Evaluates the model
    """

    #------------------------------------------------------------------------#
    def __init__(self, measurements):
        """Constructor of the timing model entry class."""

        xs = list(sorted(measurements))
        ys = [measurements[x] for x in xs]

        self.polydeg = 3

        #Because we need high resolution at the short lengths, to be able to
        #drive optimization, we will use linear optimization here. Otherwise,
        #the best that we could do would be to clip the np polynomial to zero.
        self.first_point = min(xs)
        self.t_per_nm = measurements[self.first_point] / self.first_point
        interpolated = np.polyfit(xs, ys, self.polydeg)
        self.curve = np.poly1d(interpolated)
    #------------------------------------------------------------------------#

    #------------------------------------------------------------------------#
    def evaluate(self, l):
        """Evaluates the model.

        Parameters
        ----------
        l : int
            Load wire length in nm.
        
        Returns
        -------
        float
            The modeled delay.
        """

        if l < self.first_point:
            return l * self.t_per_nm

        return self.curve(l)
    #------------------------------------------------------------------------#
##########################################################################

##########################################################################
def create_load_model():
    """Creates a load model for the given wire, depending on the length of
    the Mx wires at its output, which is the predominant factor influencing the
    observed delay (much more than the muxes themselves).

    Parameters
    ----------
    None
    
    Returns
    -------
    numpy.poly1d
        A polynomial that interpolates the few measurements that are actually performed.
    """

    #------------------------------------------------------------------------#
    def measure_load_impact(global_wire_length, sb_wire_length):
        """Measures the impact of the increased load, for the given global and intra-sb wire lengths.

        Parameters
        ----------
        global_wire_length : float
            Length of the global wire, that is being loaded, in um.
        sb_wire_length : float
            Length of the intra-sb wire exercising the load, in um.
        """
 
        txt = ".TITLE LOAD_MODEL_MEAS\n\n"
        txt += ".LIB %s\n" % (spice_model_path % (int(tech_node), int(tech_node)))
        txt += ".TRAN 0.1p 16n\n.OPTIONS BRIEF=1\n\n"
    
        txt += ".PARAM Cw=%g\n" % Cw
        txt += ".PARAM Rw=%g\n" % Rw
        txt += ".PARAM Cwy=%g\n" % Cwy
        txt += ".PARAM Rwy=%g\n" % Rwy
        txt += ".PARAM Rvia=%g\n\n" % Rvia
        txt += ".PARAM D0=%d\n" % D0
        txt += ".PARAM D1=%d\n\n" % D1
        txt += ".PARAM GL=%dn\n" % GL
        txt += ".PARAM supply_v=%.2f\n\n" % vdd
    
        #........................................................................#
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
        #........................................................................#
    
        txt += add_wire_subckt(1)
        txt += add_wire_subckt(N * 2 ** (K - 4), my_wire = True)
    
        txt += ".SUBCKT buf n_in n_out vdd STRENGTH0=D0 STRENGTH1=D1\n"
        txt += "MN1 n_mid n_in gnd gnd nmos L=GL nfin=STRENGTH0\n"
        txt += "MP1 n_mid n_in vdd vdd pmos L=GL nfin=STRENGTH0\n"
        txt += "MN2 n_out_pre_via n_mid gnd gnd nmos L=GL nfin=STRENGTH1\n"
        txt += "MP2 n_out_pre_via n_mid vdd vdd pmos L=GL nfin=STRENGTH1\n"
        txt += "Rout n_out_pre_via n_out R=Rvia\n"
        txt += ".ENDS\n\n"
    
        txt += "Vps vdd gnd supply_v\n"
        txt += "Vin n_in gnd PULSE (0 supply_v 0 0 0 2n 4n)\n\n"
        txt += "Xbuf_drv n_in n_s vdd buf\n\n"
    
        txt += "Xwire_global n_s n_wire_t_pre_via my_wire L=%.2f\n" % global_wire_length
        txt += "Rdescend_via n_wire_t_pre_via n_wire_t R=Rvia\n"
        txt += "Xwire_sb n_wire_t n_t wire L=%.2f\n" % sb_wire_length
    
        txt += ".MEASURE tfall TRIG V(n_s) VAL='supply_v/2' FALL=2\n"
        txt += "+                  TARG V(n_t) VAL supply_v/2 FALL=2\n\n"
        txt += ".MEASURE trise TRIG V(n_s) VAL='supply_v/2' RISE=2\n"
        txt += "+                  TARG V(n_t) VAL supply_v/2 RISE=2\n\n"
    
        txt += ".END"

    #........................................................................#
        def run(txt):
            """Runs HSPICE and parses the delay.
            
            Parameters
            ----------
            txt : str
                Spice netlist text.
            

            Returns
            -------
            float
                The measured delay.
            """

            with open("meas.sp", "w") as outf:
               outf.write(txt)
           
            hspice_call = os.environ["HSPICE"] + " meas.sp > meas.dump"
            os.system(hspice_call)
           
            scale_dict = {'a' : 1e-18, 'f' : 1e-15, 'p' : 1e-12, 'n' : 1e-9}
           
            with open("meas.dump", "r") as inf:
                lines = inf.readlines()
        
            os.system("rm meas.sp meas.dump")
           
            td_dict = {}
          
            get_td = lambda l : 0 if l.split()[1] == "0." else round(float(l.split()[1][:-1]), 1) * scale_dict[l.split()[1][-1]]
            for line in lines:
                if "tfall=" in line:
                    tfall = get_td(line) 
                elif "trise=" in line:
                    trise = get_td(line)
                    
            if trise < 0 or tfall < 0:
                print "Negative time!"
                raise ValueError
         
            return 0.5 * (trise + tfall)
        #........................................................................#

    
        return run(txt)
    #-----------------------------------------------------------------------#

    measure_count = 8
    max_reasonable_fanout = 16
    hpwl_coeff = 1.7709
    #Taken from VTR7, in turn from RISA, for a fanout of 16.

    active_w, active_h = get_tile_dimensions(G)
    active_h /= float(N)
    active_h *= MAX_BLE_SPAN
    max_length = (active_w + active_h) * hpwl_coeff
    #* min(max_reasonable_fanout, 2 * (sum([h[1] for h in H]) * sum([v[1] for v in V]) * N))

    high_lengths = [max_length / float(measure_count) * m for m in range(1, measure_count + 1)]
    low_lengths = [lut_height / float(2 ** (K - 4)) * m for m in range(1, measure_count + 1)]
    #NOTE: Impact of such short wires is not even measurable. We will use the wirelength in the cost function there.
    low_lengths = []

    tile_width = max(get_metal_dimensions()[0], get_tile_dimensions(G)[0])
    tile_height = max(get_metal_dimensions()[1], get_tile_dimensions(G)[1])

    model = {}

    for h in sorted(H):
        L = float(h[0]) * tile_width
        h_id = "H%d" % h[0]
        D0, D1 = H_drivers[h[0]]
        model.update({h_id : {}})
        unloaded_delay = measure_load_impact(L, 0.0)
        for l in (low_lengths + high_lengths):
            model[h_id].update({l : measure_load_impact(L, l) - unloaded_delay})

    for v in sorted(V):
        L = float(v[0]) * tile_height
        v_id = "V%d" % v[0]
        D0, D1 = V_drivers[v[0]]
        model.update({v_id : {}})
        unloaded_delay = measure_load_impact(L, 0.0)
        for l in (low_lengths + high_lengths):
            model[v_id].update({l : measure_load_impact(L, l) - unloaded_delay})

    polydeg = 3
    #Degree of the approximating polynomial
    for w in model:
        model[w] = TimingModelEntry(model[w])
        print w, model[w].evaluate(100), model[w].evaluate(500), model[w].evaluate(1000), model[w].evaluate(10000)

    return model
##########################################################################

##########################################################################
def spice_all_wires(G):
    """Returns the delay of all wires.

    Parameters
    ----------
    None

    Returns
    -------
    Dict[str, float]
        A dictionary of delays.
    """

    td_dict = {}
    global D0
    global D1

    for h in sorted(H):
        h_id = "H%d" % h[0]
        D0, D1 = H_drivers[h[0]]
        td_dict.update({h_id : measure(G, h_id)})
    for v in sorted(V):
        v_id = "V%d" % v[0]
        D0, D1 = V_drivers[v[0]]
        td_dict.update(measure(G, v_id))
        if not SEPARATE_TAPS:
            for tap in range(1, tap_M):
                try:
                    td_dict.pop(v_id + "_tap_%d" % tap)
                except:
                    pass
            td_dict[v_id + "_tap_0"] = td_dict["whole"]
            td_dict.pop("whole")

    VL = "V%d" % (max(1, K6N8_LUT4 / KN_LUT4))

    td_dict.update({"cb_h" : measure(G, h_id, get_cb_delay = True)})
    td_dict.update({"cb_v" : measure(G, v_id, get_cb_delay = True)})
    td_dict.update(measure(G, "H1", meas_lut_access = True))

    to_remove = []
    for f in os.listdir('.'):
        if not args.arc_name in f:
            continue
        if f.endswith(".sp") or f.endswith(".st0") or f.endswith(".ic0") or f.endswith(".mt0") or f.endswith(".dump"):
            to_remove.append(f)
    os.system("rm -f %s" % ' '.join(to_remove))
            
    return td_dict
##########################################################################

##########################################################################
def read_delays_from_arc(arc_filename):
    """Reads the switch delays from an architecture file.

    Parameters
    ----------
    arc_filename : str
        Name of the architecture file from which to read the delays.

    Returns
    -------
    Dict[str, float]
        A dictionary of delays.
    """

    get_td = lambda l : float(l.split()[1].split('"')[1])

    with open(arc_filename, "r") as inf:
        lines = inf.readlines()

    td_dict = {}
    for line in lines:
        if "<switch " in line:
            name = line.split("name=\"")[1].split('"')[0]
            td = float(line.split("Tdel=\"")[1].split('"')[0])
            td_dict.update({name : td})
        elif "delay_constant" in line:
            if "in_port=\"ble" in line and "out_port=\"clb" in line:
                td_dict.update({"lut_access" : get_td(line)})
            elif "in_port=\"ble" in line and "out_port=\"ble" in line:
                global local_feedback_delay
                local_feedback_delay = get_td(line)
            elif "in_port=\"lut.out" in line and "out_port=\"ble.out" in line:
                td_dict.update({"ble_mux" : get_td(line)})                

    return td_dict
##########################################################################

##########################################################################
def fill_in_template(G, cb_delay, td_sb, potential_edges = None):
    """Fills in the architecture xml template.

    Parameters
    ----------
    G : nx.MultiDiGraph
        RR-graph.
    cb_delay : float
        Delay of the connection block.
    td_sb : Dict[str, float]
        A dictionary of switch-block delays.
    potential_edges : Optional[List[str]], default = None
        List of potential edges obtained with the appropriate >>make_clique<< switches turned on.
        Optionally also a dictionary of their delays can be passed.

    Returns
    -------
    None
    """

    if INSERT_EMPTY_RING:
        template_filename = "arc_model/templates/minimal_template_independent_ff_empty_ring_ready.xml"
    else:
        template_filename = "arc_model/templates/minimal_template_independent_ff.xml"

    td_lut = tech.td_lut[tech_node][K]
    td_ff_su = tech.td_ff_su[tech_node]
    td_ff_clkq = tech.td_ff_clkq[tech_node]
    td_inpad = tech.td_inpad
    td_outpad = tech.td_outpad
    td_crossbar_in = 0
    #NOTE: All cluster-input-to-LUT delay is attributed to the connection block.
    td_crossbar_fb = local_feedback_delay
    #FIXME: default_cb_td is set to cluster feedback delay. Fix this when it no longer holds.
    td_ble_mux_lut = td_sb.get("ble_mux", tech.td_ble_mux_lut)
    td_ble_mux_ff = td_sb.get("ble_mux", tech.td_ble_mux_ff)
    td_lut_access = td_sb["lut_access"]

    lut_del_matrix = "%g\n" % td_lut + (K - 1) * (6 * indent + "%g\n" % td_lut)

    O = 2
    replacement_dict = {\
                        "%%LAYOUT%%" : export_fpga_layout(),\
                        "%%SWITCHLIST%%" : export_fpga_switches(H, V, cb_delay, td_sb, potential_edges = potential_edges),\
                        "%%SEGMENTLIST%%" : export_fpga_segments(H, V, potential_edges = potential_edges),\
                        "%%N%%" : "%d" % N,\
                        "%%IO_CAPACITY%%" : "%d" % IO_CAPACITY,\
                        "%%TD_INPAD%%" : "%g" % td_inpad,\
                        "%%TD_OUTPAD%%" : "%g" % td_outpad,\
                        "%%I%%" : "%d" % cluster_inputs,\
                        "%%N*O%%" : "%d" % (N * O),\
                        "%%K%%" : "%d" % K,\
                        "%%K-1%%" : "%d" % (K - 1),\
                        "%%O%%" : "%d" % O,\
                        "%%TD_LUT%%" : lut_del_matrix[:-1],\
                        "%%TD_FF_SU%%" : "%g" % td_ff_su,\
                        "%%TD_FF_CLKQ%%" : "%g" % td_ff_clkq,\
                        "%%TD_BLE_MUX_LUT%%" : "%g" % td_ble_mux_lut,\
                        "%%TD_BLE_MUX_FF%%" : "%g" % td_ble_mux_ff,\
                        "%%N-1%%" : "%d" % (N - 1),\
                        "%%TD_CROSSBAR_IN%%" : "%g" % td_crossbar_in,\
                        "%%TD_CROSSBAR_FB%%" : "%g" % td_crossbar_fb,\
                        "%%TD_LUT_ACCESS%%" : "%g" % td_lut_access\
                       }                    

    with open(template_filename, "r") as inf:
       txt = inf.read()

    for key in replacement_dict:
        txt = txt.replace(key, replacement_dict[key])

    with open(args.arc_name, "w") as outf:
        outf.write(txt)
##########################################################################

##########################################################################
def get_cost_changes_per_switch(G, criticalities = None):
    """Returns changes in wirelength and delay of the pattern, due to
    addition of a switch that may potentially exist. It is assumed that
    the current placement of the SB-muxes is retained.

    Parameters
    ----------
    G : nx.MultiDiGraph
        The routing-resource graph.
    criticalities : Optional[Dict[str, float]], default = None
        A dictionary of wire-type criticalities, extracted from some VPR run.
        If not specified, all wire types will be assigned a criticality of 1.
    
    Returns
    -------
    Dict[str, Tuple[float]]
        A dictionary of (timing, wirelength) cost changes, per switch.
    """

    potential_switches = [u for u in G if u.startswith("potential_edge__ble_%d" % NEUTRAL_BLE)]

    G_real = copy.deepcopy(G)
    rm_list = [u for u in G if "potential_edge" in u]
    for u in rm_list:
        G_real.remove_node(u)

    print len(potential_switches)
    max_td = 1e-12
    max_wl = 1e-12
    changes = {}
    cur_td, cur_wl = evaluate_pattern_layout_cost(G, criticalities) 
    for switch in sorted(potential_switches):
        p = list(G.pred[switch])[0]
        c = list(G[switch])[0]
        G_real.add_edge(p, c, tap = -1)
        new_td, new_wl = evaluate_pattern_layout_cost(G_real, criticalities)
        delta_td = max(0.0, (new_td - cur_td) / float(cur_td))
        max_td = max(max_td, delta_td)
        delta_wl = max(0.0, (new_wl - cur_wl) / float(cur_wl))
        max_wl = max(max_wl, delta_wl)
        #Also store the intrinsic switch delay. Could be useful to avoid further tuning parameters.
        changes.update({lut_canonical_potential_edge(switch) : (delta_td, delta_wl, max(0.0, new_td - cur_td))})
        G_real.remove_edge(p, c)

    #Now normalize the costs.
    for e in changes:
        delta_td = changes[e][0] / max_td
        delta_wl = changes[e][1] / max_wl
        changes[e] = (delta_td, delta_wl, changes[e][2])

    return changes
##########################################################################

##########################################################################
def export_base_costs(G, filename, criticalities = None):
    """Exports the base costs of all potential edges.

    Parameters
    ----------
    G : nx.MultiDiGraph
        The routing-resource graph.
    filename : str
        Name of the file into which to write the costs.
    criticalities : Optional[Dict[str, float]], default = None
        A dictionary of wire-type criticalities, extracted from some VPR run.
        If not specified, all wire types will be assigned a criticality of 1.


    Returns
    -------
    Dict[str, float]
        A dictionary of potential edge delays as modeled in >>get_cost_changes_per_switch<<.
    """

    base_cost = 1.0
    reset_cost = 1e-6
    variable_cost = 0.2
    scaling_factor = 10
    avalanche_p = -1
    avalanche_h = -1
    avalanche_d = 0
    avalanche_iter = 16
    splitter_crit_exponent = 8.0
    splitter_target_crit_cost = 1e-12
    if os.path.exists("avalanche.conf"):
        with open("avalanche.conf", "r") as inf:
            lines = inf.readlines()
        for line in lines:
            if "base_cost" in line:
                base_cost = float(line.split()[1])
            elif "reset_cost" in line:
                reset_cost = float(line.split()[1])
            elif "variable_cost" in line:
                variable_cost = float(line.split()[1])
            elif "avalanche_p" in line:
                avalanche_p = float(line.split()[1])
            elif "avalanche_h" in line:
                avalanche_h = float(line.split()[1])
            elif "avalanche_d" in line:
                avalanche_d = float(line.split()[1])
            elif "scaling_factor" in line:
                scaling_factor = float(line.split()[1])
            elif "avalanche_iter" in line:
                avalanche_iter = float(line.split()[1])
            elif "splitter_crit_exponent" in line:
                splitter_crit_exponent = float(line.split()[1])
            elif "splitter_target_crit_cost" in line:
                splitter_target_crit_cost = float(line.split()[1])

    potential_edges = sorted([u for u in G if u.startswith("potential_edge")])
    if not potential_edges:
        with open(filename, "w") as outf:
            outf.write("0 0 0 0 0\n1 %f\n0 -%d" % (splitter_target_crit_cost, scaling_factor))
        return

    pins, all_sizes = stack_muxes(G, get_pins = True)

    predetermined_costs = set()
    costs = {}
    changes = get_cost_changes_per_switch(G, criticalities = criticalities)
    for e in potential_edges:
        if G.node[e]["cost"] >= 0:
            costs.update({lut_canonical_potential_edge(e) : G.node[e]["cost"]})
            predetermined_costs.add(lut_canonical_potential_edge(e))
            continue
        u = e.split("__")[1]
        v = e.split("__")[2]
        u_ble = '_'.join(u.split('_')[:2])
        u_pin = u.replace(u_ble, "ble_%d" % NEUTRAL_BLE)
        u_ble = int(u_ble.split('_')[-1])
        u_coords = pins[u_pin]['o']
        v_ble = '_'.join(v.split('_')[:2])
        v_pin = v.replace(v_ble, "ble_%d" % NEUTRAL_BLE)
        v_ble = int(v_ble.split('_')[-1])
        v_coords = pins[v_pin]['o']
        distance = abs(v_coords[0] - u_coords[0]) * FP\
                 + (abs(v_coords[1] - u_coords[1]) + abs(v_ble - u_ble) * lut_height) * GP
        costs.update({lut_canonical_potential_edge(e) : changes[lut_canonical_potential_edge(e)][1]})
        #NOTE: Timing cost will influence the delay of the edge.


    #NOTE: Enable the switch if the avalanche parameters are not to be recalculated at each outer iteration.
    if not FIRST_CLIQUE_ITER:
        with open("vpr_stdout.log", "r") as inf:
            lines = inf.readlines()
        for line in lines:
            if "avalanche_p =" in line:
                avalanche_p = float(line.split()[-1]) / (10.0 ** (-1 * scaling_factor))
            elif "avalanche_h =" in line:
                avalanche_h = float(line.split()[-1]) / (10.0 ** (-1 * scaling_factor))
            elif "avalanche_d =" in line:
                avalanche_d = float(line.split()[-1]) / (10.0 ** (-1 * scaling_factor))

    if GREEDY_SWITCH_SEARCH:
        avalanche_p = avalanche_h = avalanche_d = reset_cost = 0.0
        avalanche_iter = 1
 
    txt = "%f %f %f %d %f\n" % (avalanche_p, avalanche_h, avalanche_d, avalanche_iter, reset_cost)
    #Avalanche cost update: proportional, historical, and differential terms.
    #All quantities are scaled by 10^(-scaling_factor).
    #NOTE: If costs are -1, they will be dynamically determined on the VPR side.

    txt += "%f %f\n" % (splitter_crit_exponent, splitter_target_crit_cost)
    #Exponent to scale net criticality when multiplying the edge-splitter cost, target cost visible to a net at criticality of 0.99.

    txt += "%d -%d\n" % (len(costs), scaling_factor)
    #Number of potential .edges and the scaling factor.

    CHANX_START_INDEX = 4
    #VPR starts storing the segment base cost information from index 4 in the array (preceded by SRC, SINK, IPIN, OPIN).
    #norm_fact = float(max(costs.values()))
    #if norm_fact == 0:
    #    norm_fact = 1.0

    ecnt = min([seg_ids[c] for c in costs]) - 1
    #NOTE: VPR assumes that both the vertical and the horizontal channel contain all segments and duplicates the datastructure.
    #This is a waste of space, but likely has no real consequences.
    for e in sorted(costs, key = lambda k : seg_ids[k]):
        ecnt += 1
        cost = costs[e] if e in predetermined_costs\
               else base_cost + (costs[e] * variable_cost) 
        txt += "%d %f\n" % (CHANX_START_INDEX + ecnt, cost if not GREEDY_SWITCH_SEARCH else 0.0)

    with open(filename, "w") as outf:
        outf.write(txt[:-1])

    max_delay = 1.5e-12
    min_delay = 0.5e-12
    #timing_costs = {e : min_delay + (max_delay - min_delay) * changes[e][0] for e in changes}
    #Use intrinsic delay instead.
    timing_costs = {e : changes[e][2] for e in changes}

    return timing_costs
##########################################################################

##########################################################################
def generate_optimal_placement_delay_matrix(cb_delay, td_sb, filename):
    """Generates an optimal placement delay matrix, with the given
    channel composition, assuming that an arbitrary SB can be constructed,
    without further delay penalties. This matrix can be then loaded in VPR,
    to avoid the large runtime cost of computing a matrix in the presence
    of potential-edge cliques.

    Parameters
    ----------
    cb_delay : float
        Delay of the connection block.
    td_sb : Dict[str, float]
        A dictionary of switch-block delays.
    filename : str
        Name of the file in which the matrix is to be written.

    Returns
    -------
    None
    """

    #------------------------------------------------------------------------#
    def write_placement_delay_matrix(matrix, cb_delay, filename):
        """Writes out a placement delay matrix, one entry per line.
        >>inf<< is represented as -1. Other values are nonnegative.
    
        Parameters
        ----------
        matrix : List[int]
            Row-major list of delays in picoseconds.
        cb_delay : float
            Delay of the connection block.
        filename : str
            Name of the file in which the matrix is to be written.
    
        Returns
        -------
        None
        """
    
    
        transform = lambda td : (td + cb_delay) * 1e12
        #An optional further delay transform.
    
        txt = ""
        for d in matrix:
            txt += "%d\n" % (transform(d) if d != float('inf') else -1)
    
        with open(filename, "w") as outf:
            outf.write(txt[:-1])
    #------------------------------------------------------------------------#

    
    #------------------------------------------------------------------------#
    def add_tracks(graph, x, y):
        """Adds track instances to the periodic graph.

        Parameters
        ----------
        graph : nx.DiGraph
            Periodic graph to be extended.
        x : int
            x-coordinate.
        y : int
            y-coordinate.
        
        Returns
        -------
        None
        """

        S = "S_(%d,%d)" % (x, y)

        for h in H:
            L = h[0]
            ur = "H%dR_(%d,%d)_s" % (L, x, y)
            vr = "H%dR_(%d,%d)_t" % (L, x, y)
            graph.add_node(ur, coords = [x, y])
            graph.add_node(vr, coords = [x + L, y])
            graph.add_edge(ur, vr, td = td_sb["H%d" % L])
            graph.add_edge(S, ur, td = 0.0)
    
            ul = "H%dL_(%d,%d)_s" % (L, x, y)
            vl = "H%dL_(%d,%d)_t" % (L, x, y)
            graph.add_node(ul, coords = [x, y])
            graph.add_node(vl, coords = [x - L, y])
            graph.add_edge(ul, vl, td = td_sb["H%d" % L]) 
            graph.add_edge(S, ul, td = 0.0)
    
        for v in V:
            L = v[0]
            uu = "V%dU_(%d,%d)_s" % (L, x, y)
            vu = "V%dU_(%d,%d)_t" % (L, x, y)
            graph.add_node(uu, coords = [x, y])
            graph.add_node(vu, coords = [x, y + L])
            graph.add_edge(uu, vu, td = td_sb["V%d_tap_0" % L]) 
            graph.add_edge(S, uu, td = 0.0)
    
            ud = "V%dD_(%d,%d)_s" % (L, x, y)
            vd = "V%dD_(%d,%d)_t" % (L, x, y)
            graph.add_node(ud, coords = [x, y])
            graph.add_node(vd, coords = [x, y - L])
            graph.add_edge(ud, vd, td = td_sb["V%d_tap_0" % L]) 
            graph.add_edge(S, ud, td = 0.0)
    #------------------------------------------------------------------------#

    periodic_graph = nx.DiGraph()
    for x in range(0, grid_w):
        for y in range(0, grid_h):
            periodic_graph.add_node("S_(%d,%d)" % (x, y), coords = (x, y))
            periodic_graph.add_node("T_(%d,%d)" % (x, y), coords = (x, y))
            add_tracks(periodic_graph, x, y)
        
    coord_dict = {}                
    for u, attrs in periodic_graph.nodes(data = True):
        try:
            coord_dict[tuple(attrs["coords"])].append(u)
        except:
            coord_dict.update({tuple(attrs["coords"]) : [u]})

    for coords in coord_dict:
        x, y = coords
        T = "T_(%d,%d)" % (x, y)
        for u in coord_dict[coords]:
            if not u.endswith("_t"):
                continue
            for v in coord_dict[coords]:
                if not v.endswith("_s"):
                    continue
                periodic_graph.add_edge(u, v, td = 0)
            periodic_graph.add_edge(u, T, td = 0)

    S = "S_(0,0)"
    delays = nx.single_source_dijkstra_path_length(periodic_graph, S, weight = "td")

    matrix = []
    for y in range(0, grid_h):
        for x in range(0, grid_w):
            if (x, y) == (0, 0):
                matrix.append(0)
                continue
            T = "T_(%d,%d)" % (x, y)
            matrix.append(delays[T])

    write_placement_delay_matrix(matrix, cb_delay, filename)
##########################################################################

G, grid = generate_rr_graph()
print("Started padding LEN-1 wires.\n")
added_H1, added_V1, G, grid = pad_LEN1(G)
if added_H1 == 0:
    H.pop()
if added_V1 == 0:
    V.pop()

if PHYSICAL_SQUARE:
    tile_width = max(get_metal_dimensions()[0], get_tile_dimensions(G)[0])
    tile_height = max(get_metal_dimensions()[1], get_tile_dimensions(G)[1])
    multiplier = (float(tile_width) / tile_height) ** 0.5
    total_tiles = grid_w * grid_h
    print("Old grid: %d X % d" % (grid_w, grid_h))
    grid_w = int(math.ceil(grid_w / multiplier))
    if grid_w < ABS_MIN_WIDTH:
        grid_w = ABS_MIN_WIDTH
    grid_h = int(math.ceil(float(total_tiles) / grid_w))
    if grid_h < ABS_MIN_HEIGHT:
        grid_h = ABS_MIN_HEIGHT
        tentative_w = int(math.ceil(float(total_tiles) / grid_h))
        if tentative_w >= ABS_MIN_HEIGHT:
            grid_w = tentative_w
    print("New grid: %d X %d" % (grid_w, grid_h))
    G, grid = generate_rr_graph()

if INSERT_EMPTY_RING:
    grid_h += 2 * max([v[0] for v in V])
    grid_w += 2 * max([h[0] for h in H])
    print("After empty ring insertion: %d X %d" % (grid_w, grid_h))
    G, grid = generate_rr_graph()

if not ONLY_PAD:
    print("Generating architecture.\n")
    G, grid = generate_rr_graph(make_sb_clique = MAKE_SB_CLIQUE)
    criticalities = None
    potential_edge_delays = None
    if args.change_grid_dimensions is None:
        if args.load_mux_stack_order is not None:
            with open(args.load_mux_stack_order, "r") as inf:
                SB_MUX_ORDER = ast.literal_eval(inf.read())
        else:
            load_model = create_load_model()
            if FIRST_CLIQUE_ITER:
                criticalities = None
            else:
                #FIXME: See why it sometimes happens. It is related to a difference between
                #packing and placement block names.
                try:
                    criticalities = parse_routing.fetch_timing_data()
                except:
                    print "WARNING: Problem parsing criticalities. Passing all 1.0."
            optimize_pattern_layout(G, criticalities = criticalities)

        #We need to populate the segment_id dictionary before exporting the base costs:
        potential_edges = sorted(set([attrs["seg"] for u, attrs in G.nodes(data = True)\
                                 if u.startswith("potential_edge")])) if MAKE_SB_CLIQUE else None
        seg_text, seg_ids = export_segments(H, V, potential_edges)
        del seg_text

        potential_edge_delays = export_base_costs(G, "base_costs.dump", criticalities = criticalities)
    export_rr_graph(G, grid, args.arc_name.rsplit('.', 1)[0] + "_rr.xml", potential_edge_delays = potential_edge_delays)

#Just print the area data again.
active_w, active_h =  get_tile_dimensions(G)
txt = "\nActive dimensions: %d X %d nm\n" % (active_w, active_h)
metal_w, metal_h = get_metal_dimensions()
txt += "Metal dimensions: %d X %d nm\n\n" % (metal_w, metal_h)

cb_sizes, sb_sizes = export_mux_sizes(G)
mux_sizes = cb_sizes
mux_sizes.update(sb_sizes)
size_indexed = {}
for mux in mux_sizes:
    size = mux_sizes[mux]
    mux_id = mux.split("ble_%d_" % NEUTRAL_BLE)[1]
    try:
        size_indexed[size].append(mux_id)
    except:
        size_indexed.update({size : [mux_id]})

txt += "Multiplexer sizes:\n\n"
for size in sorted(size_indexed):
    txt += "%d: " % size
    for mux in sorted(size_indexed[size]):
        txt += "%s " % mux
    txt = txt[:-1] + "\n"

with open("tile_stats.log", "w") as outf:
    outf.write(txt[:-1])

print("\n" + txt)
