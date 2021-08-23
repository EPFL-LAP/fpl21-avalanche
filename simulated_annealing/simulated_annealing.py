"""Explores switch-block patterns using the avalanche-cost method.

Parameters
----------
wd : str
    Working directory.

Returns
-------
None
"""

import time
import os
import argparse
import random
import math
import sys
sys.path.insert(0,'..')

import setenv
import check_rr_graph

from parallelize import Parallel

parser = argparse.ArgumentParser()
parser.add_argument("--wd")

args = parser.parse_args()

wd = os.getcwd()

os.system("mkdir ../%s" % args.wd)
os.system("cp -r . ../%s" % args.wd)
os.chdir("../%s" % args.wd)

grid_sizes = {"alu4" : 15, "ex5p" : 12, "tseng" : 13}
benchmarks = ["alu4", "ex5p", "tseng"]
seeds = [19225, 25124, 43033, 50936, 5300]

reset_cost = 1e-5
TECH = '4'

arc_name = "agilex"

FAILURE_CACHE = []
PREV_SPICED = None
RESPICE_THR = 5

########################################################################## 
def pattern_to_bit_string(pattern):
    """Converts a pattern to a bit string, by using indicators for switch presence.

    Parameters
    ----------
    pattern : str
        Textual description of a pattern.

    Returns
    -------
    Tuple[int]
        The bitstring.
    """
    
    bits = []
    for line in pattern.splitlines():
        if line.startswith('~'):
            bits.append(0)
        else:
            bits.append(1)

    return tuple(bits)
########################################################################## 

init_pattern_file = "sa_init.pattern"

init_edge_list = []
with open(init_pattern_file, "r") as inf:
    init_pattern_text = inf.read()

for line in init_pattern_text.splitlines():
    e = line.strip()
    if line[0] == '~':
        e = e[1:]
    init_edge_list.append(e)

pattern = pattern_to_bit_string(init_pattern_text)


##########################################################################
def call_arc_gen(grid_w, grid_h, pattern_file, sb_mux_order = None, in_arc_file = None):
    """Calls the architecture generation script.

    Parameters
    ----------
    grid_w : int
        Grid width. Might change due to physical length balancing.
    grid_h : int
        Grid height. Might change due to physical length balancing.
    pattern_file : str
        Name of the file holding the current pattern description.
    sb_mux_order : Optional[str], default = None
        Specifies the file from which to load the sb-mux order.         
    in_arc_file : Optional[str], default = None
        Specifies that a prior architecture is to be used,
        without any new SPICE simulations.
       
    Returns
    -------
    None
    """

    #Default parameters:
    arc_file = "agilex.xml"
    wire_file = "grid_example/agilex/agilex.wire"
    padding_file = "grid_example/agilex/agilex_padding.log"

    arc_gen_switches = ["--K 6",\
                        "--N 8",\
                        "--density 0.5",\
                        "--tech %s" % TECH,\
                        "--wire_file %s" % wire_file,\
                        "--import_padding %s" % padding_file,\
                        "--physical_square 1",\
                        "--arc_name %s.xml" % arc_name,\
                        "--grid_w %d" % grid_w,\
                        "--grid_h %d" % grid_h,\
                        "--load_pattern_from_file %s" % pattern_file,\
                        ]

    if sb_mux_order is not None:
        arc_gen_switches.append("--load_mux_stack_order %s" % sb_mux_order)

    if in_arc_file is not None:
        arc_gen_switches.append("--change_grid_dimensions %s" % in_arc_file)

    os.system(' '.join(["time python -u ../avalanche/arc_gen.py"] + arc_gen_switches))
##########################################################################

##########################################################################
def evaluate(pattern, force_respice = False, sb_mux_order = None):
    """Generates a new architecture, then runs PnR and retrieves the geomean delay.

    Parameters
    ----------
    pattern : Tuple[int]
        A bitstring encoding the pattern to be evaluated.
    force_respice : Optional[bool], default = False
        Specifies that the new architecture must be respiced, regardless of
        whether the ordinary respicing condition was met or not.
    sb_mux_order : Optional[str], default = None
        Specifies the file from which to load the sb-mux order.         

    Returns
    -------
    float
        Tile width in um.
    float
        Geomean routed delay, for the specified benchmark set.
    """

    if check_subset(pattern):
        with open("sa.log", "a") as outf:
            outf.write("cache-eliminated\n")
        return "failed", "failed"

    global PREV_SPICED
    inherit_delays = "prev_spiced.xml"
    if (PREV_SPICED is None) or (sb_mux_order is None)\
       or force_respice or (bits_diff(pattern, PREV_SPICED) > RESPICE_THR):
        inherit_delays = None
        PREV_SPICED = pattern
        with open("sa.log", "a") as outf:
            outf.write("respiced%s\n" % (" reoptimized" if sb_mux_order is None else ''))

    pattern_txt = ""
    for i, e in enumerate(init_edge_list):
        pattern_txt += "%s%s\n" % ('~' if pattern[i] == 0 else '', e)

    with open("cur.pattern", "w") as outf:
        outf.write(pattern_txt[:-1])

    os.system("rm -f %s*.lz4" % arc_name)
    grid_w = grid_h = max(grid_sizes.values())
    call_arc_gen(grid_w, grid_h, pattern_file = "cur.pattern", sb_mux_order = sb_mux_order, in_arc_file = inherit_delays)
    os.system("rm -f %s_rr.xml" % arc_name)
    os.system("lz4 -d %s_rr.xml.lz4" % arc_name)

    if inherit_delays is None:
        os.system("cp %s.xml prev_spiced.xml" % arc_name)
        os.system("cp %s_rr.xml prev_spiced_rr.xml" % arc_name)
    
    vpr_calls = []
    vpr_logs = []
    for circ in benchmarks:
        for seed in seeds:
            vpr_calls.append("time python -u run_pnr.py --arc %s --circ %s --seed %d" % (arc_name, circ, seed))
            vpr_logs.append("vpr_%s_%d.log" % (circ, seed))

    
    max_cpu = 15
    sleep_interval = 5

    runner = Parallel(max_cpu, sleep_interval)
    runner.init_cmd_pool(vpr_calls)
    runner.run()
    

    global FAILURE_CACHE
    cpds = {}
    for log in vpr_logs:
        circ = log.split('_')[1]
        try:
            with open(log, "r") as inf:
                lines = inf.readlines()
        except:
            FAILURE_CACHE.append(pattern)
            return "failed", "failed"
        cpd = None
        try:
            for line in lines:
                if "Final critical path:" in line:
                    cpd = float(line.split()[3])
                    print cpd
                    break
        except:
            FAILURE_CACHE.append(pattern)
            return "failed", "failed"
        if cpd is None:
            FAILURE_CACHE.append(pattern)
            return "failed", "failed"
        try:
            cpds[circ].append(cpd)
        except:
            cpds.update({circ : [cpd]})

    geom = 1.0
    for circ in cpds:
        cpds[circ].sort()
        geom *= cpds[circ][len(cpds[circ]) / 2]

    geom **= 1.0 / len(cpds)
        
    with open("tile_stats.log", "r") as inf:
        lines = inf.readlines()

    for line in lines:
        if "Active dimensions" in line:
            area = float(line.split()[-2]) * float(line.split()[-4]) / 1000000.0
            break
    
    return geom, area
##########################################################################

##########################################################################
def evaluate_move(prev_td, prev_area, new_td, new_area, T):
    """Evaluates the move

    Parameters
    ----------
    prev_td : float
        Previous timing cost.
    prev_area : float
        Previous area cost.
    new_td : float
        New timing cost.
    new_area : float
        New area cost.
    T : float
        Temperature.
    
    Returns
    -------
    bool
        Indicator of acceptance.
    """

    td_vs_area = 0.5

    delta_td = (new_td - prev_td) / float(prev_td)
    delta_area = (new_area - prev_area) / float(prev_area)
    delta_c = td_vs_area * delta_td + (1.0 - td_vs_area) * delta_area

    if delta_c < 0:
        return True

    energy_threshold = random.random()
    energy = math.exp(-1 * delta_c / T)

    return energy > energy_threshold
##########################################################################

##########################################################################
def flip(ind, cost):
    """Flips the usage of one edge.

    Parameters
    ----------
    ind : int
        Index of the line describing the edge in the edge store.
    cost : int
        Current cost.
    """

    if lines[ind][0] == '~':
        lines[ind] = lines[ind][1:]
        cost += 1
    else:
        lines[ind] = '~' + lines[ind]
        cost -= 1

    with open("stored_edges.save", "w") as outf:
        outf.write(''.join(lines))

    return move, cost
##########################################################################

##########################################################################
def log_move(move, T, prev_td_cost, td_cost, prev_area_cost, area_cost, status):
    """Logs the given move.

    Parameters
    ----------
    move : str
        Original line being flipped.
    T : float
        Current temperature.
    prev_cost : int
        Previous cost.
    cost : int
        Resulting cost.
    status : str
        Status of the move, from {accepted, unroutable, expensive}

    Returns
    -------
    None
    """

    if td_cost == "failed":
        td_cost = float('inf')
        status = "unroutable"
    if area_cost == "failed":
        area_cost = float('inf')
        status = "unroutable"

    with open("sa.log", "a") as outf:
        routable = "yes" if status == "accepted" else ("no" if status == "unroutable" else "NA")
        accepted = "yes" if status == "accepted" else "no"
        line = "%s\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%s\t%s\n" % (move, T, prev_td_cost, td_cost, prev_area_cost, area_cost,\
                                                               routable, accepted)
        outf.write(line)
        print line
##########################################################################

########################################################################## 
def a_subset_b(bits_a, bits_b):
    """Checks if one bitstring encodes a proper subset of another.

    Parameters
    ----------
    bits_a : Tuple[int]
        The first string.
    bits_b : Tuple[int]
        The second string.

    Returns
    -------
    bool
        True if a subset of b.
    """
    
    subset_eq =  all(bits_a[i] <= bits_b[i] for i in range(0, len(bits_a)))

    return subset_eq and (sum(bits_a) < sum(bits_b))
########################################################################## 

##########################################################################
def check_subset(cur):
    """Checks if the current set of used edges is a subset of any
    set known to be unroutable. In that case, it is also unroutable by construction.

    Parameters
    ----------
    cur : Tuple[int]
        Bitstring encoding the current set of used edges.

    Returns
    -------
    bool
        True if subset, False otherwise.
    Set[str]
        Current used set.
    """

    return any(a_subset_b(cur, f) for f in FAILURE_CACHE)
########################################################################## 

########################################################################## 
def bits_diff(bits_a, bits_b):
    """Returns the Hamming distance between two bitstrings.

    Parameters
    ----------
    bits_a : Tuple[int]
        The first string.
    bits_b : Tuple[int]
        The second string.

    Returns
    -------
    int
        The number of bits that differ.
    """

    return sum([abs(bits_a[i] - bits_b[i]) for i in range(0, len(bits_a))]) 
########################################################################## 
    
########################################################################## 
def fair_ind(pattern):
    """Returns a random index from the edge store, with an equal probability
    of flipping a used and an unused edge.

    Parameters
    ----------
    pattern : Tuple[int]
        The current pattern.

    Returns
    -------
    Tuple[int]
        Modified pattern.
    str
        Attempted move.
    """

    used = [i for i in range(0, len(pattern)) if pattern[i] == 1]
    unused = [i for i in range(0, len(pattern)) if pattern[i] == 0]
    
    take_used = random.randint(0, 1)
    flip_ind = None
    if take_used:
        try:
            flip_ind = random.choice(used)
        except:
            flip_ind =  random.choice(unused)
    else:
        try:
            flip_ind = random.choice(unused)
        except:
            flip_ind = random.choice(used)


    updated = list(pattern)
    updated[flip_ind] = 1 - updated[flip_ind]

    move = "(%d) %s : %d -> %d" % (sum(updated), init_edge_list[flip_ind], pattern[flip_ind], 1 - pattern[flip_ind])

    return tuple(updated), move
##########################################################################

max_iter = 100
init_temp = 0.02
temperature = init_temp
temperature_decrease = 0.95
moves_per_temperature = 100
accepted = 0
rejected = 0
random.seed(0)


prev_td_cost, prev_area_cost = evaluate(pattern)
log_move("Search start", init_temp, -1, prev_td_cost, -1, prev_area_cost, "accepted")


new_pattern, move = fair_ind(pattern)
new_td_cost, new_area_cost = evaluate(new_pattern, sb_mux_order = "sb_mux.order")


for icnt in range(0, max_iter):
    for mov in range(0, moves_per_temperature):
        new_pattern, move = fair_ind(pattern)
        new_td_cost, new_area_cost = evaluate(new_pattern, sb_mux_order = "sb_mux.order")
        if new_td_cost == "failed" or new_area_cost == "failed":
            log_move(move, temperature, prev_td_cost, new_td_cost, prev_area_cost, new_area_cost, "rejected")
            rejected += 1
            continue
        accept = evaluate_move(prev_td_cost, prev_area_cost, new_td_cost, new_area_cost, temperature)
        if accept:
            log_move(move, temperature, prev_td_cost, new_td_cost, prev_area_cost, new_area_cost, "accepted")
            accepted += 1
            pattern = new_pattern
            prev_td_cost = new_td_cost
            prev_area_cost = new_area_cost
        else:
            log_move(move, temperature, prev_td_cost, new_td_cost, prev_area_cost, new_area_cost, "rejected")
            rejected += 1
                
    temperature *= temperature_decrease
    prev_td_cost, prev_area_cost = evaluate(pattern)
os.chdir(wd)
