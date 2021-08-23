"""Explores switch-block patterns using the avalanche-cost method.

Parameters
----------
base_cost : float
    Initial base cost of the potential edges.
scaling_factor : int
    Negative exponent to which the base costs are to be raised.
avalanche_iter : int
    The desired number of iterations needed to reduce the cost of a potential edge to zero,
    assuming that historical and proportional factors are equal, that differential factor is zero, and that
    there is a constant utilization equal to the maximum utilization in the first routing iteration,
    when all base costs are reset to zero.
wd : str
    Working directory.
adoption_threshold : Optional[float], default = 1.2
    Threshold for utilization drop above which all switches are adopted into the pattern.
gnl : Optional[float], default = None
    Specifies that Gnl benchmarks with the given Rent's exponent should be used.
greedy : Optional[bool], default = False
    Specifies that we are performing a simple greedy search.
start_from_iter : Optional[int], default = None
    Specifies that the search should start from some already passed iteration.
    Useful to resume from crashes.

Returns
-------
None
"""

import time
import os
import argparse
import sys
sys.path.insert(0,'..')

import setenv
import check_rr_graph

parser = argparse.ArgumentParser()
parser.add_argument("--base_cost")
parser.add_argument("--scaling_factor")
parser.add_argument("--avalanche_iter")
parser.add_argument("--wd")
parser.add_argument("--adoption_threshold")
parser.add_argument("--gnl")
parser.add_argument("--greedy")
parser.add_argument("--start_from_iter")

args = parser.parse_args()

wd = os.getcwd()

os.system("mkdir ../%s" % args.wd)
os.system("cp -r . ../%s" % args.wd)
os.chdir("../%s" % args.wd)

CHECKS_ON = False

grid_sizes = {"alu4" : 15, "ex5p" : 12, "tseng" : 13}
benchmarks = ["alu4", "ex5p", "tseng"]
seeds = [19225, 25124, 43033, 50936, 5300]

gnl_grid_size = 39
GNL = args.gnl

if GNL is not None:
    benchmarks = list(sorted(["%s" % b.rsplit('.', 1)[0] for b in os.listdir("benchmarks/gnl_benchmarks/")\
                             if ("_p%s_size10000_" % GNL in b) and b.endswith(".blif")]))

    grid_sizes = {b : gnl_grid_size for b in benchmarks}

GREEDY = False
try:
    GREEDY = int(args.greedy)
except:
    pass

START_FROM_ITER = None
try:
    START_FROM_ITER = int(args.start_from_iter)
except:
    pass

reset_cost = 0.0 if GREEDY else 1e-5
TECH = '4'

lb_delays = {}

##########################################################################
def write_avalanche_conf():
    """Writes the avalanche configuration.

    Parameters
    ----------
    None

    Returns
    -------
    None
    """

    conf = {"base_cost" : 1,\
            "avalanche_p" : -1,\
            "avalanche_h" : -1,\
            "avalanche_d" : 0,\
            "scaling_factor" : int(args.scaling_factor),\
            "avalanche_iter" : int(args.avalanche_iter),\
           }

    variable_cost_target = 1e-12
    variable_cost = variable_cost_target / (10 ** (-1 * int(args.scaling_factor)))
    variable_cost = 0.0

    splitter_crit_exponent = 8.0
    splitter_target_crit_cost = 1e-12
    splitter_target_crit_cost /= 10 ** (-1 * int(args.scaling_factor))
    with open("avalanche.conf", "w") as outf:
        for k in sorted(conf):
            outf.write("%s %d\n" % (k, conf[k]))
        outf.write("variable_cost %f\n" % variable_cost)
        outf.write("reset_cost %f\n" % reset_cost)
        outf.write("splitter_crit_exponent %f\n" % splitter_crit_exponent)
        outf.write("splitter_target_crit_cost %f" % splitter_target_crit_cost)
##########################################################################

write_avalanche_conf()

arc_name = "agilex"

##########################################################################
def create_wafer(seed, first_clique_iter, replace_only = False):
    """Calls wafer creation.

    Parameters
    ----------
    seed : int
        The placement seed to be used. We use one for all circuits
        and change it at each global iteration.
    first_clique_iter : bool
        Specifies if this is the first clique iteration or not.
    replace_only : Optional[bool], default = False
        Tells the wafer generator to only replace the circuits.

    Returns
    -------
    None
    """

    if len(benchmarks) < 2:
        return

    first_clique_iter = 1 if first_clique_iter else 0
    replace_only = 1 if replace_only else 0
    circs = str([(b + ".blif", seed) for b in sorted(benchmarks,\
                                              key = lambda k : grid_sizes[k])])

    #Default parameters:
    wire_file = "arc_model/agilex/agilex.wire"
    padding_file = "arc_model/agilex/agilex_padding.log"

    size_args = ','.join([str(s) for s in sorted(grid_sizes.values())])

    wafer_switches = ["--wafer_w %d" % len(benchmarks),\
                      "--wafer_h 1",\
                      "--tech %s" % TECH,\
                      "--arc %s.xml" % arc_name,\
                      "--wire_file \"%s\"" % wire_file,\
                      "--padding_file \"%s\"" % padding_file,\
                      "--fpga_sizes \"%s\"" % size_args,\
                      "--merge_circs 1",\
                      "--circs \"%s\"" % circs,\
                      "--first_clique_iter %d" % first_clique_iter,\
                      "--replace_circs_only %d" % replace_only\
                     ]

    if args.adoption_threshold is not None:
        wafer_switches += ["--adoption_threshold %f" % float(args.adoption_threshold)]

    wafer_call =  ' '.join(["python -u create_wafer.py"] + wafer_switches)

    os.system(wafer_call)
##########################################################################

##########################################################################
def call_arc_gen(grid_w, grid_h, first_clique_iter = False, use_default = False):
    """Calls the architecture generation script.

    Parameters
    ----------
    grid_w : int
        Grid width. Might change due to physical length balancing.
    grid_h : int
        Grid height. Might change due to physical length balancing.
    first_clique_iter : Optional[bool], default = False
        Specifies that this is the first clique forming iteration and that all
        switches should be present, as well as that no VPR log should be read.
    use_default : Optional[bool], default = False
        Specifies that the default architecture is to be used,
        without any new SPICE simulations.
       
    Returns
    -------
    None
    """

    #Default parameters:
    in_arc_file = "arc_model/agilex/agilex.xml"
    arc_file = "agilex.xml"
    wire_file = "arc_model/agilex/agilex.wire"
    padding_file = "arc_model/agilex/agilex_padding.log"

    arc_gen_switches = ["--K 6",\
                        "--N 8",\
                        "--density 0.5",\
                        "--tech %s" % TECH,\
                        "--wire_file %s" % wire_file,\
                        "--import_padding %s" % padding_file,\
                        "--physical_square 1",\
                        "--make_sb_clique 1",\
                        "--arc_name %s.xml" % arc_name,\
                        "--grid_w %d" % grid_w,\
                        "--grid_h %d" % grid_h\
                       ]

    if first_clique_iter:
        arc_gen_switches.append("--first_clique_iter 1")
        if GNL and os.path.exists("gnl_init.pattern"):
            arc_gen_switches.append("--load_pattern_from_file gnl_init.pattern")
            arc_gen_switches.append("--pattern_is_init_only 1")

    if args.adoption_threshold is not None:
        if float(args.adoption_threshold) < 0:
            arc_gen_switches.append("--pick_top_only 1")
        else:
            arc_gen_switches.append("--adoption_threshold %f" % float(args.adoption_threshold))

    if GREEDY:
        arc_gen_switches.append("--greedy_switch_search 1")

    if use_default:
        arc_gen_switches.append("--change_grid_dimensions %s" % in_arc_file)
    
    call = ' '.join(["time python -u arc_gen.py"] + arc_gen_switches)

    os.system(call)
##########################################################################

##########################################################################
def call_vpr(arc, circ, seed = None):
    """Calls VPR

    Parameters
    ----------
    arc : str
        Name of the architecture.
    circ : str
        Name of the circuit.
    seed : Optional[int], default = None
        Placement seed. If none, no replacement is triggered.

    Returns
    -------
    None
    """

    vpr_arguments = ["%s.xml" % arc,\
                     "%s.blif" % circ\
                    ]
    
    vpr_switches = ["--route",\
                    "--read_rr_graph %s_rr.xml" % arc,\
                    "--route_chan_width 352",\
                    "--router_lookahead map",\
                    "--routing_failure_predictor off",\
                    "--max_router_iterations 300",\
                    "--incremental_reroute_delay_ripup %s" % ("off" if GNL else "on"),\
                    "--router_max_convergence_count 1",\
                    "--max_criticality %s" % ("0.00" if GNL else "0.99")\
                   ]

    if seed is not None:
        vpr_switches += ["--pack",\
                         "--place",\
                         "--import_place_delay_model %s_placement_delay.matrix" % arc,\
                         "--seed %d" % seed\
                        ]

    if CHECKS_ON:
        vpr_switches.append("--write_rr_graph check_rr.xml")

    run_req_files = ["%s.xml" % arc,\
                     "%s_rr.xml" % arc,\
                     "%s_placement_delay.matrix" % arc,\
                     "base_costs.dump",\
                     "%s%s%s.blif" % (("benchmarks" if circ != "wafer" else ''),\
                                       "gnl_benchmarks/" if GNL else '', circ),\
                     "%s.sdc" % circ\
                    ]

    if seed is None:
        run_req_files += ["%s.net" % circ,\
                          "%s.place" % circ\
                         ]

    sandbox = "sandbox_%s_%s" % (os.getcwd().rsplit('/', 1)[1], str(time.time()))
    os.system("mkdir %s" % (os.environ["VPR_RUN_PATH"] % sandbox))

    for f in run_req_files:
        os.system("cp %s %s/" % (f, os.environ["VPR_RUN_PATH"] % sandbox))

    vpr_call = "time %s" % (os.environ["VPR"] % (sandbox, ' '.join(vpr_arguments + vpr_switches)))
    os.system(vpr_call)

    for f in os.listdir(os.environ["VPR_RUN_PATH"] % sandbox):
        if f in run_req_files:
            continue
        os.system("cp %s/%s ./" % (os.environ["VPR_RUN_PATH"] % sandbox, f))

    os.system("rm -rf %s" % (os.environ["VPR_RUN_PATH"] % sandbox))
##########################################################################

##########################################################################
def multiply_base_costs(fac):
    """Multiplies all the base costs by a constant factor.

    Parameters
    ----------
    fac : float
        Multiplication factor.

    Returns
    -------
    None
    """

    with open("base_costs.dump", "r") as inf:
        lines = inf.readlines()

    txt = ''.join(lines[:3])
    for line in lines[3:]:
        cost = line.split()[1]
        line = line.replace(cost, "%f" % (float(cost) * fac))
        txt += line

    with open("base_costs.dump", "w") as outf:
        outf.write(txt)
##########################################################################

##########################################################################
def get_arc_delays():
    """Parses the wire delays from the architecture file.

    Parameters
    ----------
    None

    Returns
    -------
    str
        A printed delay dictionary.
    """

    with open("%s.xml" % arc_name, "r") as inf:
        lines = inf.readlines()

    get_attr = lambda line, attr : line.split("%s=\"" % attr, 1)[1].split('"', 1)[0]
    delays = {}
    for line in lines:
        if "<switch " in line:
            wire = get_attr(line, "name")
            td = float(get_attr(line, "Tdel"))
            delays.update({wire : td})

    txt = ""
    for wire in sorted(delays):
        txt += "%s %g\n" % (wire, delays[wire])
      
    return txt[:-1]
##########################################################################

##########################################################################
def trim_vpr_log(icnt):
    """Trims and saves the vpr log and copies the potential edge store.

    Parameters
    ----------
    icnt : int
        Iteration count.

    Returns
    -------
    None
    """

    with open("vpr_stdout.log", "r") as inf:
        lines = inf.readlines()
    print "About to trim VPR log. OK?"

    txt = ""
    rd = False
    for line in lines:
        if line.startswith("Confirming router"):
            rd = True
            continue
        if not rd:
            continue
        txt += line

    with open("vpr_iter_%d.log" % icnt, "w") as outf:
        outf.write(txt)

    os.system("cp stored_edges.save stored_edges_iter_%d.save" % icnt)
    
    with open("arc_delays_iter_%d.log" % icnt, "w") as outf:
        outf.write(get_arc_delays())
    os.system("cp arc_delays_iter_%d.log arc_delays.log" % icnt)

    os.system("cp tile_stats.log tile_stats_iter_%d.log" % icnt)
##########################################################################

##########################################################################
def parse_delays(log):
    """Parses delays from a VPR log.

    Parameters
    ----------
    log : str
        Name of the VPR log file.

    Returns
    -------
    Dict[str, float]
        Delay dictionary, for all clock domains.
    """
    
    with open(log, "r") as inf:
        lines = inf.readlines()

    delays = {}
    for line in reversed(lines):
        if "worst setup slack" in line and not "Int" in line:
            clk_from = line.split()[0]
            clk_to = line.split()[2]
            if clk_from != clk_to:
                print "Error parsing intradomain delays!"
                raise ValueError
            td = -1 * float(line.split()[-2])
            delays.update({clk_from : td})
        elif "Setup Worst Negative Slack" in line:
            wns = -1 * float(line.split()[-2])
            if not delays:
                if len(benchmarks) == 1:
                    delays.update({"clk" : wns})
                    return delays
                print "Clock domains not constrained properly for the wafer. Check if .sdc is copied."
                raise ValueError
            if not wns in delays.values():
                print "WNS not found among intradomain delays!"
                raise ValueError
            return delays
##########################################################################

##########################################################################
def get_delay_lower_bounds():
    """Obtains lower bounds for delays of all clock domains, over all seeds.
    Populates a global dictionary.

    Parameters
    ----------
    None

    Returns
    -------
    None
    """

    global lb_delays

    os.system("cp base_costs.dump base_costs.prev")
    multiply_base_costs(0.0)
    #Lower bound is obtained when all switches are free.

    for seed in seeds:
        if len(benchmarks) > 1:
            call_vpr("wafer", "wafer", seed = None)
        else:
            create_wafer(seed, first_clique_iter = 1, replace_only = 1)
            call_vpr(arc_name, benchmarks[0], seed = seed)
        lb_delays.update({seed : parse_delays("vpr_stdout.log" )})

    #Restore the original seed.
    if len(benchmarks) > 1:
        create_wafer(seeds[0], first_clique_iter = 1, replace_only = 1)

    os.system("mv base_costs.prev base_costs.dump")
##########################################################################

init = True
icnt = 0
if START_FROM_ITER is not None:
    icnt = START_FROM_ITER
    store = "stored_edges_iter_%d.save" % icnt
    os.system("cp %s stored_edges.save" % store)
    vpr_log = "vpr_iter_%d.log" % icnt
    os.system("cp %s vpr_stdout.log" % vpr_log)
    delay_file = "arc_delays_iter_%d.log" % icnt
    os.system("cp %s arc_delays.log" % delay_file)
    init = False

if init:
    os.system("rm stored_edges.save")
done = False
mult_fac = 1.0
mult_fac_inc = 1.0
init_reset_cost = reset_cost
while not done:
    circ = benchmarks[icnt % len(benchmarks)]
    seed = seeds[icnt % len(seeds)]
    icnt += 1
    os.system("rm -f %s*.lz4" % arc_name)

    if not GNL:
        grid_w = grid_h = min(grid_sizes.values())
    else:
        grid_w = grid_h = grid_sizes[circ]
    call_arc_gen(grid_w, grid_h, first_clique_iter = init)

    with open("base_costs.dump", "r") as inf:
        lines = inf.readlines()
    
    if lines[2][0] == '0':
        done = True

    os.system("lz4 -d %s_rr.xml.lz4" % arc_name)
    
    if not GNL:
        circs = [(b + ".blif", seed) for b in benchmarks]
        circs = str(sorted(circs, key = lambda c : grid_sizes[c[0].rsplit(".blif", 1)[0]]))
        create_wafer(seed, init)

    multiply_base_costs(mult_fac)
    reset_cost = init_reset_cost * mult_fac
    write_avalanche_conf()
    mult_fac *= mult_fac_inc

    if not GNL and (len(benchmarks) > 1):
        call_vpr("wafer", "wafer", seed = None)
    else:
        call_vpr(arc_name, circ, seed = seed)

    if CHECKS_ON:
        if len(benchmarks) == 1:
            reload(check_rr_graph)
            check_rr_graph.check_rr_all("check_rr.xml", None if init else "stored_edges.save")
        else:
            print "WARNING: RR-graph checks work only for individual FPGAs, not wafers."
            print "Runing checks on individual FPGAs. Pre-VPR RR-graphs will have to be used."
            rr_names = ["%s_W%d_H%d_rr.xml" % (arc_name, grid_sizes[b], grid_sizes[b]) for b in benchmarks]
            for rr_name in rr_names:
                print "Checking %s..." % rr_name
                reload(check_rr_graph)
                check_rr_graph.check_rr_all(rr_name, None if init else "stored_edges.save")
    
    init = False

    print "Iteration %d done." % icnt
    print "Seed", seed
    trim_vpr_log(icnt)

os.chdir(wd)
