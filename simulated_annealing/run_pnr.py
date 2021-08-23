import time
import os
import argparse
import sys
sys.path.insert(0,'..')

import setenv

parser = argparse.ArgumentParser()
parser.add_argument("--arc")
parser.add_argument("--circ")
parser.add_argument("--seed")

args = parser.parse_args()

##########################################################################
def call_vpr(arc, circ, seed):
    """Calls VPR

    Parameters
    ----------
    arc : str
        Name of the architecture.
    circ : str
        Name of the circuit.
    seed : int
        Placement seed.

    Returns
    -------
    None
    """

    vpr_arguments = ["%s.xml" % arc,\
                     "%s.blif" % circ\
                    ]
 
    vpr_switches = ["--pack",\
                    "--place",\
                    "--import_place_delay_model %s_placement_delay.matrix" % arc,\
                    "--seed %d" % seed\
                   ]
   
    vpr_switches += ["--route",\
                     "--read_rr_graph %s_rr.xml" % arc,\
                     "--route_chan_width 352",\
                     "--router_lookahead map",\
                     "--routing_failure_predictor off",\
                     "--max_router_iterations 300",\
                     "--incremental_reroute_delay_ripup on",\
                     "--router_max_convergence_count 1",\
                     "--max_criticality %s" % "0.99"\
                    ]

    run_req_files = ["%s.xml" % arc,\
                     "%s_rr.xml" % arc,\
                     "%s_placement_delay.matrix" % arc,\
                     "base_costs.dump",\
                     "%s.blif" % circ\
                    ]

    sandbox = "sandbox_%s_%s_%d_%s" % (os.getcwd().rsplit('/', 1)[1], circ, seed, str(time.time()))
    os.system("mkdir %s" % (os.environ["VPR_RUN_PATH"] % sandbox))

    for f in run_req_files:
        os.system("cp %s %s/" % (f, os.environ["VPR_RUN_PATH"] % sandbox))

    vpr_call = "time %s" % (os.environ["VPR"].replace(" -it ", " -i ") % (sandbox, ' '.join(vpr_arguments + vpr_switches)))
    print vpr_call
    os.system(vpr_call)

    os.system("cp %s/vpr_stdout.log ./vpr_%s_%d.log" % (os.environ["VPR_RUN_PATH"] % sandbox, circ, seed))

    os.system("rm -rf %s" % (os.environ["VPR_RUN_PATH"] % sandbox))
##########################################################################

call_vpr(args.arc, args.circ, int(args.seed))
