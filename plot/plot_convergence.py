import os
import argparse
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec

parser = argparse.ArgumentParser()
parser.add_argument("--data_dir")

args = parser.parse_args()


fig, ax = plt.subplots(4, 1, figsize = (8.3, 11.7), gridspec_kw={
                           'width_ratios': [1],
                           'height_ratios': [4, 4, 5, 0.5]})

vpr_logs = sorted([f for f in os.listdir(args.data_dir) if f.startswith("vpr_iter")],\
                   key = lambda l : int(l.split('_')[-1].split(".log")[0]))
print vpr_logs

##########################################################################
def parse_congestion(vpr_log):
    """Returns the congestion sequence from a VPR log.

    Parameters
    ----------
    vpr_log : str
        Path to the VPR log.

    Returns
    -------
    List[int]
        List of congested nodes, per routing iteration.
    """

    with open(vpr_log, "r") as inf:
        lines = inf.readlines()

    
    cong = []
    for line in lines:
        if "Restoring best routing" in line:
            break
        elif "%" in line:
            cong.append(int(line.split()[7]))

    return cong
##########################################################################

##########################################################################
def parse_cpd(vpr_log):
    """Returns the critical path delay sequence from a VPR log.

    Parameters
    ----------
    vpr_log : str
        Path to the VPR log.

    Returns
    -------
    List[float]
        List of critical path delays, per routing iteration.
    """

    with open(vpr_log, "r") as inf:
        lines = inf.readlines()

    
    cpd = []
    for line in lines:
        if "Restoring best routing" in line:
            break
        elif "%" in line:
            cpd.append(float(line.split()[-6]))

    for line in lines:
        if "Critical path:" in line:
            cpd.append(float(line.split()[-2]))

    return cpd
##########################################################################

##########################################################################
def parse_costs(vpr_log, curves):
    """Parses potential edge cost evolution from a VPR log.

    Parameters
    ----------
    vpr_log : str
        Path to the VPR log.
    curves : int
        Number of curves to return.
    
    Returns
    -------
    List[List[float]]
        A list of cost distributions, one per iteration.
    int
        Tail cut index.
    """

    with open(vpr_log, "r") as inf:
        lines = inf.readlines()

    costs = []
    for line in lines:
        if line.startswith("Edge-splitter costs"):
            costs.append([])
        elif "->" in line:
            cost = float(line.split('(')[-1].split(')')[0])
            cost = float(line.split('(')[2].split(',')[0])
            costs[-1].append(cost)
 
    for i in range(0, len(costs)):
        costs[i].sort()
  
    THR = 20 
    last_zero = []
    for cost in costs:
        last_zero.append(0)
        for i in range(0, len(cost)):
            if cost[i] < THR:
                last_zero[-1] = i
    last_zero = min(last_zero)
    print last_zero

    for i in range(0, len(costs)):
        costs[i] = costs[i][last_zero:]

    return [(i, costs[i]) for i in range(1, len(costs), len(costs) / curves)], last_zero
##########################################################################

colormap = plt.cm.Set1

all_cong = []
for log in vpr_logs:
    all_cong.append(parse_congestion("%s/%s" % (args.data_dir, log)))
all_cpd = []
for log in vpr_logs:
    all_cpd.append(parse_cpd("%s/%s" % (args.data_dir, log)))

max_len = max([len(c) for c in all_cong])
xs = range(0, max_len)

cpd_max_len = max([len(d) for d in all_cpd])
cpd_xs = range(0, cpd_max_len)
for i, cpd in enumerate(all_cpd):
    all_cpd[i] += [cpd[-1] for j in range(len(cpd), cpd_max_len)]

ax[0].set_color_cycle([colormap(i) for i in np.linspace(0, 1, len(all_cong))])
ax[1].set_color_cycle([colormap(i) for i in np.linspace(0, 1, len(all_cpd))])

for i, c in enumerate(all_cong):
    print i
    all_cong[i] += [0 for j in range(0, max_len - len(c))]
    ax[0].plot(xs, all_cong[i], label = "G. iter. #%d" % i)
    all_cpd[i] += [all_cpd[i][-1] for j in range(0, max_len - len(all_cpd[i]))]
    ax[1].plot(cpd_xs, all_cpd[i], label = "G. iter. #%d" % i)

ax[0].set_title("Congestion Convergence")
ax[0].set_xlabel("Routing iteration")
ax[0].set_ylabel("# Oversused nodes")

ax[1].set_title("Delay Convergence")
ax[1].set_xlabel("Routing iteration")
ax[1].set_ylabel("CPD [ns]")

CURVES = 10
costs, last_zero = parse_costs("%s/%s" % (args.data_dir, vpr_logs[0]), CURVES)
ax[2].set_color_cycle([colormap(i) for i in np.linspace(0, 1, CURVES)])
for i, c in enumerate(costs):
    ax[2].plot(range(0, len(c[1])), c[1], label = "Iter. #%d/%d" % (c[0], len(all_cong[0])))

ticks = ax[2].get_xticks()
ax[2].set_xticklabels([str(int(i) + last_zero) for i in ticks])
ax[2].set_title("Cost Evolution in Global Iteration #1")
ax[2].set_xlabel("Potential edge index")
ax[2].set_ylabel("Utilization [#SBs out of 1614]")

ax[0].legend(loc = "upper right")
#ax[1].legend(loc = "upper right")
ax[2].legend(loc = "upper left")

##########################################################################
def get_description():
    """Returns the description of the run.

    Parameters
    ----------
    None

    Returns
    -------
    str
        Run description.
    int
        Total number of accepted edges.
    """

    scaling_factor = int(args.data_dir.split("sf")[1].split('_')[0])
    iter_to_zero = int(args.data_dir.split('i')[1].split('/')[0])

    with open("%s/vpr_iter_1.log" % args.data_dir, "r") as inf:
        lines = inf.readlines()

    for line in lines:
        if line.startswith("avalanche_p"):
            avalanche_p = float(line.split()[-1])
        elif line.startswith("avalanche_h"):
            avalanche_h = float(line.split()[-1])
        elif line.startswith("avalanche_d"):
            avalanche_d = float(line.split()[-1])

    accepted = {}
    rejected = {}
    for f in os.listdir(args.data_dir):
        if not f.startswith("stored_edges_iter"):
            continue
        with open("%s/%s" % (args.data_dir, f), "r") as inf:
            lines = inf.readlines()
        itr = int(f.split('_')[-1].split('.')[0])
        accepted.update({itr : 0})
        rejected.update({itr : 0})
        for line in lines:
            if line.startswith('~'):
                rejected[itr] += 1
            else:
                accepted[itr] += 1

    accepted.update({1 : 0})
    rejected.update({1 : 0})
    acceptance_stats = "Accepted/Rejected(iter): "
    for i in sorted(accepted):
        if i == 1:
            continue
        acceptance_stats += "%d/%d, " % (accepted[i] - accepted[i - 1], rejected[i] - rejected[i - 1])

    txt = "\nStart cost = 10^(-%d); iterations to zero = %d; total accepted = %d\n\n"\
        % (scaling_factor, iter_to_zero, max(accepted.values()))
    txt += "avalanche_p, _h, _d = (%.2g, %.2g, %.2g)\n\n" % (avalanche_p, avalanche_h, avalanche_d)
    txt += acceptance_stats[:-2]

    return txt, max(accepted.values())
##########################################################################         

desc, total_accepted = "", 0 #get_description()
ax[3].text(0, 0, desc, fontsize = 10)
ax[3].set_ylim(0, 0.6)
ax[3].axis("off")

plt.subplots_adjust(wspace=0, hspace=0)

fig.tight_layout()

plt.savefig("%s.pdf" % args.data_dir[:-1])

print "Total accepted: %d" % total_accepted
