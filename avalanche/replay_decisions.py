"""Replays the decisions to check if they are correct.
"""

import os

##########################################################################
def check_iter(icnt):
    """Checks whether the decision in a given iteration is consistent with
    the logs and the previous delays.

    Parameters
    ----------
    icnt : int
        Iteration count.

    Returns
    -------
    int
        Minimum utilization at acceptance.

    Raises
    ------
    ValueError
        If the choice is not consistent.
    """

    stored = []
    store = "stored_edges_iter_%d.save" % (icnt + 1)
    try:
        with open(store, "r") as inf:
            lines = inf.readlines()

        stored = [line.strip() for line in lines if line[0] != '~']
    except:
        pass

    stored.sort()

    prev_stored = []
    prev_store = "stored_edges_iter_%d.save" % (icnt)
    try:
        with open(prev_store, "r") as inf:
            lines = inf.readlines()

        prev_stored = [line.strip() for line in lines if line[0] != '~']
    except:
        pass

    chosen = [e for e in stored if not e in prev_stored]

    vpr_log = "vpr_iter_%d.log" % icnt

    usages = {}
    with open(vpr_log, "r") as inf:
        lines = inf.readlines()

    rd = 0
    for line in lines:
        if "Total routing area" in line:
            rd += 1
            continue
        if rd == 1 and "Edge-splitter" in line:
            rd += 1
            continue
        if rd < 2:
            continue

        try:
            e = line.split('(', 1)[1].split(')', 1)[0]
            cost = float(line.split(" -> ", 1)[1].split('(', 1)[0])
            usage = int(line.split("): (", 1)[1].split(',', 1)[0])
            usages.update({e : (usage, cost)})
        except:
            break

    free = [e for e in usages if usages[e][1] == 0.0]
    free.sort()

    if free:
        if free != chosen:
            print "Expected decision =", free
            print "Stored decisions =", chosen
            raise ValueError
        else:
            min_util = min(free, key = lambda k : usages[k][0])
            return min_util, usages[min_util]

    max_util = max(usages.values(), key = lambda u : u[0])[0]
    if max_util == 0:
        candidates = []
    else:
        candidates = [e for e in usages if usages[e][0] == max_util]

    if len(candidates) > 1:
        print "Tied", len(candidates)
        td_dict = {}
        delay_file = "arc_delays_iter_%d.log" % icnt
        with open(delay_file, "r") as inf:
            lines = inf.readlines()

        for line in lines:
            td_dict.update({line.split()[0] : float(line.split()[1])})
        
        choice = [min(candidates, key = lambda c : td_dict[c])]
    else:
        choice = candidates

    if choice != chosen:
        print "Expected decision =", choice
        print "Stored decisions =", chosen
        raise ValueError
    elif max_util == 0:
            return [], -1
    else:
        return choice[0], usages[choice[0]]
##########################################################################

max_iter = 85

for icnt in range(1, max_iter):
    print icnt, check_iter(icnt)
    
