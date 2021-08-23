"""Checks the avalanche costs from a VPR log.

Parameters
----------
log : str
    Name of the log file

Returns
-------
bool
    True if legal False otherwise.
"""

import os
import argparse
from collections import namedtuple

parser = argparse.ArgumentParser()
parser.add_argument("--log")

args = parser.parse_args()

Cost = namedtuple("Cost", ["index", "name", "cur_usage", "hist_usage",\
                           "prev_usage", "stored_cost", "derived_cost", "start_cost"])
##########################################################################
def parse_line(line):
    """Parses a line holding a cost.

    Parameters
    ----------
    line : str
        A line from the VPR log.

    Returns
    -------
    Cost
        Parsed cost.
    """

    if "avalanche_p" in line:
        global avalanche_p
        avalanche_p = float(line.split()[-1])
    elif "avalanche_h" in line:
        global avalanche_h
        avalanche_h = float(line.split()[-1])
    elif "avalanche_d" in line:
        global avalanche_d
        avalanche_d = float(line.split()[-1])

    if not " -> " in line:
        return None

    index = int(line.split('(', 1)[0])
    name = line.split('(', 1)[1].split(')')[0]

    usages = line.split(": (", 1)[1].split(')', 1)[0].split(", ")
    cur_usage = int(usages[0])
    hist_usage = int(usages[1])
    prev_usage = int(usages[2])

    costs = line.split("-> ", 1)[1]
    stored_cost = float(costs.split('(', 1)[0])
    derived_cost = float(costs.split('(', 1)[1].split(')', 1)[0])
    start_cost = float(costs.split('/')[1].strip())

    cost = Cost(index = index, name = name, cur_usage = cur_usage, hist_usage = hist_usage,\
                prev_usage = prev_usage, stored_cost = stored_cost, derived_cost = derived_cost,\
                start_cost = start_cost)

    return cost
##########################################################################

avalanche_p = avalanche_h = avalanche_d = -1


with open(args.log, "r") as inf:
    lines = inf.readlines()

cost_dict = {}
for line in lines:
    cost = parse_line(line)
    if cost is None:
        continue
    try:
        cost_dict[cost.index].append(cost)
    except:
        cost_dict.update({cost.index : [cost]})

##########################################################################
def check_costs(cost_dict):
    """Checks the cost dictionary for legality.

    Parameters
    ----------
    cost_dict : Dict[int, Cost]
        Cost dictionary.

    Returns
    -------
    bool
        True if legal else False.

    Raises
    ------
    ValueError
        If anything is not legal.
    """

    #------------------------------------------------------------------------#
    def compute_cur_cost(cost):
        """Computes the current cost based on the usage statistics.

        Parameters
        ----------
        cost : Cost
            Cost record.

        Returns
        -------
        float
            The computed current cost.

        Raises
        ------
        ValueError
            If the computed cost differs from the one in >>cost<<.
        """

        computed = cost.start_cost - avalanche_p * cost.cur_usage\
                                   - avalanche_h * cost.hist_usage\
                                   - avalanche_d * (cost.cur_usage - cost.prev_usage)

        if abs(computed - cost.derived_cost) / min(computed, cost.derived_cost) > 0.001:
            print "Derived cost incorrect."
            print cost
            print computed
            raise ValueError
        
        return computed
    #------------------------------------------------------------------------#

    for ind in sorted(cost_dict):
        costs = cost_dict[ind]

        """Test 1: Zero initial costs. In the first routing iteration, all costs
                   must be reset to zero, to give the router enough flexibility.
                   More generally, it can be some very small constant.
        """
        vpr_iter = 1
        try:
            vpr_iter = int(args.log.split('.', 1)[0].split('_')[-1])
        except:
            pass
 
        reset_cost = 1e-12 * (1.5 ** max(0, vpr_iter - 2))
        if abs(costs[0].stored_cost - reset_cost) > 0.1 * reset_cost:
            print "Cost differs from the reset cost in the first routing iteration."
            print costs[0]
            raise ValueError

        """Test 2: No change after routing. After the routing completes, before
                   the costs are finally printed, no change can occur. Because
                   the final update is only scheduled and not committed in the last
                   routing iteration, the exception can be in the stored costs.
                   For this, we first check if the scheduled value is correct and
                   if it is, commit the change. Then we remove the last cost in
                   the sequence, because it includes no historical cost updates
                   (routing is already complete) and would hence raise an exception
                   in that subsequent test.
        """
        computed = compute_cur_cost(costs[-2])
        costs[-2] = costs[-2]._replace(stored_cost = costs[-2].derived_cost)
        if costs[-1] != costs[-2]:
            print "Costs changed after routing end."
            print len(costs) - 1, costs[-1]
            print len(costs) - 2, costs[-2]
            raise ValueError
        costs = costs[:-1]

        hist_usage = costs[0].cur_usage
        for i, cost in enumerate(costs[1:], 1):

            """Test 3: Previous usage. Previous usage must be the current usage of the
                       preceding iteration.
            """
            if cost.prev_usage != costs[i - 1].cur_usage:
                print i, "Storing of previous usage incorrect."
                print costs[i-1]
                print cost
                raise ValueError
           
            """Test 4: Historical usage. Historical usage must be appropriately updated
                       from iteration to iteration.
            """ 
            if cost.hist_usage != hist_usage:
                print i, "Accumulation of historical usage incorrect."
                print cost
                print hist_usage
                raise ValueError
            hist_usage += cost.cur_usage
           
            """Test 5: Derived cost. The cost scheduled for update must be derivable from
                       the usage statistics and the avalanche parameters.
            """ 
            computed_cur_cost = compute_cur_cost(cost)
          
            """Test 6: Stored cost. Stored cost must equal the derived cost, when prints
                       come after the commits, which is the case here.
            """
            red_to_zero = cost.stored_cost == cost.derived_cost == 0.0
            if not red_to_zero and (abs(cost.stored_cost - cost.derived_cost) / min(cost.stored_cost, cost.derived_cost)\
                                    > avalanche_p):
                print "Stored cost and derived cost differ."
                print cost
                raise ValueError
            
            """Test 7: Current cost does not exceed the starting cost. Current cost can only
                       exceed the starting one when there are large swings in usage between iterations
                       and the derivative factor is high. With this undone, the cost must never exceed
                       the starting cost.
            """
            if (cost.stored_cost + avalanche_d * (cost.cur_usage - cost.prev_usage)) / cost.start_cost > 1.001:
                print "Stored cost with derivative factor undone still larger than the start cost."
                print cost
                print avalanche_d * (cost.cur_usage - cost.prev_usage)
                raise ValueError

    return True
##########################################################################

print check_costs(cost_dict)
