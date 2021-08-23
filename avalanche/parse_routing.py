"""Parses the final routing and converts it into a Networkx graph.
"""

import copy
import networkx as nx

##########################################################################
def parse_global_routing(filename):
    """Parses .route global routing files.

    Parameters
    ----------
    filename : str
        Name of the .route file.

    Returns
    -------
    Dict[str, nx.DiGraph]
        A dictionary of route trees, indexed by the net names.
    """

    

    with open(filename, "r") as inf:
        lines = inf.readlines()
    
    trees = {}
    net = None
    prev_switch = -1
    for lcnt, line in enumerate(lines):
        if line.startswith("Net"):
            #if trees.get(net, None) is not None:
                #for e in sorted(trees[net].edges(), key = lambda e : trees[net][e[0]][e[1]]["ecnt"]):
                #    print e
                #raw_input()
            net = line.split('(')[1].split(')')[0]
            tree = nx.DiGraph()
            trees.update({net : tree})
            prev_node = None
            ecnt = 0
        elif line.startswith("Node:"):
            node_num = int(line.split()[1])
            node_type = line.split()[2]
            node_loc = (int(line.split('(')[1].split(',')[0]), int(line.split(',')[1].split(')')[0]))
            node_pin = line.split()[-3] if node_type in ("IPIN", "OPIN") else None
            node_switch = int(line.split()[-1])
            if not tree.has_node(node_num):
                tree.add_node(node_num, node_type = node_type,\
                                        node_loc = node_loc,\
                                        node_pin = node_pin,\
                                        node_switch = prev_switch)
                if prev_node == None and node_type != "SOURCE":
                    print "No previous node and the current node is not a source."
                    print lcnt, line
                    raise ValueError
                if node_type != "SOURCE":
                    tree.add_edge(prev_node, node_num, ecnt = ecnt)
                    ecnt += 1
            prev_node = node_num
            prev_switch = node_switch
            #FIXME: Check this again!

    return trees
##########################################################################

##########################################################################
def parse_local_routing(filename):
    """Parses local routing information from a .net packing file.

    Parameters
    ----------
    filename : str
        Name of the packing file.

    Returns
    -------
    nx.DiGraph
        Graph representing the local routing.
    nx.DiGraph
        Block tree.
    """

    get_attr = lambda attr, line : line.split("%s=\"" % attr)[1].split('"')[0]

    with open(filename, "r") as inf:
        lines = inf.readlines()

    block_tree = nx.DiGraph()
    block_stack = []
    lvl = 0
    parent_dict = {}
    for lcnt, line in enumerate(lines):
        if "<block name=" in line:
            name = get_attr("name", line)
            instance = get_attr("instance", line)
            mode = None
            start = lcnt
            try:
                mode = get_attr("mode", line)
            except:
                pass
            lvl += 1
            block_stack.append((name, instance, mode, start, lvl))
            parent_dict.update({lvl : len(block_stack) - 1})
            if "/>" in line:
                name, instance, mode, start, lvl = block_stack.pop(-1)
                block_tree.add_node(start, name = name, instance = instance, mode = mode,\
                                    start = start, end = lcnt, lvl = lvl, parent  = None)
                lvl -= 1
                if lvl > 0:
                    parent = block_stack[parent_dict[lvl]]
                    block_tree.node[start]["parent"] = parent[3] #start
        elif "</block>" in line:
            name, instance, mode, start, lvl = block_stack.pop(-1)
            block_tree.add_node(start, name = name, instance = instance, mode = mode,\
                                start = start, end = lcnt, lvl = lvl, parent  = None)
            lvl -= 1
            if lvl > 0:
                parent = block_stack[parent_dict[lvl]]
                block_tree.node[start]["parent"] = parent[3] #start

    for u, attrs in block_tree.nodes(data = True):
        if attrs["parent"] is not None:
            assert attrs["parent"] in block_tree, attrs["parent"]
            block_tree.add_edge(attrs["parent"], u)

    #TODO: Check what happens with the LUTs in wire mode.
    
    netlist = nx.DiGraph()
    for b, attrs in block_tree.nodes(data = True):
        start = attrs["start"]
        end = attrs["end"]
        state = "idle"
        for line in lines[start:end]:
            if "<inputs>" in line:
                state = "rd_inputs"
                continue
            if "</inputs>" in line:
                state = "idle"
                continue
            if "</outputs>" in line:
                #Parsing of outputs of the top-level model will not be initiated, but we do not care for that in any case.
                state = "idle"
                break
            if "<outputs>" in line:
                state = "rd_outputs"
                continue
            if state == "idle":
                continue
            if state in ("rd_inputs", "rd_outputs"):
                if "<port " in line:
                    name = get_attr("name", line)
                    port_map = line.split('>', 1)[1].rsplit('<', 1)[0].split()
                    for p, pin in enumerate(port_map):
                        if pin == "open":
                            driver, through = None, None
                        elif (state == "rd_inputs" and attrs["lvl"] < 3)\
                             or (state == "rd_outputs" and block_tree.out_degree(b) == 0):
                            driver, through = pin, None
                        else:
                            driver, through = pin.split("-&gt;")
                        netlist.add_node("%d.%s[%d]" % (b, name, p), block = b,\
                                         node_type = state.split('_')[1][:-1], driver = driver, through = through)

    for u, attrs in netlist.nodes(data = True):
        if attrs["node_type"] == "input":
            if attrs["through"] is not None:
                parent = list(block_tree.pred[attrs["block"]])[0]
                driver_instance = attrs["driver"].split('.')[0]
                if driver_instance == block_tree.node[parent]["instance"].split('[')[0]:
                    driver = "%d.%s" % (parent, attrs["driver"].split('.')[1])
                else:
                    for sibling in block_tree[parent]:
                        if block_tree.node[sibling]["instance"] == driver_instance:
                            driver = "%d.%s" % (sibling, attrs["driver"].split('.')[1])
                assert netlist.has_node(driver), driver
                netlist.add_edge(driver, u)
        else:
            if attrs["through"] is not None:
                driver_instance = attrs["driver"].split('.')[0]
                for child in block_tree[attrs["block"]]:
                    if block_tree.node[child]["instance"] == driver_instance:
                        driver = "%d.%s" % (child, attrs["driver"].split('.')[1])
                assert netlist.has_node(driver), driver
                netlist.add_edge(driver, u)
                #print driver, u, "out"
                #print attrs["through"]

                       
    #Now annotate the CLB outputs with the primitives that drive them, so that merging of global routing is possible.
    nodes_per_block = {}
    for u, attrs in netlist.nodes(data = True):
        try:
            nodes_per_block[attrs["block"]].add(u)
        except:
            nodes_per_block.update({attrs["block"] : set([u])})

    for u, attrs in netlist.nodes(data = True):
        #IO-pads are treated specially.
        if "inpad" in u:
            if attrs["through"] is None and attrs["driver"] is not None:
                netlist.node[u]["primitive"] = attrs["driver"]
        if u.split('.', 1)[1].split('[')[0] != 'O':
            continue
        if attrs["driver"] is None:
            continue
        stack = [attrs["driver"]]
        block = attrs["block"]
        while stack:
            driver = stack.pop()
            found = False
            for child_block in block_tree[block]:
                if block_tree.node[child_block]["instance"] == driver.split('.', 1)[0]:
                    found = True
                    break
            assert found, "Child block containing the instance of the driver %s not found in the block tree." % driver
            block = child_block
            for v in nodes_per_block[block]:
                if v.split('.', 1)[1] == driver.split('.', 1)[1]:
                    driver = netlist.node[v]["driver"]
                    if block_tree.out_degree[block] == 0:
                        netlist.node[u]["primitive"] = driver
                        print "%s <- %s" % (u, driver)
                    else:
                        stack.append(driver)
                    break 

    #Now annotate the placement locations to be able to merge the global routing.
    with open(filename.replace(".net", ".place"), "r") as inf:
        lines = inf.readlines()

    block_name_dict = {attrs["name"] : u for u, attrs in block_tree.nodes(data = True) if attrs["lvl"] == 2}
    for line in lines[5:]:
        name = line.split()[0]
        x = int(line.split()[1])
        y = int(line.split()[2])
        block = block_name_dict[name]
        stack = [block]
        while stack:
            block = stack.pop(-1)
            block_tree.node[block]["loc"] = (x, y)
            for child in block_tree[block]:
                stack.append(child)

    for u, attrs in netlist.nodes(data = True):
        block = attrs["block"]
        if block_tree.node[block]["lvl"] == 1:
            continue
        netlist.node[u]["loc"] = block_tree.node[block]["loc"]
            
    return netlist, block_tree        
##########################################################################

##########################################################################
def merge_global_and_local_routing(global_routing_trees, local_routing_netlist):
    """Merges the global routing trees into the local routing netlist.

    Parameters
    ----------
    global_routing_trees : Dict[str, nx.DiGraph]
        A dictionary of routing trees, one for each net.
    local_routing_netlist : nx.DiGraph
        A local routing netlist.

    Returns
    -------
    nx.DiGraph
        A single merged netlist.
    """

    netlist = copy.deepcopy(local_routing_netlist)
    #TODO: Connect also the primary IOs.
    #NOTE: Probably stale TODO. Check this!

    #Group local routing input pins by the driving net
    ipins = {}
    for u, attrs in local_routing_netlist.nodes(data = True):
        #Skip the top module.
        if u.split('.', 1)[0] == '1':
            continue
        if u.split('.', 1)[1].split('[', 1)[0] in ('I', "outpad"):
            driver = attrs["driver"]
            if driver is None:
                continue
            try:
                ipins[driver].add(u)
            except:
                ipins.update({driver : set([u])})

    #Group local routing output pins by the driving net
    opins = {}
    for u, attrs in local_routing_netlist.nodes(data = True):
        #Skip the top module.
        if u.split('.', 1)[0] == '1':
            continue
        if u.split('.', 1)[1].split('[', 1)[0] in ('O', "inpad") :
            driver = attrs.get("primitive", None)
            if driver is None:
                continue
            try:
                opins[driver].add(u)
            except:
                opins.update({driver : set([u])})

    for net in global_routing_trees:
        tree = global_routing_trees[net]
        for u, attrs in tree.nodes(data = True):
            merged_node = "global__%d" % u 
            netlist.add_node(merged_node)
            netlist.node[merged_node].update(copy.deepcopy(attrs))
            loc = attrs["node_loc"]
            if attrs["node_type"] == "IPIN":
                for i in ipins[net]:
                    #print net, i
                    if loc == local_routing_netlist.node[i]["loc"]:
                        netlist.add_edge(merged_node, i)
            elif attrs["node_type"] == "OPIN":
                for o in opins[net]:
                    if loc == local_routing_netlist.node[o]["loc"]:
                        netlist.add_edge(o, merged_node)
        for u, v, attrs in tree.edges(data = True):
            merged_u = "global__%d" % u 
            merged_v = "global__%d" % v 
            netlist.add_edge(merged_u, merged_v)
            netlist[merged_u][merged_v].update(attrs)

    #Remove the sink and the source nodes as they are irrelevant and create issues in
    #total slack computation.

    rm_list = []
    for u, attrs in netlist.nodes(data = True):
        if u.startswith("global__") and attrs["node_type"] in ("SOURCE", "SINK"):
            rm_list.append(u)

    for u in rm_list:
        netlist.remove_node(u)

    return netlist
##########################################################################

##########################################################################
def annotate_delays(netlist, block_tree, arc_filename, rr_filename):
    """Parses the delays from the architecture file and annotates the netlist accordingly.

    Parameters
    ----------
    netlist : nx.DiGraph
        Parsed routing graph (both local and global merged).
    block_tree : nx.DiGraph
        Tree of parsed packed blocks.
    arc_filename : str
        Name of the architecture file used to produce the routing.
    rr_filename : str
        Name of the rr-graph file. Necessary for matching the switch types.

    Returns
    -------
    nx.DiGraph
        Delay-annotated netlist.
    """

    get_attr = lambda attr, line : line.split("%s=\"" % attr)[1].split('"')[0]

    with open(rr_filename, "r") as inf:
        lines = inf.readlines()

    switch_delays = {}
    rd_delays = False
    for line in lines:
        if "</switches>" in line:
            break
        if "<switch id" in line:
            switch_id = int(get_attr("id", line))
            rd_delays = True
            continue
        if rd_delays:
            td = float(get_attr("Tdel", line))
            switch_delays.update({switch_id : td})
            rd_delays = False
    switch_delays.update({-1 : 0.0})

    with open(arc_filename, "r") as inf:
        lines = inf.readlines()

    local_delays = {}
    for lcnt, line in enumerate(lines):
        if "delay_constant" in line:
            in_port = get_attr("in_port", line)
            in_port = in_port.split('[')[0].split('.', 1)[0]
            out_port = get_attr("out_port", line)
            out_port = out_port.split('[')[0].split('.', 1)[0]
            local_delays.update({(in_port, out_port) : float(get_attr("max", line))})
        elif "<delay_matrix" in line:
            local_delays.update({"lut" : float(lines[lcnt + 1].strip())})
        elif "<T_setup" in line:
            local_delays.update({"tsu" : float(get_attr("value", line))})
        elif "<T_clock_to_Q" in line:
            local_delays.update({"tclkq" : float(get_attr("max", line))})    

    for u, attrs in netlist.nodes(data = True):
        if attrs.get("node_switch", None) is not None:
            netlist.node[u]["td"] = switch_delays[attrs["node_switch"]]
        block = attrs.get("block", None)
        if block is None:
            continue
        target_type = block_tree.node[block]["instance"].split('[', 1)[0]
        if target_type == "lut" and u.split('.', 1)[1].startswith("out"):
            netlist.node[u]["td"] = local_delays["lut"]
            continue
        if target_type == "ff":
            if u.split('.', 1)[1][0] == 'D':
                netlist.node[u]["td"] = local_delays["tsu"]
                netlist.node[u]["FF_D"] = True
            elif u.split('.', 1)[1][0] == 'Q':
                netlist.node[u]["td"] = local_delays["tclkq"]
                netlist.node[u]["FF_Q"] = True
            continue
        driver = attrs.get("driver", None)
        if driver is None:
            continue
        driver = driver.split('[', 1)[0].split('.', 1)[0]
        td = local_delays.get((driver, target_type), 0.0)
        netlist.node[u]["td"] = td

    return netlist
##########################################################################

##########################################################################
def insert_primitive_edges(netlist, block_tree):
    """Inserts edges through primitives, to enable correct timing analysis
    (FFs are opened later).

    Parameters
    ----------
    netlist : nx.DiGraph
        Timing-annotated routing netlist.
    block_tree : nx.DiGraph
        Tree of parsed packed blocks.

    Returns
    -------
    nx.DiGraph
        Netlist updated with the cross-primitive edges.
    """

    ins = {}
    outs = {}
    for u in netlist:
        if "global" in u:
            continue
        instance = int(u.split('.', 1)[0])
        if block_tree.out_degree(instance):
            continue
        if ".in" in u:
            try:
                ins[instance].add(u)
            except:
                ins.update({instance : set([u])})
        elif ".out" in u:
            try:
                outs[instance].add(u)
            except:
                outs.update({instance : set([u])})

    for instance in ins:
        for i in ins[instance]:
            if netlist.node[i]["driver"] is None:
                continue
            for o in outs.get(instance, []):
                if netlist.node[o]["driver"] is None:
                    continue
                netlist.add_edge(i, o)
##########################################################################

##########################################################################
def strip_unused(netlist):
    """Strips unused nodes.

    Parameters
    ----------
    netlist : nx.DiGraph
        Timing-annotated routing netlist.

    Returns
    -------
    None
    """

    rm_list = []
    for u, attrs in netlist.nodes(data = True):
        if attrs.get("driver", "drv") is None:
            rm_list.append(u)
        elif attrs.get("node_type", "") == "output" and netlist.out_degree(u) == 0:
            rm_list.append(u)

    for u in rm_list:
        netlist.remove_node(u)
##########################################################################

##########################################################################
def sta(netlist):
    """Performs a simple static timing analysis.

    Parameters
    ----------
    netlist : nx.DiGraph
        Timing-annotated routing netlist.
    
    Returns
    -------
    float
        Critical path delay.
    nx.DiGraph
        Netlist copy with annotated slacks on the edges.
    """

    #------------------------------------------------------------------------#
    def open_ffs(netlist):
        """Opens FFs.

        Parameters
        ----------
        netlist : nx.DiGraph
            Timing-annotated routing netlist.

        Returns
        -------
        nx.DiGraph
            A copy of the netlist, with FFs opened.
        """

        opened = copy.deepcopy(netlist)
        for u, v in netlist.edges():
            if netlist.node[u].get("FF_D", False) and netlist.node[v].get("FF_Q", False):
                opened.remove_edge(u, v)

        return opened
    #------------------------------------------------------------------------#

    #------------------------------------------------------------------------#
    def split_clock_domains(netlist):
        """Splits the netlist into different clock domains.
        This must be done when wafers are being used.

        Parameters
        ----------
        netlist : nx.DiGraph
            Timing-annotated routing netlist.

        Returns
        -------
        List[nx.DiGraph]
            A collection of graphs that represent each clock domain.
            It may be that some of them are still disconnected.
        """

        underlying = netlist.to_undirected()
        components = list(sorted(nx.connected_components(underlying), key = len))

        clock_domains = {}
        for i, c in enumerate(components):
            for u in c:
                primitive = netlist.node[u].get("primitive", "")
                if "circ" in primitive:
                    words = primitive.split('_')
                    domain = None
                    for j, w in enumerate(words):
                        if w == "circ":
                            domain = '_'.join([w, words[j + 1]])
                            break
                    try:
                        clock_domains[domain].append(i)
                    except:
                        clock_domains.update({domain : [i]})
                    break
                driver = netlist.node[u].get("driver", "")
                if "circ" in driver:
                    words = driver.split('_')
                    domain = None
                    for j, w in enumerate(words):
                        if w == "circ":
                            domain = '_'.join([w, words[j + 1]])
                            break
                    try:
                        clock_domains[domain].append(i)
                    except:
                        clock_domains.update({domain : [i]})
                    break

        domain_graphs = []
        for domain in clock_domains:
            nodes = set()
            for c in clock_domains[domain]:
                nodes |= components[c]
            domain_graphs.append(copy.deepcopy(netlist))
            rm_list = [u for u in netlist if not u in nodes]
            #NOTE: This can probably be done more efficiently, but this is robust.
            for u in rm_list:
                domain_graphs[-1].remove_node(u)

        return domain_graphs            
    #------------------------------------------------------------------------#

    netlist_cp = open_ffs(netlist)

    domain_graphs = split_clock_domains(netlist_cp)

    TNS = 0.0
    WNS = 0.0
    for dg in domain_graphs:
        topo_nodes = list(nx.algorithms.dag.topological_sort(dg))
        rev_topo_nodes = reversed(topo_nodes)
    
        cp = nx.algorithms.dag.dag_longest_path(dg)
        for u in cp:
            print u, netlist_cp.node[u]
    
        for u in topo_nodes:
            td = netlist_cp.node[u].get("td", 0)
            tar = td + max([0] + [netlist_cp.node[p]["tar"] for p in netlist_cp.pred[u]])
            netlist_cp.node[u]["tar"] = tar
            dg.node[u]["tar"] = tar
    
        cpd = max([attrs["tar"] for u, attrs in dg.nodes(data = True)])
        tns = sum([attrs["tar"] for u, attrs in dg.nodes(data = True) if dg.out_degree(u) == 0])
    
        TNS += tns
        WNS = max(WNS, cpd)
    
        for u in rev_topo_nodes:
            td = netlist_cp.node[u].get("td", 0)
            treq = min([cpd] + [netlist_cp.node[c]["treq"] for c in netlist_cp[u]])
            netlist_cp.node[u]["treq"] = treq
            netlist_cp.node[u]["slack"] = treq - netlist_cp.node[u]["tar"]
            assert treq >= netlist_cp.node[u]["tar"], "Negative slack on node %s %s" % (str(u), str(netlist_cp.node[u]))
    
        for u, v in dg.edges():
            netlist_cp[u][v]["slack"] = netlist_cp.node[v]["treq"] - netlist_cp.node[u]["tar"]

    print "WNS =", WNS
    print "TNS =", TNS

    return cpd, netlist_cp
##########################################################################

##########################################################################
def get_setup_histogram(netlist, vpr_log):
    """Returns a setup histogram, with bucket thresholds read from a VPR log.
    Useful for verification.

    Parameters
    ----------
    netlist : nx.DiGraph
        Timing-annotated routing netlist.
    vpr_log : str
        File name of the VPR log from which the thresholds are to be read.
    """

    with open(vpr_log, "r") as inf:
        lines = inf.readlines()

    thresholds = set()
    rd = False
    for line in lines:
        if "Setup slack histogram:" in line:
            rd = True
            continue
        if not rd:
            continue
        if not line.startswith('['):
            break
        up = float(line.split()[1][1:-1])
        down = float(line.split()[2][1:-1])
        thresholds.add(up)
        thresholds.add(down)

    thresholds = list(sorted(thresholds))

    sink_tars = sorted([attrs["tar"] for u, attrs in netlist.nodes(data = True) if netlist.out_degree(u) == 0])
    print sink_tars
    histogram = []
    for t, thr in enumerate(thresholds[1:], 1):
        bucket = len([s for s in sink_tars if thresholds[t - 1] <= s <= thr])
        histogram.append(bucket)

    print thresholds
    print histogram

    for u, attrs in netlist.nodes(data = True):
        if netlist.out_degree(u) == 0:
            print u, attrs
            raw_input()
##########################################################################

##########################################################################
def get_wire_type_criticalities(netlist):
    """Calculates the criticalities of all wire types.

    Parameters
    ----------
    netlist : nx.DiGraph
        Timing-annotated routing netlist.

    Returns
    -------
    Dict[str, float]
        A criticality dictionary.
    """

    slacks = {}
    cpd = max([attrs["tar"] for u, attrs in netlist.nodes(data = True) if netlist.out_degree(u) == 0])

    for u, attrs in netlist.nodes(data = True):
        if not "node_switch" in attrs:
            continue
        switch = attrs["node_switch"]
        slacks.update({switch : min(slacks.get(switch, float('inf')), attrs["slack"])})

    crits = {s : 1 - slacks[s] / cpd for s in slacks}
    
    print slacks        
    print crits

    return crits
##########################################################################

##########################################################################
def fetch_timing_data():
    """Top function for parsing the timing data of a circuit.

    Parameters
    ----------
    None

    Returns
    -------
    Dict[str, float]
        A dictionary of criticalities for each wire type.
    """

    with open("vpr_stdout.log", "r") as inf:
        lines = inf.readlines()


    #If max_criticality is 0.0, timing was not taken into account during routing, so
    #any criticalities would only create noise (this is the case for e.g., Gnl).
    #In that case, return None, to default to all 1.0 criticalities.

    for line in lines:
        if "--max_criticality" in line:
            max_criticality = float(line.split("--max_criticality", 1)[1].split()[0])
            if max_criticality == 0.0:
                return None

    rd = False
    for line in lines:
        if "VPR was run with the following command-line:" in line:
            rd = True
            continue
        if not rd:
            continue
        arc_name = line.split()[1].split(".xml", 1)[0]
        circ = line.split()[2].split(".blif")[0]
        break

    global_routing_trees = parse_global_routing("%s.route" % circ)
    local_routing_netlist, packing_block_tree = parse_local_routing("%s.net" % circ)
    netlist = merge_global_and_local_routing(global_routing_trees, local_routing_netlist)
    strip_unused(netlist)
    annotate_delays(netlist, packing_block_tree, "%s.xml" % arc_name, "%s_rr.xml" % arc_name)
    insert_primitive_edges(netlist, packing_block_tree)
    cpd, annotated_netlist = sta(netlist)

    return get_wire_type_criticalities(annotated_netlist)
#########################################################################
