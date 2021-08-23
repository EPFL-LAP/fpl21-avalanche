import argparse
import networkx as nx
from collections import namedtuple

get_attr = lambda attr, line : line.split("%s=\"" % attr, 1)[1].split('"', 1)[0]

max_x = max_y = -1
min_x = min_y = float("inf")

##########################################################################
def parse_switches(lines):
    """Parses the switch ids from an rr-graph.

    Parameters
    ----------
    lines : List[str]
        All lines of the rr-graph file.

    Returns
    -------
    Dict[str, int]
        A dictionary of switch ids, for all switch types.
    """

    switches = {}
    for line in lines:
        if "<switch " in line:
            name = get_attr("name", line)
            sw_id = int(get_attr("id", line))
            switches.update({name : sw_id})
        elif "</switches>" in line:
            break

    return switches
##########################################################################

##########################################################################
def parse_segments(lines):
    """Parses the segment ids from an rr-graph.

    Parameters
    ----------
    lines : List[str]
        All lines of the rr-graph file.

    Returns
    -------
    Dict[str, int]
        A dictionary of segment ids, for all segment types.
    """

    segmentes = {}
    for line in lines:
        if "<segment " in line:
            name = get_attr("name", line)
            sw_id = int(get_attr("id", line))
            segmentes.update({name : sw_id})
        elif "</segments>" in line:
            break

    return segmentes
##########################################################################

Coord = namedtuple("Coord", ['x', 'y'])
Wire = namedtuple("Wire", ["seg", "start", "end"])
##########################################################################
def parse_wires(lines):
    """Parses wire nodes from an rr-graph. 

    Parameters
    ----------
    lines : List[str]
        All lines of the rr-graph file.

    Returns
    -------
    nx.DiGraph
        A graph with wires instantiated as nodes.
    """

    global max_x
    global max_y
    global min_x
    global min_y

    G = nx.DiGraph()
    for lcnt, line in enumerate(lines[:-1]):
        if "<node " in line and "CHAN" in line:
            node = int(get_attr("id", line))
            inc = get_attr("direction", line) == "INC_DIR"
            low = Coord(int(get_attr("xlow", lines[lcnt + 1])), int(get_attr("ylow", lines[lcnt + 1])))
            high = Coord(int(get_attr("xhigh", lines[lcnt + 1])), int(get_attr("yhigh", lines[lcnt + 1])))
            seg = int(get_attr("segment_id", lines[lcnt + 3]))
            if inc:
                wire = Wire(seg, low, high)
            else:
                wire = Wire(seg, high, low)
            G.add_node(node, wire = wire)
            max_x = max(max_x, high.x)
            max_y = max(max_y, high.y)
            min_x = min(min_x, low.x)
            min_y = min(min_y, low.y)
           
        elif "</rr_nodes>" in line:
            break

    return G
##########################################################################

##########################################################################
def parse_edges(lines, G):
    """Parses wire nodes from an rr-graph. 

    Parameters
    ----------
    lines : List[str]
        All lines of the rr-graph file.
    G : nx.DiGraph
        A graph with wires instantiated as nodes.

    Returns
    -------
    None
    """

    for line in lines:
        if "<edge " in line:
            u = int(get_attr("src_node", line))
            if not G.has_node(u):
                continue
            v = int(get_attr("sink_node", line))
            if not G.has_node(v):
                continue
            sw_id = int(get_attr("switch_id", line))
            G.add_edge(u, v, switch = sw_id, persistance = "real")
##########################################################################

##########################################################################
def get_all_wires_ending_at_coord(G, inv_seg_dict, x, y):
    """Returns all wires that are ending at the particular coordinate.

    Parameters
    ----------
    G : nx.DiGraph
        A parsed routing graph (subgraph on wires).
    inv_seg_dict : Dict[int, str]
        Mapping from segment ids to names.
    x : int
        Hoizontal coordinate.
    y : int
        Vertical coordinate.
    """

    wire_types = {}
    for u, attrs in G.nodes(data = True):
        wire = attrs["wire"]
        x_u, y_u = wire.end
        if (x_u, y_u) != (x, y):
            continue
        wire_type = inv_seg_dict[wire.seg]
        if "potential_edge" in wire_type:
            continue
        
        if 'H' in wire_type:
            if wire.start < wire.end:
                d = 'R'
            else:
                d = 'L'
        elif wire.start < wire.end:
            d = 'U'
        else:
           d = 'D'
        wire_type += "_%s" % d

        try:
            wire_types[wire_type] += 1
        except:
            wire_types.update({wire_type : 1})

    for t in sorted(wire_types):
        print t, wire_types[t]
##########################################################################

##########################################################################
def flatten_potential(G, inv_seg_dict, inv_switch_dict):
    """Removes the edge-splitter nodes and adds actual edges between
    their parents and children, checking for consistency in the meantime.

    Parameters
    ----------
    G : nx.DiGraph
        A parsed routing graph (subgraph on wires).
    inv_seg_dict : Dict[int, str]
        Mapping from segment ids to names.
    inv_switch_dict : Dict[int, str]
        Mapping from switch ids to names.
           
    Returns
    -------
    Set[Coord]
        A set of coordinates that are influenced by boundary effects.

    Raises
    ------
    ValueError
        If the name of the edge-splitter does not match its neighbors.
        If the switches do not correspond to target nodes.
    """

    edge_splitters = [u for u, attrs in G.nodes(data = True) if "potential_edge" in inv_seg_dict[attrs["wire"].seg]]

    boundary_affected = set()
    for u in edge_splitters:
        u_seg = G.node[u]["wire"].seg
        x, y = G.node[u]["wire"].end
        fanin = list(G.pred[u])
        if len(fanin) != 1:
            understood_specificity = False
            if len(fanin) == 0:
                d_src = inv_seg_dict[u_seg].split("__")[1].split('_')[1]
                d_sink = inv_seg_dict[u_seg].split("__")[-1].split('_')[1]
                if len(G[u]) == 0:
                    #TODO: This can be made more rigorous, but in principle, the isolated egde-splitters
                    #should be the effect of reflections.
                    if x in (0, 1, max_x) or y in (0, 1, max_y):
                        understood_specificity = True
                    
                elif x <= 1 and d_src == 'R':
                    #No right-going horizontal wire can end to the right of a block at x = 0 (see VPR rendering).
                    understood_specificity = True
                elif x >= max_x and d_src == 'L':
                    #No left-going horizontal wire can end to the right of a block at x = max_x (see VPR rendering).
                    understood_specificity = True
                elif y <= 1 and d_src == 'U':
                    #There are only horizontal tracks at y = 0, so nothing going up can end at y = 1 (see VPR rendering).
                    understood_specificity = True
                elif y >= max_y and d_src == 'D':
                    #There are only horizontal tracks at y = grid_h - 1, so nothing going up can end at y = grid_h - 2
                    #(see VPR rendering).
                    understood_specificity = True

            if len(G[u]) > 1:
                understood_specificity = False
                #We must not tolerate multiple children under any circumstance.
            if understood_specificity:
                boundary_affected.add(Coord(x, y))
                G.remove_node(u)
                continue
            print "Node %d (%s) has %d parents." % (u, str(G.node[u]), len(fanin))
            print inv_seg_dict[u_seg], G.node[u]["wire"].end
            raise ValueError

        p = fanin[0]
        in_switch = G[p][u]["switch"]
        if inv_switch_dict[in_switch] != inv_seg_dict[u_seg]:
            print "In switch mismatch for edge (%d, %d)." % (p, u)
            raise ValueError

        tentative_p_label = inv_seg_dict[u_seg].split("__")[1]
        tentative_p_wire = tentative_p_label.split('_')[0]
        if tentative_p_wire.startswith('V'):
            tentative_p_wire += "_tap_0"
        p_seg = inv_seg_dict[G.node[p]["wire"].seg]
        if p_seg != tentative_p_wire:
            print "Edge-splitter/parent name mismatch for edge (%d, %d)." % (p, u)
            print p_seg, tentative_p_wire
            raise ValueError

        fanout = list(G[u])
        if len(fanout) != 1:
            print "Node %d (%s) has %d children." % (u, str(G.node[u]), len(fanout))
            raise ValueError

        c = fanout[0]
        c_seg = inv_seg_dict[G.node[c]["wire"].seg]
        out_switch = G[u][c]["switch"]
        if inv_switch_dict[out_switch] != c_seg:
            print "Out switch mismatch for edge (%d, %d)." % (u, c)
            raise ValueError

        tentative_c_label = inv_seg_dict[u_seg].split("__")[-1]
        tentative_c_wire = tentative_c_label.split('_')[1]
        if tentative_c_wire.startswith('V'):
            tentative_c_wire += "_tap_0"
        if c_seg != tentative_c_wire:
            print "Edge-splitter/child name mismatch for edge (%d, %d)." % (u, c)
            print c_seg, tentative_c_wire
            raise ValueError

        G.remove_node(u)
        ble_span = inv_seg_dict[u_seg].split("__")[-1].split('_')[0]
        G.add_edge(p, c, persistance = "potential", ble_span = ble_span)

    return boundary_affected
##########################################################################
       
##########################################################################
def get_switches(G, inv_seg_dict, x, y):
    """Returns a list of switches for a given location on the grid.
    This can be used for checking correctness of pattern assembly.
    - All patterns (at least sufficiently far from the edge of the FPGA)
      must be have an identical signature, which must in turn correspond
      to the data stored in >>stored_edges.save<<.

      Only the data about the cross-ble twists for accepted switches is
      not available in the routing graph at the moment. In principle, this
      can be added too, in the graph generation, by expanding the wire type set.

    Parameters
    ----------
    G : nx.DiGraph
        A parsed routing graph (subgraph on wires).
    inv_seg_dict : Dict[int, str]
        Mapping from segment ids to names.
    x : int
        x-coordinate.
    y : int
        y-coordinate.
    """

    nodes = [u for u in G if G.node[u]["wire"].end == (x, y)]

    include_ble_span = False
    #Since we do not have this for accepted switches, better not to introduce clutter anywhere.
    #TODO: Change this if the RR-graph is expanded.
    edges = {}
    for u in nodes:
        u_name = inv_seg_dict[G.node[u]["wire"].seg]
        for c in G[u]:
            c_name = inv_seg_dict[G.node[c]["wire"].seg]
            e = "%s_edge__%s__%s%s" % (G[u][c]["persistance"], u_name,\
                                      ("%s_" % G[u][c]["ble_span"]) if include_ble_span and G[u][c]["persistance"] == "potential"\
                                      else '', c_name)
            try:
                edges[e] += 1
            except:
                edges.update({e : 1})
    
       
    return edges
##########################################################################

##########################################################################
def canonical_switch_pattern_print(s):
    """Canonically prints a parsed switch pattern.

    Parameters
    ----------
    s : Dict[str, int]
        A parsed switch pattern.
    
    Returns
    -------
    str
        Canonical print out.
    """

    txt = ""
    for e in sorted(s):
        txt += "%s %d\n" % (e, s[e])

    return txt[:-1]
##########################################################################

##########################################################################
def calculate_maximum_potential_edge_counts(channel_composition, N, max_ble_span):
    """Computes the maximum number of possible occurrences per potential edge type.

    Parameters
    ----------
    channel_composition : Dict[str, int]
        Channel composition description.
    N : int
        Number of BLEs in the cluster.
    max_ble_span : int
        Maximum BLE span in the pattern.
    
    Returns
    -------
    Dict[str, int]
        Maximum number of occurrences of each edge type. 
    """

    back_dir = {'L' : 'R', 'R' : 'L', 'U' : 'D', 'D' : 'U'}

    counts = {}
    for src_ble in range(0, N):
        for sink_ble in range(max(0, src_ble - max_ble_span),\
                              min(N - 1, src_ble + max_ble_span) + 1):
            for w_src in channel_composition:
                src_dirs = ('L', 'R')
                if w_src[0] == 'V':
                    src_dirs = ('U', 'D')
                for src_dir in src_dirs:
                    for w_sink in channel_composition:
                        sink_dirs = ('L', 'R')
                        if w_sink[0] == 'V':
                            sink_dirs = ('U', 'D')
                        for sink_dir in sink_dirs:
                            if sink_dir == back_dir[src_dir]:
                                continue
                            inc = channel_composition[w_src] * channel_composition[w_sink]
                            try:
                                counts[(w_src, w_sink)] += inc 
                            except:
                                counts.update({(w_src, w_sink) : inc})

    e_str = lambda e : "potential_edge__%s%s__%s%s"\
          % (e[0], "_tap_0" if e[0][0] == 'V' else '',\
             e[1], "_tap_0" if e[1][0] == 'V' else '')

    return {e_str(e) : counts[e] for e in counts}
##########################################################################

##########################################################################
def get_accepted_switch_counts(store_file, N):
    """Counts the accepted switches per abstract type, as recorded in the
    store file.

    Parameters
    ----------
    store_file : str
        Filename of the file storing the edge adoption statistics.
        None indicates that no switches have been accepted until now.
    N : int
        Cluster size.

    Returns
    -------
    Dict[str, int]
        Number of used switches, per used abstract type.
    """

    if store_file is None:
        return {}

    with open(store_file, "r") as inf:
        lines = inf.readlines()

    used_edges = []
    for line in lines:
        if line.startswith('~'):
            break
        used_edges.append(line.strip())

    #------------------------------------------------------------------------#
    def strip_edge(edge):
        """Strips the BLE offset, direction, and index from the edge.

        Parameters
        ----------
        edge : str
            Edge description.

        Returns
        -------
        str
            Edge description, with the second-order information removed.
        int
            BLE offset.
        """

        ble_offset = int(edge.split("__")[-1].split("lut")[1].split('_')[0][1:])\
                   * (-1 if edge.split("__")[-1].startswith("lutm") else 1)

        u = edge.split("__")[1].split('_')[0]
        if u.startswith('V'):
            u += "_tap_0"
        v = edge.split("__")[2].split('_')[1]
        if v.startswith('V'):
            v += "_tap_0"
        
        base_edge = "potential_edge__%s__%s" % (u, v)

        return base_edge, ble_offset
    #------------------------------------------------------------------------#
    
    counts = {}
    for edge in used_edges:
        base_edge, ble_offset = strip_edge(edge)
        count = 0
        for src_ble in range(0, N):
            sink_ble = src_ble + ble_offset
            if sink_ble >= N or sink_ble < 0:
                continue
            count += 1
        try:
            counts[base_edge] += count
        except:
            counts.update({base_edge : count})

    return counts
##########################################################################

##########################################################################
def get_expected_switch_pattern(channel_composition, N, max_ble_span, store_file, is_final_iter = False):
    """Returns the expected switch pattern.

    Parameters
    ----------
    channel_composition : Dict[str, int]
        Channel composition description.
    N : int
        Number of BLEs in the cluster.
    max_ble_span : int
        Maximum BLE span in the pattern.
    store_file : str
        Filename of the file storing the edge adoption statistics.
    is_final_iter : Optional[bool], default = False
        Tells that the rr-graph was produced in the final iteration, hence there are no potential edges.

    Returns
    -------
    str
        A textual description of the switch pattern.
    """

    pattern = {}
    maximum_potential_switch_counts = calculate_maximum_potential_edge_counts(channel_composition, N, max_ble_span)
    accepted_switch_counts = get_accepted_switch_counts(store_file, N)

    for e in accepted_switch_counts:
        maximum_potential_switch_counts[e] -= accepted_switch_counts[e]
        pattern.update({e.replace("potential", "real") : accepted_switch_counts[e]})

    if not is_final_iter:
        pattern.update(maximum_potential_switch_counts)

    return canonical_switch_pattern_print(pattern)
##########################################################################  

##########################################################################  
def check_rr_all(rr_graph_file, store_file):
    """The top function performing all the implemented checks.

    Parameters
    ----------
    rr_graph_file : str
        Filename of the rr-graph.
    store_file : str
        Filename of the file storing the edge adoption statistics.

    Returns
    -------
    None

    Raises
    ------
    ValueError
        If any of the checks fails.
    """
    
    #Architecture generation parameters:
    #########################################
    N = 8
    MAX_BLE_SPAN = 3
    channel_composition = {"H1" : 2, "H2" : 1, "H4" : 1, "H6" : 1, "V1" : 2, "V4" : 1}

    with open(rr_graph_file, "r") as inf:
        lines = inf.readlines()
    #TODO: Read this in from the other files.
    #########################################

    switch_dict = parse_switches(lines)
    inv_switch_dict = {i : s for s, i in switch_dict.items()}
    seg_dict = parse_segments(lines)
    inv_seg_dict = {i : s for s, i in seg_dict.items()}
    
    G = parse_wires(lines)
    parse_edges(lines, G)

    total_graph_size = G.number_of_nodes()    
    boundary_affected = flatten_potential(G, inv_seg_dict, inv_switch_dict)
    real_graph_size = G.number_of_nodes()
    is_final_iter = total_graph_size == real_graph_size
    
    complete_switches = {}
    #Switch-blocks not influenced by the boundary effects.
    partial_switches = {}
    #Switch-blocks influenced by the boundary effects.
    
    for x in range(min_x, max_x + 1):
        for y in range(min_y, max_y + 1):
            switches = get_switches(G, inv_seg_dict, x, y)
            loc = Coord(x, y)
            if loc in boundary_affected or x in (min_x, max_x, 1) or y in (min_y, max_y, 1):
                #For x = 1 V y = 1, reflections happen on D and L, as there are no vertical tracks in y = 0
                #nor horizontal in x = 0.
                partial_switches.update({loc : switches})
            else:
                complete_switches.update({loc : switches})
 
    complete_switch = canonical_switch_pattern_print(complete_switches[Coord(3, 3)])#list(complete_switches)[0]])
    expected_switch = get_expected_switch_pattern(channel_composition, N, MAX_BLE_SPAN, store_file, is_final_iter = is_final_iter)

    if expected_switch != complete_switch:
        print "Expected and observed switches differ!"
        print "Expected:"
        print expected_switch
        print "Observed:"
        print complete_switch
        raise ValueError
 
    if not all(canonical_switch_pattern_print(s) == complete_switch for s in complete_switches.values()):
        print  "Nonboundary switches differ!"
 
        different_complete_switches = {}
        for loc in complete_switches:
            switch_txt = canonical_switch_pattern_print(complete_switches[loc])
            try:
                different_complete_switches[switch_txt].append(loc)
            except:
                different_complete_switches.update({switch_txt : [loc]})
          
        for switch in different_complete_switches:
            print sorted(different_complete_switches[switch])
            print switch

        raise ValueError

    print "Routing resource graph consistency check passed."
##########################################################################  
