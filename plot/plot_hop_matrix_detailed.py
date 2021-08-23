import os
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from itertools import izip

store = "inputs/stored_edges.save"
init_delay_file = "inputs/arc_delays_iter_1.log"
final_delay_file = "inputs/arc_delays.log"

N = 8
W = 30

tracks = ("H1_L_0",\
          "H1_L_1",\
          "H1_R_0",\
          "H1_R_1",\
          "H2_L_0",\
          "H2_R_0",\
          "H4_L_0",\
          "H4_R_0",\
          "H6_L_0",\
          "H6_R_0",\
          "V1_U_0",\
          "V1_D_0",\
          "V1_U_1",\
          "V1_D_1",\
          "V4_U_0",\
          "V4_D_0")

init_delays = {}
final_delays = {}

compose_node = lambda wire, ble, coords : '_'.join([wire] + [str(ble)] + [str(coords[0]), str(coords[1])])

get_len = lambda u : int(u.split('_', 1)[0][1:])
offset_dict = {'L' : (-1, 0), 'R' : (1, 0), 'U' : (0, 1), 'D' : (0, -1)}
get_offset = lambda u : offset_dict[u.split('_')[1]]
get_index = lambda u : int(u.split('_')[2])
get_ble = lambda u : int(u.split('_')[3])
get_coords = lambda u : (int(u.split('_')[4]), int(u.split('_')[5]))
get_wire = lambda u : '_'.join(u.split('_')[:3])
get_wire_type = lambda u : u.split('_')[0]
get_init_td = lambda u : init_delays[get_wire_type(u)]
get_final_td = lambda u : final_delays[get_wire_type(u)]

##########################################################################
def parse_delays(delay_file):
    """Parses wire delays from an architecture delay file.

    Parameters
    ----------
    delay_file : str
        Name of the delay file.
    
    Returns
    -------
    Dict[str, int]
        A dictionary of delays in picoseconds for all wire types.
    """

    delays = {}
    with open(delay_file, "r") as inf:
        lines = inf.readlines()

    for line in lines:
        if line.startswith("cb"):
            break
        wire = line.split()[0].replace("_tap_0", '')
        td = int(round(float(line.split()[1]) * 1e12))
        delays.update({wire : td})

    return delays
##########################################################################

##########################################################################
def parse_edges(edge_store):
    """Parses the stored edges from a file.

    Parameters
    ----------
    edge_store : str
        Name of the edge store file.

    Returns
    -------
    List[Tuple[str, int]]
        A list of parsed edges.
    """


    store_ble_offset_dict = {"lutm1" : -1, "lutp0" : 0, "lutp1" : 1}
    strip_store_line = lambda line : line.replace("potential_edge__", '').replace("_tap_0", '').strip()
    parse_store_line = lambda line : (line.split("__")[0], '_'.join(line.split("__")[1].split('_')[1:]),\
                                      store_ble_offset_dict[line.split("__")[1].split('_', 1)[0]])


    parsed_edges = []
    with open(edge_store, "r") as inf:
        lines = inf.readlines()
    
    for line in lines:
        if line[0] == '~':
            break
        parsed_edges.append(parse_store_line(strip_store_line(line)))

    return parsed_edges
##########################################################################

##########################################################################
def construct_graph(parsed_edges):
    """Constructs a graph on which distances can be measured.

    Parameters
    ----------
    parsed_edges : List[Tuple[str, int]]
        A list of parsed edges.

    Returns
    -------
    nx.DiGraph
        The constructed graph.
    """

    G = nx.DiGraph()
    
    nodes_per_coord = {}
    for x in range(0, W):
        for y in range(0, W):
            nodes_per_coord.update({(x, y) : set()})
            src = "src_%d_%d" % (x, y)
            G.add_node(src, final_td= 0.0, init_td = 0.0)
            nodes_per_coord[(x, y)].add(src)
            sink = "sink_%d_%d" % (x, y)
            G.add_node(sink, final_td = 0.0, init_td = 0.0)
            nodes_per_coord[(x, y)].add(sink)
            for ble in range(0, N):
                for t in tracks: 
                    t_node = compose_node(t, ble, (x, y))
                    G.add_node(t_node, final_td = get_final_td(t),\
                                       init_td = get_init_td(t))
                    G.add_edge(src, t_node, final_td = 0.0, init_td = 0.0)
                    nodes_per_coord[(x, y)].add(t_node)
    
    for x in range(0, W):
        for y in range(0, W):
            for node in nodes_per_coord[(x, y)]:
                if node[0] not in ('H', 'V'):
                    continue
                offset = get_offset(node)
                L = get_len(node)
                dest = (x + L * offset[0], y + L * offset[1])
                sink = "sink_%d_%d" % (dest[0], dest[1])
                if not G.has_node(sink):
                    continue
                G.add_edge(node, sink, final_td = G.node[node]["final_td"], init_td = G.node[node]["init_td"])
    
                src_ble = get_ble(node)
                src_wire = get_wire(node)
                for u, v, ble_offset in parsed_edges:
                    if u != src_wire:
                        continue
                    dest_node = compose_node(v, src_ble + ble_offset, dest)
                    if G.has_node(dest_node):
                        G.add_edge(node, dest_node, final_td = G.node[node]["final_td"], init_td = G.node[node]["init_td"])

    return G
##########################################################################

##########################################################################
def get_map(G, key = "hop"):
    """Produces a distance map.

    Parameters
    ----------
    G : nx.DiGraph
        The graph on which to seek a map.
    key : Optional[str], default = hop
        Specifies which map to produce (final_delay, init_delay, or hop).
    
    Returns
    -------
    Dict[Tuple[int], int]
        The assembled matrix.
    """
                  
    src = "src_%d_%d" % (W / 2, W / 2)

    if key == "hop":
        distances = nx.single_source_dijkstra_path_length(G, src)
        sub = 1
        #-1 for src to wire
    elif key == "final_td":
        distances = nx.single_source_dijkstra_path_length(G, src, weight = "final_td")
        sub = 0
    elif key == "init_td":
        distances = nx.single_source_dijkstra_path_length(G, src, weight = "init_td")
        sub = 0
    else:
        print "Unrecognized matrix type", key
        raise ValueError

    matrix = {}
    for y in range(0, W):
        for x in range(0, W):
            sink = "sink_%d_%d" % (x, y)
            try:
                matrix.update({(x, y) : distances[sink] - sub})
            except:
                matrix.update({(x, y) : -1})

    return matrix
##########################################################################

##########################################################################
def fill_graph(G):
    """Fills the graph with all missing allowed edges, which enables lower
    bound calculation.

    Parameters
    ----------
    G : nx.DiGraph
        The graph to fill.
    
    Returns
    -------
    nx.DiGraph
        Filled graph.
    """

    G = G.copy()

    for src_track in tracks:
        src_offset = get_offset(src_track)
        src_L = get_len(src_track)
        for dest_track in tracks:
            dest_offset = get_offset(dest_track)
            if dest_offset == (-1 * src_offset[0], -1 * src_offset[1]):
                continue
    
            for x in range(0, W):
                for y in range(0, W):
                    src_node = compose_node(src_track, 0, (x, y))
                    dest_node = compose_node(dest_track, 0, (x + src_L * src_offset[0], y + src_L * src_offset[1]))
                    if G.has_node(dest_node):
                        G.add_edge(src_node, dest_node, final_td = G.node[src_node]["final_td"], init_td = G.node[src_node]["init_td"])
    
    return G
##########################################################################

##########################################################################
def plot_matrix(matrix, min_val, max_val, title, include_bar, filename):
    """Plots the given matrix.

    Parameters
    ----------
    matrix : Dict[Tuple[int], int]
        Matrix to be plotted.
    min_val : int
        Minimum value to use. If None, minmum found in the data will be used.
    max_val : int
        Maximum value to use. If None, maximum found in the data will be used.
    title : str
        Title of the plot.
    include_bar : bool
        Include the color bar.
    filename : str
        Name of the file into which to save the plot.

    Returns
    -------
    None
    """

    min_val_observed = float('inf')
    max_val_observed = -1000000
    data = np.zeros(shape = (W, W), dtype = np.int)
    for key, val in matrix.items():
        data[key] = int(round(val))
        min_val_observed = min(min_val_observed, data[key])
        max_val_observed = max(max_val_observed, data[key])
    data = data.transpose(1, 0)
    #Pyplot rotates the axes.

    xlabel= "X-offset"
    ylabel="Y-offset"

    fontsize = 6
  
    if min_val is None:
        min_val = min_val_observed
    if max_val is None:
        max_val = max_val_observed 
         
    plt.figure(figsize = (W / 2, W / 2))
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    c = plt.pcolor(data, edgecolors = 'k', linewidths = 4, cmap = 'RdBu', vmin = min_val, vmax = max_val)

    #------------------------------------------------------------------------#
    def show_values(pc, fmt="%d", **kw):
        """Prints the values inside the heatmap cells.

        Parameters
        ----------
        pc : plt.pcolor
            Plotted heatmap.
        fmt : Optional[str], default = %d
            String formatter.
        **kw : void
            Other arguments to be passed to the text writter.
        
        Returns
        -------
        None
        """

        pc.update_scalarmappable()
        ax = pc.get_axes()
        for p, color, value in izip(pc.get_paths(), pc.get_facecolors(), pc.get_array()):
            x, y = p.vertices[:-2, :].mean(0)
            if np.all(color[:3] > 0.5):
                color = (0.0, 0.0, 0.0)
            else:
                color = (1.0, 1.0, 1.0)
            ax.text(x, y, fmt % value, ha = "center", va = "center", color=color, size = fontsize, **kw)
    #------------------------------------------------------------------------#

    show_values(c)

    if include_bar:
        plt.colorbar(c)

    plt.savefig(filename)
    plt.close()
##########################################################################

##########################################################################
def relative_matrix(matrix_a, matrix_b):
    """Creates a relative matrix from two given matrices ((a-b)/b*100%).

    Parameters
    ----------
    matrix_a : Dict[Tuple[int], int]
        New matrix.
    matrix_b : Dict[Tuple[int], int]
        Base matrix.

    Returns
    -------
    Dict[Tuple[int], int]
        Relative matrix of the two given ones.
    """

    return {coords : int(round(float((matrix_a[coords] - matrix_b[coords])) / matrix_b[coords] * 100)) for coords in matrix_b}
##########################################################################

init_delays = parse_delays(init_delay_file)
final_delays = parse_delays(final_delay_file)
parsed_edges = parse_edges(store)

G = construct_graph(parsed_edges)
hop_matrix = get_map(G, key = "hop")
delay_matrix = get_map(G, key = "final_td")

G_full = fill_graph(G)
lb_hop_matrix = get_map(G_full, key = "hop")
lb_delay_matrix = get_map(G_full, key = "init_td")

normalized_hop_matrix = relative_matrix(hop_matrix, lb_hop_matrix)
normalized_delay_matrix = relative_matrix(delay_matrix, lb_delay_matrix)

plot_matrix(hop_matrix, -1, W, "Hop-Distance Matrix", True, "hop_matrix_a_at_1.5.pdf")
plot_matrix(normalized_hop_matrix, 0, 300, "Lower Bound-Normalized Hop-Distance Matrix", True,\
            "hop_matrix_a_at_1.5_normalized.pdf")

plot_matrix(delay_matrix, 0, None, "Delay-Distance Matrix", True, "delay_matrix_a_at_1.5.pdf")
plot_matrix(normalized_delay_matrix, 0, None, "Lower Bound-Normalized Delay-Distance Matrix", True,\
            "delay_matrix_a_at_1.5_normalized.pdf")
