"""Functions for reading in files from Auckland format and storing them
   in a network object"""
import numpy as np

def export_nodes_to_Manc_fmt(node_array, path):
    #write node array into format for C++ code
    node_file = open(path + '.nodes','w')
    for n in node_array:
        node_file.write('%d,%f,%f,%f\n'%(int(n[0]),n[1],n[2],n[3]))
    node_file.close()

def export_elems_to_Manc_fmt(elem_array, path, radius):
    #write edge array into format for C++ code
    elem_file = open(path + '.branches','w')
    i = 0
    for e in elem_array:
        elem_file.write('%d,%d,%d,%f\n,'%(e[0],e[1],e[2],radius[i]))
        i += 1
    elem_file.close()

def export_terminals_to_Manc_fmt(tnode_array, path, pressure=None):
    #write terminal array into format for C++ code
    i = 0
    tnode_file = open(path + '.termnodes','w')
    for n in tnode_array:
        tnode_file.write('%d'%int(n[0]))
        if pressure != None:
            tnode_file.write(',\'g\',%f'%pressure[i])
        i += 1
    tnode_file.close()

def read_network(exnode, exelem, term_exnode, radius=None, pressure=None, flow=None):
    #read network from ABI format

    #read in nodes
    nfile = open(exnode)
    nsr = nfile.readlines()

    #each node is a 3d position vector, guess initial number of nodes by filesize
    nodes = np.zeros((3,int(len(nsr)/4)))
    i = 0
    k = 0
    while i < len(nsr):
        if 'Node:' in nsr[i]:
            #read in x,y,z parts
            for j in np.arange(3):
                nodes[j,k] = float(nsr[i+j+1])
            i += 3
            k += 1
        i += 1
    #delete node entries not used
    nodes = np.delete(nodes,np.arange(k,int(len(nsr)/4)),1)
    nfile.close()

    #read in edge file
    efile = open(exelem)
    esr = efile.readlines()
    elems = np.zeros((2,int(len(esr)/4)),dtype=int)
    i = 0
    k = 0
    while i < len(esr):
        #parse nodes in and out for edge
        if 'Nodes:' in esr[i]:
            i += 1
            nn = esr[i].split()
            elems[0,k] = int(nn[0])-1  #node no.s start at 1 in file
            elems[1,k] = int(nn[1])-1
            k += 1
        i += 1
    elems = np.delete(elems,np.arange(k,int(len(esr)/4)),1)
    efile.close()

    #read in termnodes -- each entry is a node representing a terminal villous
    tfile = open(term_exnode)
    tsr = tfile.readlines()
    term_nodes = np.zeros(int(len(tsr)/4),dtype=int)
    i = 0
    k = 0
    while i < len(tsr):
        if 'Node:' in tsr[i]:
            ls = tsr[i].split()
            term_nodes[k] = int(ls[1])-1   #node no.s start at 1 in file
            i += 3
            k += 1
        i += 1
    term_nodes = np.delete(term_nodes,np.arange(k,int(len(tsr)/4)),0)
    tfile.close()

    #find edges in and out of each node
    edges_in = np.zeros(nodes.shape[1],dtype=object)
    edges_out = np.zeros(nodes.shape[1],dtype=object)
    for i in np.arange(len(edges_in)):
        edges_in[i] = np.array([],dtype=int)
        edges_out[i] = np.array([],dtype=int)
    for i in np.arange(elems.shape[1]):
        node = elems[1,i]
        edges_in[node] = np.append(edges_in[node],i)
        node = elems[0,i]
        edges_out[node] = np.append(edges_out[node],i)

    #count number of edges in an out of each node
    n_edges_in = np.zeros(nodes.shape[1],dtype=int)
    n_edges_out = np.zeros(nodes.shape[1],dtype=int)
    for i in np.arange(len(n_edges_in)):
        n_edges_in[i] = len(edges_in[i])
        n_edges_out[i] = len(edges_out[i])

    nn = np.arange(nodes.shape[1])
    #collect entry (those with no edges in) and exit nodes (those with no edges out)
    entry_nodes = nn[n_edges_in==0]
    exit_nodes = nn[n_edges_out==0]

    #read in (edge) radii -- if file has been given
    elem_rads = np.zeros(elems.shape[1])
    if radius != None:
        rfile = open(radius)
        rsr = rfile.readlines()
        i = 0
        k = 0
        while i < len(rsr):
            if 'Values:' in rsr[i]:
                i += 1
                rr = rsr[i].split()
                elem_rads[k] = 0.5*(float(rr[0]) + float(rr[1]))
                k += 1
            i += 1
        rfile.close()

    #read in (node) pressures -- if file has been given
    node_press = np.zeros(nodes.shape[1])
    if pressure != None:
        pfile = open(pressure)
        psr = pfile.readlines()
        i = 0
        k = 0
        while i < len(psr):
            if 'Node:' in psr[i]:
                i += 1
                node_press[k] = float(psr[i])
                k += 1
            i += 1
        pfile.close()

    #read in (edge) flow -- if file has been given
    elem_flows = np.zeros(elems.shape[1])
    if flow != None:
        ffile = open(flow)
        fsr = ffile.readlines()
        i = 0
        k = 0
        while i < len(fsr):
            if 'Values:' in fsr[i]:
                i += 1
                ff = fsr[i].split()
                elem_flows[k] = 0.5*(float(ff[0]) + float(ff[1]))
                k += 1
            i += 1
        ffile.close()


    #return all tree information in the form of a dict
    return {"Nodes": nodes, "Edges": elems, "Term_nodes": term_nodes,\
            "Edges_in": edges_in, "Edges_out": edges_out, "N_edges_in": n_edges_in,\
            "N_edges_out": n_edges_out, "Entry_nodes": entry_nodes, \
            "Exit_nodes": exit_nodes, "Radius": elem_rads, "Pressure": node_press,\
            "Flows": elem_flows}



def get_terminal_properties(network):
    netsplit = split_trees(network)
    Pin = network["Pressure"][network["Term_nodes"]]
    Pout = network["Pressure"][network["Term_nodes_ven"]]
    Pdrop = Pin - Pout
    TermFlows = network["Flows"][network["Term_edges"]]
    TermEdgesA = np.zeros(len(network["Term_nodes"]),dtype=int)
    TermEdgesV = np.zeros(len(network["Term_nodes_ven"]),dtype=int)
    i = 0
    for t in network["Term_nodes"]:
        TermEdgesA[i] = network["Edges_in"][t][0]
        i += 1

    i = 0
    for t in network["Term_nodes_ven"]:
        TermEdgesV[i] = network["Edges_out"][t][0]
        i += 1

    TermRadsA = network["Radius"][TermEdgesA]
    TermRadsV = network["Radius"][TermEdgesV]

    return Pdrop, TermFlows, TermRadsA, TermRadsV, Pin, Pout


def split_trees(network):
    network["Term_edges"] = np.zeros(len(network["Term_nodes"]),dtype=int)
    i = 0
    for t in network["Term_nodes"]:
        network["Term_edges"][i] = network["Edges_out"][t][0]
        i += 1
    network["Term_nodes_ven"] = network["Edges"][1,network["Term_edges"]]

    return network

def separate_trees(network):
    pass
