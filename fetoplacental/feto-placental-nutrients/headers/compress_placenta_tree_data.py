import numpy as np
import pandas as pd
import bz2
import pickle
import _pickle as cPickle
import bz2
from os import path

import placenta_tree_extract_process_results_updated as pl

TreeLabels = ['output-normal']
TreeFilenames = ['output-normal']

def process_and_compress_trees(node_fname, elem_fname, term_fname, rad_fname, \
                               press_fname, flow_fname, flux_vals,  Ntrees, \
                               paths_to_try, outfile_path):
    """Read in networks from placenta simulations and compress data."""
    flowbc_networks = np.zeros((len(flux_vals),Ntrees),dtype=object)
    flowbc_term_props = np.zeros((len(flux_vals),Ntrees),dtype=object)
    for n in np.arange(Ntrees):   #loop over trees
        for i in np.arange(len(flux_vals)):   #loop over fluxes
            if path.isdir(paths_to_try[i,n]):     #check it exists
                #read and extract properties
                flowbc_networks[i,n] = pl.read_network(paths_to_try[i,n] + node_fname, \
                            paths_to_try[i,n] + elem_fname, paths_to_try[i,n] + term_fname,\
                            paths_to_try[i,n] + rad_fname, paths_to_try[i,n] + press_fname,\
                            paths_to_try[i,n] + flow_fname)
                flowbc_term_props[i,n] = pl.get_terminal_properties(flowbc_networks[i,n])
            else:
                print(paths_to_try[i,n],' not found')

    with bz2.BZ2File(outfile_path,'w') as f:
        #write
        cPickle.dump([flux_vals,flowbc_networks], f)
    
    op_head, op_ext= path.splitext(outfile_path)
    outfile_path_term = op_head + '_terminal' + op_ext
    with bz2.BZ2File(outfile_path_term,'w') as f:
        #write term nodes file
        cPickle.dump([flowbc_term_props], f)

def compress_trees(PathToTrees, node_fname, elem_fname, term_fname, rad_fname, \
                         press_fname, flow_fname, flux_vals):

    tree_paths = np.array([TreeLabels[0]+'/output_flowbc_'])
    Ntrees = len(tree_paths)
    full_paths = np.zeros((len(flux_vals),Ntrees),dtype=object)
    for i in np.arange(len(flux_vals)):
        for j in np.arange(Ntrees):
            full_paths[i,j] = '%s/%s%.1fml_min_0.0/'%(PathToTrees,\
                                                tree_paths[j],flux_vals[i])

    process_and_compress_trees(node_fname, elem_fname, term_fname, rad_fname, \
                               press_fname, flow_fname, flux_vals, Ntrees, full_paths, \
                               'tree_flux_sweep_series.pbz2')

