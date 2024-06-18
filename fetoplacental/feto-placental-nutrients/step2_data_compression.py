import sys
import numpy as np
sys.path.append('./headers')
import compress_placenta_tree_data as cmpl

#This script converts UoA format outputs to UoM format outputs

PathToTrees = './'

flux_vals = np.array([20,40,60,80,100,150,200,250,300,350,400,450,500],dtype=float)
node_fname = 'full_tree.exnode'
elem_fname = 'full_tree.exelem'
term_fname = 'terminal.exnode'
rad_fname = 'radius_perf.exelem'
press_fname = 'pressure_perf.exnode'
flow_fname = 'flow_perf.exelem'

cmpl.compress_trees(PathToTrees, node_fname, elem_fname, term_fname, rad_fname, \
                      press_fname, flow_fname, flux_vals)