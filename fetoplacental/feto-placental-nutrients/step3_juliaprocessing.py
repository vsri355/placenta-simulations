import pyjulia
include("./headers/calculate_villous_fluxes.jl")

import _pickle as cPickle
import pickle
import bz2

#To run this script, the files tree_flux_sweep_series_terminal.pbz2 are needed

TreeFilenames = ["healthy_gen", "FGR_gen", "uct_tree1", "uct_tree2"]

flux_vals = Vector{Float64}([20,40,60,80,100,150,200,250,300,350,400,450,500]);

term_props = load_pickle(decompress("./output-compressed_data/tree_flux_sweep_series_terminal.pbz2"))

Nparallel = 3        #Number of convolute units in parallel
Nseries = 1          #Number of terminal villi in a row from a single mature intermediate villous
Nparallel_cap = 1    #Number of parallel capillaries in an imaged convolute (leiser)
Nconv = 10       #as per leiser 10 terminal conduits in a single feeding vessel
Ngens = 3

calculate_and_save_allfluxes("../Compressed_data/Series_trees/",
                    term_props[1], SerialTreeFilenames, flux_vals,
                    Nparallel, Nseries, Nparallel_cap, Nconv, Ngens; use_full_res=true)