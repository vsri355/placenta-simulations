from julia import Main
Main.include("./scripts/calculate_villous_fluxes.jl")
import numpy as np

#To run this script, the files tree_flux_sweep_series_terminal.pbz2 are needed

TreeFilenames = ["normal"]
outputdir = 'output-normal/'
inputdir = outputdir

flux_vals = np.array([20,40,60,80,100,150,200,250,300,350,400,450,500],dtype=float)

term_props = Main.load_pickle(decompress(inputdir+"tree_flux_sweep_series_terminal.pbz2"))


Nparallel = 3        #Number of convolute units in parallel
Nseries = 1          #Number of terminal villi in a row from a single mature intermediate villous
Nparallel_cap = 1    #Number of parallel capillaries in an imaged convolute (leiser)
Nconv = 10       #as per leiser 10 terminal conduits in a single feeding vessel
Ngens = 3

Main.calculate_and_save_allfluxes(outputdir,
                    term_props[0], TreeFilenames, flux_vals,
                    Nparallel, Nseries, Nparallel_cap, Nconv, Ngens, use_full_res=True)