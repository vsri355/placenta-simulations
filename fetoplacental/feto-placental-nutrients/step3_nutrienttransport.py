from julia import Main
Main.include("./scripts/calculate_villous_fluxes.jl")
import numpy as np
import sys
sys.path.append('./scripts')
import placenta_calculations as pcalc
#import shutil
import glob


#To run this script, the files tree_flux_sweep_series_terminal.pbz2 are needed

TreeFilenames = ["normal"]
outputdir = 'output-normal/'
inputdir = outputdir
flux_vals = np.array([20,40,60,80,100,150,200,250,300,350,400,450,500],dtype=float)
Nparallel = 1        #Number of convolute units in parallel
Nseries = 3          #Number of terminal villi in a row from a single mature intermediate villous
Nparallel_cap = 1    #Number of parallel capillaries in an imaged convolute (leiser), typically 1
Nconv = 10       #as per leiser 10 terminal conduits in a single feeding vessel (typically 10)
Ngens = 3   #Typically 3
tree_path = outputdir+'/tree_flux_sweep.pbz2'
term_path = outputdir+'/tree_flux_sweep_terminal.pbz2'
villous_flux_dir = outputdir

term_props = Main.load_pickle(decompress(term_path))

Main.calculate_and_save_allfluxes(outputdir,
                    term_props[0], TreeFilenames, flux_vals,
                    Nparallel, Nseries, Nparallel_cap, Nconv, Ngens, use_full_res=True)


vtk_print_fluxes = np.array([250])    #flux value to produce vtk file for
pcalc.calculate_and_save_measures(outputdir+"/flux_sweep_outcomes.pkl", \
        tree_path, term_path, villous_flux_dir, ['normal'], ['normal'], \
        pcalc.SoluteNames, pcalc.DC, pcalc.B, pcalc.Dt, Nparallel, Nseries, Nparallel_cap, Nconv, Ngens, \
        vtk_print_fluxes)


#for file in glob.glob(r'./*.vtk'):
#    shutil.move(file,outputdir)