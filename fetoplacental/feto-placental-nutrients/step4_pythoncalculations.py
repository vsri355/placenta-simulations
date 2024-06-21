import sys
sys.path.append('./scripts')
import placenta_calculations as pcalc
import shutil
import glob
import numpy as np

outputdir = './output-normal'

tree_path = outputdir+'/tree_flux_sweep_series.pbz2'
term_path = outputdir+'/tree_flux_sweep_series_terminal.pbz2'

villous_flux_dir = outputdir
Nparallel = 3        #Number of convolute units in parallel
Nseries = 1          #Number of terminal villi in a row from a single mature intermediate villous
Nparallel_cap = 1    #Number of parallel capillaries in an imaged convolute (leiser)
Nconv = 10       #as per leiser 10 terminal conduits in a single feeding vessel
Ngens = 3
vtk_print_fluxes = np.array([250])    #flux value to produce vtk file for
pcalc.calculate_and_save_measures(outputdir+"/series_flux_sweep_outcomes.pkl", \
        tree_path, term_path, villous_flux_dir, ['normal'], ['normal'], \
        pcalc.SoluteNames, pcalc.DC, pcalc.B, pcalc.Dt, Nparallel, Nseries, Nparallel_cap, Nconv, Ngens, \
        vtk_print_fluxes)


for file in glob.glob(r'./*.vtk'):
    shutil.move(file,outputdir)