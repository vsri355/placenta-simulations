import sys
sys.path.append('./scripts')
import placenta_calculations as pcalc
import compress_placenta_tree_data as cmpl
import placenta_sim_plots as psp
import pickle
import numpy as np

outputdir = 'output-normal/'
TreeLabels = ['normal']
InputFluxFiles = [ outputdir+ 'tv_fluxes_normal_flow250.pkl']

Nparallel = 1        #Number of convolute units in parallel
Nseries = 3          #Number of terminal villi in a row from a single mature intermediate villous
Nparallel_cap = 1    #Number of parallel capillaries in an imaged convolute (leiser), typically 1
Nconv = 10       #as per leiser 10 terminal conduits in a single feeding vessel (typically 10)
Ngens = 3   #Typically 3

psp.plot_flux_histograms(outputdir+'/terminal_flow.pdf', InputFluxFiles,TreeLabels,
                         Nseries, Nparallel)

dataset = {}
with open(outputdir+ 'flux_sweep_outcomes.pkl','rb') as f:
    dataset = pickle.load(f)


psp.plot_flow_vs_absolute_uptake(outputdir+'/flow_vs_absolute_uptake_oxygen.pdf',\
                            dataset,TreeLabels,[(0.0,40.0)], SolutesToOmit=[1,2,3,4,5,6,7],
                            plot_ideal=True)