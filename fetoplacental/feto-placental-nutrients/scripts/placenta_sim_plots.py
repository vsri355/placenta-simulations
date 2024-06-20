import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import placenta_tree_extract_process_results_updated as pl
import placenta_calculations as pcalc
import compress_placenta_tree_data as cmpl
import scipy.stats as stats
from os import path
import os
import pickle

os.environ["PATH"] += os.pathsep + '/usr/local/texlive/2015/bin/x86_64-darwin'
plt.rcParams.update({'font.size': 15, 'text.usetex': False})
titles = ["(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)","(i)"]
styles = ['-o','-s','-x','-+','--o','--s','--x','--+',':o']
linestyles = ['-','--',':','-.']

def plot_grid(filename, x, y_mat, solute_numbers, tree_numbers, solute_names,
              tree_labels, ylims, xlabel, ylabel, y2_mat = np.zeros((0,0,0))):
    Nsolutes = len(solute_numbers)
    Nrows = int(1 + np.floor((Nsolutes-1)/3))
    Ncols = np.min([3,Nsolutes])
    fig, ax = plt.subplots(Nrows, Ncols, figsize=(4*Ncols,4*Nrows))
    for (isol, s) in enumerate(solute_numbers):
        #get current axis
        if Nrows > 1:
            ax_obj = ax[int(np.floor(isol/3)),isol%3]
        else:
            if Ncols > 1:
                ax_obj = ax[isol%3]
            else:
                ax_obj = ax
        
        for (itree, n) in enumerate(tree_numbers):
            if np.shape(x) == np.shape(y_mat):
                xh = x[:,n,s]
            else:
                xh = x
            yh = y_mat[:,n,s]
            x1h = xh[np.isnan(yh) == False]
            y1h = yh[np.isnan(yh) == False]

            if np.shape(y2_mat) == np.shape(y_mat):
                ax_obj.plot(x1h, y1h, styles[n], c=('C%d'%n), label=tree_labels[n])
                y2h = y2_mat[:,n,s]
                x2h = xh[np.isnan(y2h) == False]
                y2h = yh[np.isnan(y2h) == False]
                ax_obj.plot(x2h, y2h, styles[5+n], c=('C%d'%n))
            else:
                ax_obj.plot(x1h, y1h, styles[n], c=('C%d'%n), label=tree_labels[n])
        if isol == len(solute_numbers)-1:
            ax_obj.legend()
        ax_obj.set_ylim(ylims[isol])
        ax_obj.set_xlabel(xlabel)
        ax_obj.set_ylabel(ylabel)
        ax_obj.set_title('%s %s'%(titles[isol],solute_names[s]))
    fig.tight_layout()
    fig.savefig(filename)
    return fig, ax

def plot_flow_vs_absolute_uptake(filename, dataset, tree_names, ylims, SolutesToOmit=[],\
                                 plot_ideal = False):
    BoolSolutesToInclude = np.ones(len(dataset['Solutes']),dtype=bool)
    BoolSolutesToInclude[SolutesToOmit] = False
    SolutesToInclude = np.arange(len(dataset['Solutes']))[BoolSolutesToInclude]
    Nsolutes = len(SolutesToInclude)
    Ntrees = len(dataset['Tree Labels'])
    tree_numbers = np.arange(Ntrees)
    solute_names = np.array(pcalc.SoluteNames)
    xlabel = 'Flow (ml/min)'
    ylabel = 'Uptake $N$ (ml/min)'
    if plot_ideal:
        plot_grid(filename, dataset['Fluxes'], dataset['Total Exchange']*60, SolutesToInclude,\
              tree_numbers, solute_names, tree_names, ylims, xlabel, ylabel, dataset['Ideal Exchange']*60)
    else:
        plot_grid(filename, dataset['Fluxes'], dataset['Total Exchange']*60, SolutesToInclude,\
              tree_numbers, solute_names, tree_names, ylims, xlabel, ylabel)

def plot_flow_vs_relative_uptake(filename, dataset, tree_names, ylims, SolutesToOmit=[],\
                                 plot_ideal = False):
    BoolSolutesToInclude = np.ones(len(dataset['Solutes']),dtype=bool)
    BoolSolutesToInclude[SolutesToOmit] = False
    SolutesToInclude = np.arange(len(dataset['Solutes']))[BoolSolutesToInclude]
    Nsolutes = len(SolutesToInclude)
    Ntrees = len(dataset['Tree Labels'])
    tree_numbers = np.arange(Ntrees)
    solute_names = np.array(pcalc.SoluteNames)
    xlabel = 'Flow (ml/min)'
    ylabel = 'Relative uptake $N/N_{\\rm max}$'
    if plot_ideal:
        plot_grid(filename, dataset['Fluxes'], dataset['Total Exchange']/dataset['Max Exchange'], \
            SolutesToInclude, tree_numbers, solute_names, tree_names, ylims, xlabel, ylabel,\
            dataset['Ideal Exchange']/dataset['Max Exchange'])
    else:
        plot_grid(filename, dataset['Fluxes'], dataset['Total Exchange']/dataset['Max Exchange'], \
            SolutesToInclude, tree_numbers, solute_names, tree_names, ylims, xlabel, ylabel)

def plot_flow_vs_scaled_uptake(filename, dataset, tree_names, ylims, Neffs, SolutesToOmit=[],\
                                 plot_ideal = False):
    BoolSolutesToInclude = np.ones(len(dataset['Solutes']),dtype=bool)
    BoolSolutesToInclude[SolutesToOmit] = False
    SolutesToInclude = np.arange(len(dataset['Solutes']))[BoolSolutesToInclude]
    Nsolutes = len(SolutesToInclude)
    Ntrees = len(dataset['Tree Labels'])
    tree_numbers = np.arange(Ntrees)
    solute_names = np.array(pcalc.SoluteNames)
    xlabel = 'Flow per functional unit (nl/min)'
    ylabel = 'Uptake per functional unit (nl/min)'
    nflows = np.shape(dataset['Total Exchange'])[0]
    nsols = np.shape(dataset['Total Exchange'])[2]
    scaled_exchange = np.zeros(np.shape(dataset['Total Exchange']))
    scaled_fluxes = np.zeros(np.shape(dataset['Total Exchange']))
    for i in np.arange(nflows):
        for k in np.arange(nsols):
            scaled_exchange[i,:,k] = dataset['Total Exchange'][i,:,k]/Neffs
            scaled_fluxes[i,:,k] = dataset['Fluxes'][i]/Neffs
    plot_grid(filename, scaled_fluxes*1E6, scaled_exchange*60*1E6, SolutesToInclude, \
                  tree_numbers, solute_names, tree_names, ylims, xlabel, ylabel)

def plot_invDa_vs_absolute_uptake(filename, dataset, tree_names, ylims, SolutesToOmit=[],\
                                 plot_ideal = False):
    BoolSolutesToInclude = np.ones(len(dataset['Solutes']),dtype=bool)
    BoolSolutesToInclude[SolutesToOmit] = False
    SolutesToInclude = np.arange(len(dataset['Solutes']))[BoolSolutesToInclude]
    Nsolutes = len(SolutesToInclude)
    Ntrees = len(dataset['Tree Labels'])
    tree_numbers = np.arange(Ntrees)
    solute_names = np.array(pcalc.SoluteNames)
    xlabel = 'Inverse Global Damk\\xf6hler $\\bar{\\mathrm{Da}}^{-1}$'
    ylabel = 'Uptake $N$ (ml/min)'
    if plot_ideal:
        plot_grid(filename, dataset['Inverse Damkohler Global'], \
            dataset['Total Exchange']*60, SolutesToInclude, tree_numbers, \
            solute_names, tree_names, ylims, xlabel, ylabel, dataset['Ideal Exchange']*60)
    else:
        plot_grid(filename, dataset['Inverse Damkohler Global'], \
            dataset['Total Exchange']*60, SolutesToInclude, tree_numbers, \
            solute_names, tree_names, ylims, xlabel, ylabel)
        
def plot_invDa_vs_relative_uptake(filename, dataset, tree_names, ylims, SolutesToOmit=[],\
                                 plot_ideal = False):
    BoolSolutesToInclude = np.ones(len(dataset['Solutes']),dtype=bool)
    BoolSolutesToInclude[SolutesToOmit] = False
    SolutesToInclude = np.arange(len(dataset['Solutes']))[BoolSolutesToInclude]
    Nsolutes = len(SolutesToInclude)
    Ntrees = len(dataset['Tree Labels'])
    tree_numbers = np.arange(Ntrees)
    solute_names = np.array(pcalc.SoluteNames)
    xlabel = 'Inverse Global Damk\\xf6hler $\\bar{\\mathrm{Da}}^{-1}$'
    ylabel = 'Relative uptake $N/N_{\\rm max}$'
    fig, ax = plot_grid(filename, dataset['Inverse Damkohler Global'], \
            dataset['Total Exchange']/dataset['Max Exchange'], SolutesToInclude, tree_numbers, \
            solute_names, tree_names, ylims, xlabel, ylabel)
    if plot_ideal:
        for (i,s) in enumerate(SolutesToInclude):
            xlim = ax[i].get_xlim()
            InvDa = xlim[0] + (np.arange(101)/100)*(xlim[1] - xlim[0])
            RelEx = np.zeros(len(InvDa))
            Nmax = pcalc.Dt[s] * pcalc.L_vil_default * pcalc.DC[s] * Nseries
            InvDa2Q = 1E9 * pcalc.Dt[s] * pcalc.L_vil_default / (pcalc.B[s] * Nparallel)
            for j in np.arange(len(InvDa)):
                Qh = InvDa[j] * InvDa2Q
                Nh =  pcalc.Net_flux_calc(InvDa[j], Qh, pcalc.Dt[s],\
                            pcalc.B[s], pcalc.Dt[s], pcalc.DC[s], Nseries, Nparallel)
                RelEx[j] = Nh/Nmax
            ax[i].plot(InvDa,RelEx,':',c='k',label='Ideal')
        ax[len(SolutesToInclude)-1].legend()
        fig.tight_layout()
        fig.savefig(filename)

def plot_uptake_vs_L(Filename, dataset, SolutesToOmit=[], plot_ideal=False):
    BoolSolutesToInclude = np.ones(len(dataset['Solutes']),dtype=bool)
    BoolSolutesToInclude[SolutesToOmit] = False
    SolutesToInclude = np.arange(len(dataset['Solutes']))[BoolSolutesToInclude]
    Nsolutes = len(SolutesToInclude)
    Ntrees = len(dataset['Tree Labels'])
    tree_numbers = np.arange(Ntrees)
    solute_names = np.array(pcalc.SoluteNames)[SolutesToInclude]
    xlabel = 'Inverse Damk\\xf6hler Da$^{-1}$'
    ylabel = 'Relative uptake $N/N_{\\rm max}$'
    L_vals = dataset["LValues"]
    fig,ax = plt.subplots(1,2,figsize=(10,4))
    min_invda = 1E6
    max_invda = 0
    for it in tree_numbers:
        ax[0].plot(L_vals, dataset['Total Exchange'][:,it,Solutes[0]]*60, linestyles[it], \
                     label=dataset['Tree Labels'][it], c = 'k') 
        for (isol,s) in enumerate(SolutesToInclude):
            ax[0].plot(L_vals, dataset['Total Exchange'][:,it,s]*60, linestyles[it], \
                     label='', c = pcalc.SoluteColours[s])

    for it in tree_numbers:
        for (isol,s) in enumerate(SolutesToInclude):
            invDa = dataset['Inverse Damkohler Global'][:,it,s]*pcalc.L_vil_default/L_vals
            if np.min(invDa) < min_invda:
                min_invda = np.min(invDa)
            if np.max(invDa) > max_invda:
                max_invda = np.max(invDa)
            ax[1].plot(invDa, dataset['Total Exchange'][:,it,s]/dataset['Max Exchange'][:,it,s], \
                     linestyles[it], c = pcalc.SoluteColours[s])
            if it == 0:
                ax[1].plot(invDa, dataset['Total Exchange'][:,it,s]/dataset['Max Exchange'][:,it,s], \
                     linestyles[it], label=pcalc.SoluteNames[s], c = pcalc.SoluteColours[s])
            else:
                ax[1].plot(invDa, dataset['Total Exchange'][:,it,s]/dataset['Max Exchange'][:,it,s], \
                     linestyles[it],label='',c = pcalc.SoluteColours[s])

 
    ax[0].set_xlabel('Diffusion lengthscale $\\mathcal{L} (mm)$') 
    ax[0].set_ylabel('Uptake $N$ (ml/min)') 
    ax[0].legend() 
    ax[0].set_yscale('log') 
    ax[0].set_xlim(L_vals[0],L_vals[-1])
    ax[1].set_xlabel('Inverse Damk\\xf6hler Da$^{-1}$') 
    ax[1].set_ylabel('Relative uptake $N/N_{\\rm max}$')
    ax[1].set_xscale('log') 
    ax[1].set_xlim(min_invda,max_invda)
    ax[1].set_ylim(0,1) 
    ax[1].legend() 

    fig.tight_layout()
    fig.savefig(Filename)

def plot_flux_histograms(Filename, InputFluxFiles, Labels, Nseries, Nparallel):
    fig, ax = plt.subplots(figsize=(6,4.25))
    max_flux = 180 #in nl/min
    bh = np.arange(0,max_flux,2)
    
    for n in np.arange(len(InputFluxFiles)):
        FlowDict = pickle.load(open(InputFluxFiles[n],'rb'))
        ax.hist(FlowDict["TV_fluxes"]*60000, bins=bh, label=Labels[n], color='C%d'%n,\
                alpha=0.3, density=True, weights=FlowDict["Nweights"])
        kernel = stats.gaussian_kde(FlowDict["TV_fluxes"]*60000, weights=FlowDict["Nweights"])
        ax.plot(bh,kernel(bh),c='C%d'%n)
        
    ax.set_xlabel('Terminal villus flux (nl min$^{-1}$)')
    ax.set_ylabel('Probability density')
    
    ax.set_xlim((0,  max_flux))
    ax.set_ylim((0,  0.1))
    xticks = np.arange(0,max_flux+1,20)
    ax.set_xticks(xticks)
    ax.legend()
    #create inv Da axis too
    ax2 = ax.twiny()
    maxinvDa = pcalc.inv_damkohler(max_flux/60000,pcalc.Dt_O2[0], pcalc.B_O2[0], Nseries, Nparallel)
    #round downto nearest multiple of 2
    maxinvDaround = np.floor(maxinvDa/2)*2
    invDa_ticks = np.arange(0,maxinvDaround+1,2)
    invDa_flux_ticks = (invDa_ticks/maxinvDa)*max_flux
    def tick_function(flux_nl_min):
        flux = flux_nl_min/60000
        invDa = pcalc.inv_damkohler(flux, pcalc.Dt_O2[0], pcalc.B_O2[0], Nseries, Nparallel)
        #FLux in mm3/s
        return ["%.2f" % z for z in invDa]
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xticks(invDa_flux_ticks)
    ax2.set_xticklabels(tick_function(invDa_flux_ticks))
    ax2.set_xlabel("O$_2$ Inverse Damk" + "\xf6" + "hler (Da$^{-1}$)")
    #ax.set_ylim(0.5,10000)
    #ax.set_yscale('log')
    #ax.set_title('(a)')

    fig.tight_layout()
    fig.savefig(Filename)

def plot_relative_uptake_histograms(Filename, TVFluxFilePaths, Solutes, TreeLabels, Nseries, Nparallel):
    Ntrees = len(TVFluxFilePaths)
    Nrows = int(1 + np.floor((Ntrees-1)/2))
    Ncols = np.min([2,Ntrees])
    fig, ax = plt.subplots(Nrows,Ncols,figsize=(Ncols*5,Nrows*3.5))
    bh = np.arange(0,1.001,0.01)
    for n in np.arange(Ntrees):
        bottom_row = True
        left_col = True
        if Nrows > 1:
                ax_obj = ax[int(n/2),n%2]
                if int(n/2) < Nrows-1:
                    bottom_row = False
        else:
            if Ncols > 1:
                ax_obj = ax[n]
            else:
                ax_obj = ax
        FlowDict = pickle.load(open(TVFluxFilePaths[n],'rb'))
        for (isol, s) in enumerate(Solutes):
            Nf_spec = pcalc.Net_flux(FlowDict["TV_fluxes"], pcalc.Dt[s], pcalc.B[s], pcalc.Dt[s], \
                                  pcalc.DC[s], Nseries, Nparallel)
            Nmax = pcalc.Max_diff_cap(pcalc.Dt[s], pcalc.DC[s], Nseries, Nparallel)
            ax_obj.hist(Nf_spec/Nmax, bins=bh, label=pcalc.SoluteNames[s], alpha=0.5, \
                color=pcalc.SoluteColours[s], density=True, weights=FlowDict["Nweights"])
        
        ax_obj.set_title('%s %s'%(titles[n],TreeLabels[n]))
        ax_obj.set_ylim((0,20))
        ax_obj.set_xlim((0,1))
        if bottom_row:
            ax_obj.set_xlabel('$N/N_{\\rm max}$')
        if (n == 0 & isol == 0):
            ax_obj.legend()
        if n%2 > 0:
            ax_obj.set_ylabel('Probability density') 
    fig.tight_layout()
    fig.savefig(Filename)

def plot_varying_meso_scale_heterogeneity(filename, dataset, sdvals, SolutesToOmit=[], 
                                          plot_ideal=True):
    BoolSolutesToInclude = np.ones(len(dataset['Solutes']),dtype=bool)
    BoolSolutesToInclude[SolutesToOmit] = False
    SolutesToInclude = np.arange(len(dataset['Solutes']))[BoolSolutesToInclude]
    print(SolutesToInclude)
    NSolutes = len(SolutesToInclude)
    fig,ax = plt.subplots(1,NSolutes,figsize=(12,4))
    for (ns,s) in enumerate(SolutesToInclude):
        IDaGh = dataset['Inverse Damkohler Global'][:,:,s]
        Ntot_flowbc = dataset['Total Exchange'][:,:,s]
        Nmaxtot_flowbc = dataset['Max Exchange'][:,:,s]
        Ideal_Ntot = dataset['Ideal Exchange'][:,:,s]
        for n in np.arange(len(sdvals)):
            #skip odd numbers
            if(n%2==0):
                Label = '$\\sigma = $%s$\\times$'%(sdvals[n])
                ax[ns].plot(IDaGh[:,n], Ntot_flowbc[:,n]/Nmaxtot_flowbc[:,n], styles[int(n/2)],label=Label)
        ax[ns].plot(IDaGh[:,0],Ideal_Ntot[:,0]/ Nmaxtot_flowbc[:,0],'--',c='k'%n,label='Ideal')
        ax[ns].set_xlim(0,np.max(IDaGh))
        ax[ns].set_ylim(0,1)
        ax[ns].set_xlabel('Da$^{-1}$')
        ax[ns].set_title('%s %s'%(titles[ns],dataset['Solutes'][s]))
    ax[0].set_ylabel('$N/N_{max}$')
    ax[2].set_ylim(0,0.1)
    ax[2].legend()
    fig.tight_layout()
    fig.savefig(filename)

def plot_series_vs_parallel(OutFilename, series_summary_filepath, parallel_summary_filepath):
    sfdata = pickle.load(open(series_summary_filepath,'rb'))
    pfdata = pickle.load(open(parallel_summary_filepath,'rb'))
    #plot series vs parallel for all trees O2 OR  one tree all solutes?
    SolutesToInclude = [0]
    Nsolutes = len(SolutesToInclude)
    Ntrees = len(sfdata['Tree Labels'])
    fig,ax = plt.subplots(1,figsize=(6,5))
    styles = ['-s','--o',':+','-.x']

    for i in np.arange(Ntrees):
        for s in SolutesToInclude:
            yh = sfdata['Total Exchange'][:,i,s]*60
            xh = sfdata['Fluxes']
            x1h = xh[np.isnan(yh) == False]
            y1h = yh[np.isnan(yh) == False]
            if(i==0):
                ax.plot(x1h,y1h,'-',c='C0',label='Series')
            ax.plot(x1h,y1h,styles[i],c='C0')

            yh = pfdata['Total Exchange'][:,i,s]*60
            xh = pfdata['Fluxes']
            x1h = xh[np.isnan(yh) == False]
            y1h = yh[np.isnan(yh) == False]
            if(i==0):
                ax.plot(x1h,y1h,'-',c='C1',label='Parallel')
            ax.plot(x1h,y1h,styles[i],c='k',label=sfdata['Tree Labels'][i])
            ax.plot(x1h,y1h,styles[i],c='C1')
    ax.set_xlabel('Flow (ml/min)')
    ax.set_ylabel('Uptake $N$ (ml/min)')
    ax.legend()
    fig.tight_layout()
    fig.savefig(OutFilename)
    return fig, ax


def calc_number_functional_units(InputFluxFiles):
    Nfu = np.zeros(len(InputFluxFiles))
    n = 0
    for f in InputFluxFiles:
        FlowDict = pickle.load(open(f,'rb'))
        mean_flux = np.mean(FlowDict["TV_fluxes"])
        Nfu[n] = np.sum(FlowDict["TV_fluxes"] > 0.001*mean_flux)
        n+=1

    return Nfu
                    
