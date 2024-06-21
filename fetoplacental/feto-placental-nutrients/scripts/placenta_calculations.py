import bz2
import os
import pickle
from itertools import islice

import _pickle as cPickle
import compress_placenta_tree_data as cmpl
import matplotlib.pyplot as plt
import numpy as np
import placenta_tree_extract_process_results_updated as pl
import scipy as sp
import scipy.sparse as sparse
import scipy.stats as stats
from sklearn import linear_model 
import pandas as pd

L_vil_default = 0.016   #length scale m
Lc_vil = 0.002  #in m
alphac = 5.5

#O2 #for now we just use the first value, but might be better do do param uncertainty
B_O2 = [140,140]     #advection constant
Dt_O2 = [2E-9,2E-9]   #diffusivity in m^2/s
Sol_O2 = (1000/760)*0.0238 #solubility converted from bunsen sol to ml/L from Power 1968
DPP_O2 = 77  #partial pressure difference Hill 1973
DC_O2calc = DPP_O2*Sol_O2*1000  #concentration difference in mL / m3
DC_O2 = [DC_O2calc,DC_O2calc]

#CO2
B_CO2 = [1.0,10.0]     #advection constant
Dt_CO2 = [1.9E-9,1.9E-9]  #diffusivity in m^2/s
DPP_CO2 = 40  #partial pressure difference Hill 1973
Sol_CO2 = 0.7  #ml/L
DC_CO2calc = DPP_CO2*Sol_CO2*1000  #concentration difference in mL / m3
DC_CO2 = [DC_CO2calc,DC_CO2calc]

#CO
B_CO = [1E4,1E4]
Dt_CO = Dt_O2  #diffusivity in m^2/s
Sol_CO = (1000/760)*0.0189    #solubility converted from bunsen sol to ml/L from Power 1968
DPP_CO = 0.025 #partial pressure difference from wiki but no source - absolute PP (not drop)
DC_CO = [DPP_CO*Sol_CO*1000,DPP_CO*Sol_CO*1000]   #concentration difference in ml / m3

#Heat
B_heat = [1.0,1.0]
Dt_heat = [1.4E-7,1.4E-7]  #diffusivity in m^2/s
DC_heat = [1.0,1.0]

#Urea
B_urea = [1.0,1.0]
Dt_urea = [1.4E-9,1.4E-9]
DCUcalc = 2.5*0.06006*1000/1.32  #urea 2.5 mmol/L & 60.06 g/mol & 1.32 g/ml -> ml/m3
DC_urea = [DCUcalc,DCUcalc]

#Lactic acid
B_lactic = [1.0,1.0]
Dt_lactic = [1.0E-9, 1.0E-9]
DCLcalc = 1.0*0.009008*1000/1.21  #LA 1.0 mmol/L & 90.08 g/mol & 1.21 g/ml -> ml/m3
DC_lactic = [DCLcalc,DCLcalc]   #ml/m3

#Fructose
B_fructose = [1.0,1.0]
Dt_fructose = [1.0E-13,1.0E-12]
DC_fructose = [1.0,1.0]   #edit

#Glucose
B_glucose = [1.0,1.0]
Dt_glucose= [1.0E-12,1.0E-11]
DC_glucose = [1.0,1.0] #edit

SoluteNames = ['O$_2$','CO','Urea','CO$_2$', 'Fructose', 'Glucose', 'Lactic Acid', 'Heat']
vtkSoluteNames = ['O2','CO','Urea','CO2', 'Fructose', 'Glucose', 'Lactic_Acid', 'Heat']
SoluteColours = ['r','k','y','b','c','m','g','C1']
DC = np.array([np.mean(DC_O2), np.mean(DC_CO),  np.mean(DC_urea), \
            np.mean(DC_CO2), np.mean(DC_fructose), np.mean(DC_glucose),\
            np.mean(DC_lactic), np.mean(DC_heat)])

B = np.array([np.mean(B_O2), np.mean(B_CO),  np.mean(B_urea), \
            np.mean(B_CO2), np.mean(B_fructose), np.mean(B_glucose),\
            np.mean(B_lactic), np.mean(B_heat)])

Dt = np.array([np.mean(Dt_O2), np.mean(Dt_CO),  np.mean(Dt_urea), \
            np.mean(Dt_CO2), np.mean(Dt_fructose), np.mean(Dt_glucose),\
            np.mean(Dt_lactic), np.mean(Dt_heat)])

def mu(Dt, Dp, L_vil = L_vil_default):
    return Dt*L_vil / (Dp*Lc_vil)

def inv_damkohler(Flux, Dt, B, Nseries, Nparallel, L_vil = L_vil_default):
    """inverse Damkohler for each element in series"""
    #FLux in mm3/s
    return Flux*B/(1E9*Dt*Nparallel*L_vil)

def inv_damkohler_alt(Flux, Dt, B, Nseries, Nparallel, L_vil = L_vil_default):
    """inverse Damkohler of all elements in series assuming this scales L_vil"""
    #FLux in mm3/s
    return Flux*B/(1E9*Dt*Nseries*Nparallel*L_vil)

def Max_net_flux(Flux, B, DC):
    return Flux * B * DC / 1E9

def Max_diff_cap(Dt, DC, Nseries, Nparallel, L_vil = L_vil_default):
    return Dt * L_vil * Nseries * Nparallel * DC

def Net_flux_calc(Dam1, Flux, Dt, B, Dp, DC, Nseries, Nparallel, L_vil = L_vil_default):
    Da = 1 / Dam1
    muh = mu(Dt, Dp, L_vil)
    C = (muh / Da)**(2/3) / alphac
    Nmax = Max_net_flux(Flux, B, DC)
    lam = 1 - 1.0/((1.0 - np.exp(-Da))**(-1) + C)
    Nh = Nmax * (1 - lam**Nseries)

    return Nh

def Net_flux(Flux, Dt, B, Dp, DC, Nseries, Nparallel, L_vil = L_vil_default):
    """Net flux where series units act as if flux fully recovers, i.e.
       each results in the same relative change in concentration"""
    Dam1 = inv_damkohler(Flux, Dt, B, Nseries, Nparallel, L_vil)
    if isinstance(Dam1,np.ndarray):   #array
        Nh = np.zeros(len(Dam1))
        Nh[Dam1 > 0] = Net_flux_calc(Dam1[Dam1 > 0], Flux[Dam1 > 0], Dt, B, Dp, DC, \
                                     Nseries, Nparallel, L_vil)
    else:   #scalar
        if Dam1 > 0:
            Nh =  Net_flux_calc(Dam1, Flux, Dt, B, Dp, DC, Nseries, Nparallel, L_vil)
        else:
            Nh = 0


    return Nh

def Net_flux_alt(Flux, Dt, B, Dp, DC, Nseries, Nparallel, L_vil = L_vil_default):
    """Net flux where series units act as if all connected in one unit, i.e.
       total vessel length is multiplied by Nseries"""
    Da = 1 / inv_damkohler_alt(Flux, Dt, B, Nseries, Nparallel, L_vil)
    muh = mu(Dt, Dp, L_vil)
    C = (muh / Da)**(2/3) / alphac
    Nmax = Max_net_flux(Flux, B, DC)
    lam = 1 - 1.0/((1.0 - np.exp(-Da))**(-1) + C)
    Nh = Nmax * (1 - lam)
    return Nh

def all_pressure_drops(TermPdrop, TermFlows, Pdrop_meso, Q_meso):
    """calculate relative meso-tree pressures for each terminal unit"""
    N = len(Pdrop_meso)
    Pall = np.zeros(N*len(TermPdrop))
    Qall = np.zeros(N*len(TermPdrop))
    Qmtot = np.sum(Q_meso)      #sum of flows with 1.0 Pa pressure drop

    for i in np.arange(len(TermPdrop)):
        Pall[i*N:(i+1)*N] = Pdrop_meso * TermPdrop[i] #scale to total pressure drop
        Qall[i*N:(i+1)*N] = (Q_meso / Qmtot) * TermFlows[i] #scale to actual flow

    return Pall, Qall

def convolute_with_lognormal(TermPdrop, TermFlows, N, ustd, offset, CapRes):
    """"""
    umean = 1/N - offset  #this gets multiplied by flux in
    mu = np.log(umean**2/np.sqrt(ustd**2 + umean**2))
    sigma = np.sqrt(np.log(1 + ustd**2/umean**2))

    x = (0.5 + np.arange(N))/N
    v = stats.lognorm.ppf(x,sigma,loc=offset,scale=np.exp(mu))
    #v = (1/N)*v/np.mean(v)

    vtotal = np.zeros(N*len(TermFlows))

    for j in np.arange(len(TermFlows)):
        vtotal[j*N:(j+1)*N] = v*TermFlows[j]

    return CapRes*vtotal, vtotal

def transport_up_tree(cterm, network):
    """calculate q*adjacency to calculate steady state conc everywhere on tree"""
    Nnodes = network['Nodes'].shape[1]
    row = []
    col = []
    entries = []
    for n in np.arange(Nnodes):     #flux in - flux out = 0
        if n in network['Term_nodes']:
            row.append(n)
            col.append(n)
            entries.append(1.0)
        else:
            flux_in = 0
            flux_out = 0
            for j in network['Edges_in'][n]:
                if network['Flows'][j] > 0:
                    row.append(n)
                    col.append(network['Edges'][0,j])
                    entries.append(network['Flows'][j])
                    flux_in += network['Flows'][j]
                else:
                    flux_out += -network['Flows'][j]
            for j in network['Edges_out'][n]:
                if network['Flows'][j] < 0:
                    row.append(n)
                    col.append(network['Edges'][1,j])
                    entries.append(-network['Flows'][j])
                    flux_in += -network['Flows'][j]
                else:
                    flux_out += network['Flows'][j]
            if flux_out == 0:
                flux_out = flux_in
            row.append(n)
            col.append(n)
            entries.append(-flux_out)

    qadjacency = sparse.csc_matrix((entries, (row, col)), shape=(Nnodes,Nnodes))
    b = np.zeros(Nnodes)
    b[network['Term_nodes']] = cterm
    cnodes = sparse.linalg.spsolve(qadjacency,b)
    return cnodes

def weighted_median(f, weights):
    iorder = np.argsort(f)
    wsum = np.concatenate((np.array([0]),np.cumsum(weights[iorder])))
    fsort = f[iorder]
    if wsum[-1]%2 == 1:
        imed = (wsum[-1] + 1)/2 #median is a single value
        b1 = (wsum[0:-1] < imed)
        b2 = (wsum[1:] >= imed)
        return fsort[b1 & b2][0]   
    else:
        ilow = wsum[-1]/2
        iupp = ilow + 1
        b1low = (wsum[0:-1] < ilow)
        b2low = (wsum[1:] >= ilow)
        b1upp = (wsum[0:-1] < iupp)
        b2upp = (wsum[1:] >= iupp)
        return 0.5*(fsort[b1low & b2low][0] + fsort[b1upp & b2upp][0])

def create_vessel_superedges(network):
    current_nodes = network["Entry_nodes"]  #get nodes to begin with
    node_assigned = np.zeros(np.shape(network["Nodes"])[1],dtype=bool)   #keep track of which nodes have been asssigned to superedges
    SEs = []   #list of superedges
    current_SEs = np.arange(len(current_nodes)) 
    for (k,se) in enumerate(current_SEs):
        SEs.append(np.array([current_nodes[k]],dtype=int))
        node_assigned[current_nodes[k]] = True

    while len(current_nodes) > 0:
        next_nodes = np.zeros(0,dtype=int)
        next_SEs = np.zeros(0,dtype=int)
        for (ik, k) in enumerate(current_nodes):
            eout = network["Edges_out"][k]  #could be more than one
            if len(eout) > 0:
                kout = network["Edges"][1,eout]
                if len(kout) == 1:   #continue current superedge
                    if node_assigned[kout[0]] == False:   #node not already assigned to a superedge
                        SEs[current_SEs[ik]] = np.append(SEs[current_SEs[ik]],kout[0])  #add to current SE
                        next_nodes = np.append(next_nodes,kout[0])   #add to next loop
                        next_SEs = np.append(next_SEs,current_SEs[ik])
                        node_assigned[kout[0]] = True
                else:
                    for ikout in np.arange(len(kout)):   #
                        if node_assigned[kout[ikout]] == False:
                            next_SEs = np.append(next_SEs,len(SEs))   #add index of new SE
                            SEs.append(np.array([k,kout[ikout]],dtype=int))   #add start of new SE
                            next_nodes = np.append(next_nodes,kout[ikout])  #add to next loop
                            node_assigned[kout[ikout]] = True
        current_nodes = next_nodes
        current_SEs = next_SEs
    
    network["VesselSuperedges"] = SEs

def plot_tree_vtk(filehead, network, *, NodeData = np.array([]), EdgeData = np.array([]),\
                  TermData = np.array([]), NodeNames = np.array([]), EdgeNames = np.array([]), \
                  TermNames = np.array([])):
    #output network (dict format) into .vtk format
    #produces two files "filehead_tree.vtk" (a network of edges and nodes)
    #               and "filehead_term.vtk" (terminal nodes containing terminal properties)

    create_vessel_superedges(network)
    filename = ("%s_tree.vtk"%filehead)
    f  = open(filename,'w')
    #write vtk legacy header
    f.write("# vtk DataFile Version 2.0\nTree data %s\nASCII\nDATASET POLYDATA\n\n"%filehead)
    Npts = len(network['Nodes'][0])
    f.write("POINTS %d float\n"%Npts)
    for k in np.arange(Npts):
        f.write("%f %f %f\n"%(network['Nodes'][0][k],network['Nodes'][1][k],network['Nodes'][2][k]))

    NSedges = len(network['VesselSuperedges'])
    lensum = 0
    for j in np.arange(NSedges):
        lensum += 1 + len(network['VesselSuperedges'][j])
    f.write("\nLINES %d %d\n" % (NSedges,lensum))
    for j in np.arange(NSedges):
        f.write("%d"%len(network['VesselSuperedges'][j]))
        for k in network['VesselSuperedges'][j]:
            f.write(" %d"%k)
        f.write('\n')

    f.write("\nPOINT_DATA %d\n" % (Npts))
    f.write("\nSCALARS Pressure float 1\nLOOKUP_TABLE default\n")
    for k in np.arange(Npts):
        f.write("%.3e\n"%network['Pressure'][k])

    #extra node data (if given) is printed here
    if len(NodeData) > 0:
        if len(np.shape(NodeData)) > 1:
            ncols = np.shape(NodeData)[1]
        else:
            ncols = 1
            NodeData = np.reshape(NodeData,(np.shape(NodeData)[0],1))
        for nc in np.arange(ncols):
            f.write("\nSCALARS ")
            if len(NodeNames) > nc:
                f.write("%s "%NodeNames[nc])
            else:
                f.write("NodeData%d "%nc)
            f.write("float 1\nLOOKUP_TABLE default\n")
            for k in np.arange(Npts):
                f.write("%.3e\n"%NodeData[k,nc])

    #need to convert all the edge data to node data
    f.write("\nSCALARS Flux float 1\nLOOKUP_TABLE default\n")
    for k in np.arange(Npts):
        tot_flow = 0
        count = 0
        for jin in network["Edges_in"][k]:
            tot_flow += network['Flows'][jin]
            count += 1
        for jout in network["Edges_out"][k]:
            tot_flow += network['Flows'][jout]
            count += 1 
        f.write("%.3e\n"%(tot_flow/count))

    f.write("\nSCALARS Radius float 1\nLOOKUP_TABLE default\n")
    for k in np.arange(Npts):
        tot_rad = 0
        count = 0
        for jin in network["Edges_in"][k]:
            tot_rad += network['Radius'][jin]
            count += 1
        for jout in network["Edges_out"][k]:
            tot_rad += network['Radius'][jout]
            count += 1 
        f.write("%.3e\n"%(tot_rad/count))

    if len(EdgeData) > 0:
        if len(np.shape(EdgeData)) > 1:
            ncols = np.shape(EdgeData)[1]
        else:
            ncols = 1
            EdgeData = np.reshape(EdgeData,(np.shape(EdgeData)[0],1))
        for nc in np.arange(ncols):
            f.write("\nSCALARS ")
            if len(EdgeNames) > nc:
                f.write("%s "%EdgeNames[nc])
            else:
                f.write("EdgeData%d "%nc)
            f.write("float 1\nLOOKUP_TABLE default\n")
        for k in np.arange(Npts):
            tot = 0
            count = 0
            for jin in network["Edges_in"][k]:
                tot += EdgeData[jin,nc]
                count += 1
            for jout in network["Edges_out"][k]:
                tot += EdgeData[jout,nc]
                count += 1 
            f.write("%.3e\n"%(tot/count))

    # f.write("\nCELL_DATA %d\n" % (Nedges))

    # f.write("\nSCALARS Flux float 1\nLOOKUP_TABLE default\n")
    # for j in np.arange(Nedges):
    #     f.write("%.3e\n"%network['Flows'][j])

    # f.write("\nSCALARS Radius float 1\nLOOKUP_TABLE default\n")
    # for j in np.arange(Nedges):
    #     f.write("%.3e\n"%network['Radius'][j])

    # #extra edge data (if given) is printed here
    # if len(EdgeData) > 0:
    #     if len(np.shape(EdgeData)) > 1:
    #         ncols = np.shape(EdgeData)[1]
    #     else:
    #         ncols = 1
    #         EdgeData = np.reshape(EdgeData,(np.shape(EdgeData)[0],1))
    #     for nc in np.arange(ncols):
    #         f.write("\nSCALARS ")
    #         if len(EdgeNames) > nc:
    #             f.write("%s "%EdgeNames[nc])
    #         else:
    #             f.write("EdgeData%d "%nc)
    #         f.write("float 1\nLOOKUP_TABLE default\n")
    #         for k in np.arange(Nedges):
    #             f.write("%.3e\n"%EdgeData[k,nc])

    # f.write("\nVECTORS EdgeVec float\nLOOKUP_TABLE default\n")
    # for j in np.arange(Nedges):
    #     kout = network['Edges'][1][j] 
    #     kin = network['Edges'][0][j]
    #     direction = np.zeros(3)
    #     for n in np.arange(3):
    #         direction[n] = network['Nodes'][n][kout] - network['Nodes'][n][kin]
    #     f.write("%f %f %f\n"%(direction[0],direction[1],direction[2]))

    # f.close()

    #terminal nodes file
    filenamet = ("%s_term.vtk"%filehead)
    ft = open(filenamet,'w')
    ft.write("# vtk DataFile Version 2.0\nTerm data %s\nASCII\nDATASET UNSTRUCTURED_GRID\n\n"%filehead)
    Ntnodes = len(network['Term_nodes'])
    ft.write("POINTS %d float\n"%Ntnodes)
    for k in network['Term_nodes']:
        ft.write("%f %f %f\n"%(network['Nodes'][0][k],network['Nodes'][1][k],network['Nodes'][2][k]))


    ft.write("\nPOINT_DATA %d\n" % (Ntnodes))
    ft.write("\nSCALARS Pressure float\nLOOKUP_TABLE default\n")
    for k in network['Term_nodes']:
        ft.write("%.3e\n"%network['Pressure'][k])

    #extra data (if given) for terminal nodes
    if len(TermData) > 0:
        if len(np.shape(TermData)) > 1:
            ncols = np.shape(TermData)[1]
        else:
            ncols = 1
            TermData = np.reshape(TermData,(np.shape(TermData)[0],1))
        for nc in np.arange(ncols):
            ft.write("\nSCALARS ")
            if len(TermNames) > nc:
                ft.write("%s "%TermNames[nc])
            else:
                ft.write("TermData%d "%nc)
            ft.write("float\nLOOKUP_TABLE default\n")
            for k in np.arange(Ntnodes):
                ft.write("%.3e\n"%TermData[k,nc])

    ft.close()


def calc_all_uptakes(Network, TermData, Flux, PathToFlux, Solutes, TreeFilename, \
    TreeLabel, Ngens, Nconv, Nparallel ,Nparallel_cap, Nseries, Dt, DC, B,\
    vtk_print = False, L_vil_val = L_vil_default):
    NSolutes = len(Solutes)
    Nmeso = Ngens*Nconv
    TermFlux = TermData[1]
    AllTVFluxes = pickle.load(open(PathToFlux,'rb'))
    NonZero = (AllTVFluxes['Nweights'] > 0)
    TVFluxes = AllTVFluxes['TV_fluxes'][NonZero]
    TVWeights = AllTVFluxes['Nweights'][NonZero]
    InvDaGlobal = np.zeros(NSolutes)
    Ntot = np.zeros(NSolutes)
    Nmaxtot = np.zeros(NSolutes)
    Nideal = np.zeros(NSolutes)
    if vtk_print:   #write VTK
        Tdata = np.zeros((len(TermFlux),1+3*len(Solutes)))
        Ndata = np.zeros((Network['Nodes'].shape[1],len(Solutes)))
        Names = np.zeros(1+3*len(Solutes),dtype=object)
        NNames = np.zeros(len(Solutes),dtype=object)
        Names[0] = 'Flux'
        Tdata[:,0] = TermFlux
    for s in np.arange(NSolutes):
        NTermTotal = np.sum(TVWeights)
        Qtot = np.sum(TVFluxes*TVWeights)
        InvDaGlobal[s] = inv_damkohler(Qtot, Dt[s], B[s], \
                            Nseries, Nparallel, L_vil_val)/(NTermTotal)
        Qideal = Qtot/(NTermTotal)
        Nideal[s] = NTermTotal*Net_flux(Qideal, Dt[s], B[s], Dt[s], DC[s], Nseries, Nparallel, L_vil_val)
        InvDa = inv_damkohler(TVFluxes, Dt[s], B[s], Nseries, Nparallel, L_vil_val)
        Nf = Net_flux(TVFluxes, Dt[s], B[s], Dt[s], DC[s], Nseries, Nparallel, L_vil_val)
        Conc = Nf*1E9/(TVFluxes*B[s])
        Ntot[s] = np.sum(Nf*TVWeights)
        print('Flux check:', Qtot*60/1000, ' = ', Flux)

        Nmax = Max_diff_cap(Dt[s], DC[s], Nseries, Nparallel, L_vil_val)
        Nmaxtot[s] = Nmax*NTermTotal
        if vtk_print:
            vtksolname = ''
            for smatch in np.arange(len(vtkSoluteNames)):
                if SoluteNames[smatch] == Solutes[s]:
                    vtksolname = vtkSoluteNames[smatch]
                    break
            Names[1+s] = ('Inv_Da_' + vtksolname + '_median')
            Names[1+len(Solutes)+s] = ('Uptake_' + vtksolname)
            Names[1+2*len(Solutes)+s] = ('Conc_' + vtksolname + '_median')
            for k in np.arange(len(TermFlux)):
                Tdata[k,1+s] = np.median(InvDa[(k*Nmeso):((k+1)*Nmeso)])
                Tdata[k,1+len(Solutes)+s] = 60*np.sum(Nf[(k*Nmeso):((k+1)*Nmeso)])   #ml/min
                Tdata[k,1+2*len(Solutes)+s] = np.median(Conc[(k*Nmeso):((k+1)*Nmeso)])
            NNames[s] = ('Conc' + vtksolname)
            ch = np.transpose(Tdata[:,1+2*len(Solutes)+s])
            Ndata[:,s] = transport_up_tree(ch, Network)
    if vtk_print:   #write VTK
        filehead = TreeFilename + "_np%d"%(Nparallel) + "ns%d"%(Nseries) +  "_flux_%d"%(Flux)
        plot_tree_vtk(filehead, Network, NodeData = Ndata,\
                NodeNames = NNames, TermData = Tdata, TermNames = Names)
        
    return Ntot, Nmaxtot, InvDaGlobal, Nideal

def calculate_and_save_measures(OutputName, TreePath, TermPath, FluxDir, TreeFileNames, TreeLabels, \
                                Solutes, DC, B, Dt, Nparallel, Nseries, Nparallel_cap, \
                                Nconv, Ngens, vtk_print = np.array([])):

    [Fluxes, Networks] = pickle.load(bz2.BZ2File(TreePath,'rb'))
    [TermData] = pickle.load(bz2.BZ2File(TermPath,'rb'))
    Networks = np.array(Networks)
    Fluxes = np.array(Fluxes)
    TermData = np.array(TermData)
    [NFluxes,NTrees] = np.shape(Networks)
    NSolutes = len(Solutes)
    Ntot = np.zeros((NFluxes,NTrees,NSolutes),dtype=float)
    Nmaxtot = np.zeros((NFluxes,NTrees,NSolutes),dtype=float)
    Nideal = np.zeros((NFluxes,NTrees,NSolutes),dtype=float)
    InvDaGlobal = np.zeros((NFluxes,NTrees,NSolutes),dtype=float)
    print(NFluxes,NTrees)
    for j in np.arange(NTrees):
        for i in np.arange(NFluxes):
            if Networks[i,j] != 0:
                bool_vtk = (Fluxes[i] in vtk_print)
                print(TreeLabels[j], Fluxes[i], bool_vtk)
                TVPath = os.path.join(FluxDir,'tv_fluxes_' + TreeFileNames[j] + '_flow' + str(int(Fluxes[i])) + '.pkl')
                print(TVPath)
                Ntoth, Nmaxtoth, InvDaGlobalh, Nidealh = calc_all_uptakes(Networks[i,j], TermData[i,j],\
                        Fluxes[i], TVPath, Solutes, FluxDir + TreeFileNames[j], TreeLabels[j], Ngens, Nconv,\
                        Nparallel, Nparallel_cap, Nseries, Dt, DC, B, bool_vtk)
                Ntot[i,j,:] = Ntoth
                Nmaxtot[i,j,:] = Nmaxtoth
                InvDaGlobal[i,j,:] = InvDaGlobalh
                Nideal[i,j,:] = Nidealh
            else:
                for s in np.arange(len(Solutes)):
                    InvDaGlobal[i,j,:] = np.nan
                    Ntot[i,j,:] = np.nan
                    Nmaxtot[i,j,:] = np.nan
                    Nideal[i,j,:] = np.nan

    Output =  {'Fluxes': Fluxes, 'Solutes': Solutes, 'Total Exchange': \
            Ntot, 'Max Exchange': Nmaxtot, 'Ideal Exchange': Nideal, \
            'Inverse Damkohler Global': InvDaGlobal, 'Tree Labels': TreeLabels}
            #other stuff

    with open(OutputName,'wb') as f:
        pickle.dump(Output, f)

def calculate_L_sensitivity(Nparallel,Nseries,Nparallel_cap,Nconv,Ngens):
    #loop over L, keep flow fixed
    L_scales = np.arange(0.2,5,0.1)
    L_vals = L_scales*L_vil_default
    Ns = len(SoluteNames)
    NL = len(L_vals)
    TreePath = '../Compressed_data/tree_flux_sweep_series.pbz2'
    FluxDir = '../Compressed_data/Series_trees'
    [Fluxes, Networks, TermData] = pickle.load(bz2.BZ2File(TreePath,'rb'))
    [NFluxes,NTrees] = np.shape(Networks)
    Ntot = np.zeros((NL,NTrees,Ns))
    Nmaxtot = np.zeros((NL,NTrees,Ns))
    Nideal= np.zeros((NL,NTrees,Ns))
    InvDaGlobal= np.zeros((NL,NTrees,Ns))
    iFluxes = np.arange(len(Fluxes))
    i = iFluxes[np.abs(Fluxes - 250) < 1E-3][0]
    for j in np.arange(NTrees):
        TVPath = os.path.join(FluxDir,'tv_fluxes_' + cmpl.SerialTreeFilenames[j] +\
                 '_flow' + str(int(Fluxes[i])) + '.pkl')
        FlowDict = pickle.load(open(TVPath,'rb'))
        NTermTotal = np.sum(FlowDict['Nweights'])
        for il in np.arange(NL):
            L_vil = L_vals[il]
            Ntoth, Nmaxtoth, InvDaGlobalh, Nidealh = calc_all_uptakes(Networks[i,j], TermData[i,j],\
                    Fluxes[i], TVPath, SoluteNames, cmpl.SerialTreeFilenames[j], cmpl.SerialTreeLabels[j],\
                    Ngens, Nconv, Nparallel, Nparallel_cap, Nseries, Dt, DC, B, L_vil_val = L_vil)
            Ntot[il,j,] = Ntoth
            Nmaxtot[il,j,:] = Nmaxtoth
            InvDaGlobal[il,j,:] = InvDaGlobalh
            Nideal[il,j,:] = Nidealh


    Output =  {'LValues': L_vals, 'Solutes': SoluteNames, 'Total Exchange': \
                Ntot, 'Max Exchange': Nmaxtot, 'Ideal Exchange': Nideal, \
                'Inverse Damkohler Global': InvDaGlobal, 'Tree Labels': cmpl.SerialTreeLabels}

    #print to file
    with open('../Compressed_data/Lsweep_outcomes.pkl','wb') as f:
        pickle.dump(Output, f)


def do_minimisation(obs,guess):

    def MLERegression(params):
        offset, std = params[0], params[1] # inputs are guesses at our parameters
        mean = 1/len(obs) - offset
        mu = np.log(mean**2/np.sqrt(std**2 + mean**2))
        sigma = np.sqrt(np.log(1 + std**2/mean**2))
        negLL = -np.sum( stats.lognorm.logpdf(obs,sigma,loc=offset,scale=np.exp(mu)))# return negative LL

        return(negLL)
    
    return sp.optimize.minimize(MLERegression, guess, method ='Nelder-Mead', options={'disp': True})

def do_minimisation(obs,guess,xlims):

    def MLERegression(params):
        if (params[0] < xlims[0,0]) or (params[0] > xlims[0,1]) or (params[1] < xlims[1,0]) or (params[1] > xlims[1,1]):
            return 1E6
        offset, std = params[0], params[1] # inputs are guesses at our parameters
        mean = 1.0/len(obs) - offset
        mu = np.log(mean**2/np.sqrt(std**2 + mean**2))
        sigma = np.sqrt(np.log(1 + std**2/mean**2))
        negLL = -np.sum( stats.lognorm.logpdf(obs,sigma,loc=offset,scale=np.exp(mu)))# return negative LL

        return(negLL)
    
    return sp.optimize.minimize(MLERegression, guess, method ='Nelder-Mead')

def calculate_meso_hetero_sensitivity(input_path):
    MesoDict = pickle.load(open(input_path,'rb'))
    Keys = MesoDict["LoopOrder"]
    Nx = np.array([len(MesoDict[Keys[0]]),len(MesoDict[Keys[1]]),\
                len(MesoDict[Keys[2]]),len(MesoDict[Keys[3]])])
    X = np.zeros((Nx[0]*Nx[1]*Nx[2]*Nx[3],4))
    y1 = np.zeros(Nx[0]*Nx[1]*Nx[2]*Nx[3])
    y2 = np.zeros(Nx[0]*Nx[1]*Nx[2]*Nx[3])
    ymu = np.zeros(Nx[0]*Nx[1]*Nx[2]*Nx[3])
    ysigma = np.zeros(Nx[0]*Nx[1]*Nx[2]*Nx[3])
    for (i1,v1) in enumerate(MesoDict[Keys[0]]):
        for (i2,v2) in enumerate(MesoDict[Keys[1]]):
            for (i3,v3) in enumerate(MesoDict[Keys[2]]):
                for (i4,v4) in enumerate(MesoDict[Keys[3]]):
                    Qvalsh = MesoDict["Qmeso"][i1,i2,i3,i4][:-1]
                    Qweightsh = MesoDict["Qweights"][i1,i2,i3,i4][:-1]
                    Qweightsh = (np.rint(Qweightsh)).astype(int)
                    Qtot = np.sum(Qvalsh*Qweightsh)
                    Qvalsnorm = Qvalsh/Qtot
                    
                    obs = np.repeat(Qvalsnorm,Qweightsh)
                    # guesses and optimize
                    guess = np.array([0.99*np.min(obs), np.std(obs)])
                    xlims = np.zeros((2,2))
                    xlims[0,:] = np.array([1E-6,1.0/len(obs)])
                    xlims[1,:] = np.array([1E-6,10.0])
                    results = do_minimisation(obs,guess,xlims)
                    n = i1*Nx[1]*Nx[2]*Nx[3] + i2*Nx[2]*Nx[3] + i3*Nx[3] + i4
                    X[n,:] = np.transpose([v1,v2,v3,v4])
                    y1[n] = results.x[0]
                    y2[n] = results.x[1]
                    mean = 1.0/len(obs) - y1[n]
                    ymu[n] = np.log(mean**2/np.sqrt(y2[n]**2 + mean**2))
                    ysigma[n] = np.sqrt(np.log(1 + y2[n]**2/mean**2))
                    # print(results.x)

    Z = {'offset':y1,'std':y2}
    for i in np.arange(4):
        Z[Keys[i]] = X[:,i]
    df = pd.DataFrame(Z)
    df.to_csv('lognormal_fitting_outcomes.csv')
     

if __name__ == "__main__":
    pass
    #put parameter upper and lower values?
    


    # calculate_meso_hetero_sensitivity()
    # meso hetero trees