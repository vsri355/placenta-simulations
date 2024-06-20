using PyCall
using SparseArrays
using SpecialFunctions
using Distributions

py"""
import _pickle as cPickle
import pickle
import bz2

def load_pickle(file):
    return pickle.load(file)

def decompress(fpath):
    return bz2.BZ2File(fpath,'rb')

def pickle_dump(data,fpath):
    with open(fpath,'wb') as f:
        pickle.dump(data,f)

def compress_pickle_dump(data,fpath):
    with bz2.BZ2File(fpath,'w') as f:
        cPickle.dump(data, f)         #write
"""

load_pickle = py"load_pickle"
decompress = py"decompress"
pickle_dump = py"pickle_dump"
compress_pickle_dump = py"compress_pickle_dump"

function calc_visc_factor(radius::Float64, hb::Float64)   #radius in um
    controlhb = 0.45
    beta = 4.0/(1.0 + exp(-0.0593*(2.0*radius - 6.74)))
    visc_factor = (1.0 + (exp(hb*beta) - 1.0)/(exp(controlhb*beta) - 1.0)*
                  (110.0*exp(-2.848*radius) + 3.0 - 3.45*exp(-0.07*radius)))/4.0

    return visc_factor
end

function solve_mesoscale_tree(numgens::Int, numconvolutes::Int, numparallel::Int,
    numparallel_cap::Int, numseries::Int, art_rad_in::Float64, ven_rad_in::Float64,
    press_in::Float64 = 1.0, press_out::Float64 = 0.0; update_resistance::Bool = true)

    int_length = 1.5      #mm %Length of each intermediate villous
    int_radius = 0.030/2.0    #%radius of each intermediate villous
    seg_length = int_length/numconvolutes #; %lengh of each intermediate villous segment

    cap_length = 3.0/numparallel_cap   #%mm %length of individual capillary
    cap_rad = (0.0144)/2.0   # ! %radius of individual capillary

    visc_factor = 1.0
    mu = 0.400e-02     #0.33600e-02 #; %viscosity Pa s

    int_rad_ain = art_rad_in    #!mm Unstrained radius of inlet villous
    int_rad_vin = ven_rad_in   #!mm radius of ouutlet intermediate villous
    int_rad_aout =  0.03/2.0   #! mm radius of mature intermediate villous artery
    int_rad_vout = 0.03

    # total_cap_volume = 0.0
    # total_art_volume = 0.0
    # total_vein_volume = 0.0
    # total_cap_surface_area = 0.0
    # total_art_surface_area = 0.0
    # total_vein_surface_area = 0.0

    #!In the Interface 2015 model the resistance of a "terminal capillary convolute) is defined anatomically by Leiser
    #!Or, we can use the measured resistance of a "terminal capillary convolute" from Erlich
      # if(capillary_model_type.eq.3)then
    cap_resistance = 5.6e5
    R_cap =  cap_resistance*numseries/numparallel
    
      # else if(capillary_model_type.eq.2)then
    # R_cap=(8.0*mu*visc_factor*cap_length)/(np.pi*cap_rad**4.0)*numseries/\
    #             (numparallel_cap*numparallel)    #!%resistance of each capillary convolute segment
    # cap_resistance = ((8.0*mu*visc_factor*cap_length)/(np.pi*cap_rad**4.0))/numparallel_cap
      # endif

    p_unknowns = 2*numconvolutes*numgens
    q_unknowns = numgens*numconvolutes+1 #!Only need to calculate convolute and total flows
    MatrixSize = p_unknowns+q_unknowns
    NonZeros = Int(q_unknowns*4 +  p_unknowns*q_unknowns/2 + p_unknowns*2 - 5)
    irows = zeros(NonZeros)
    icols = zeros(NonZeros)
    entries = zeros(NonZeros)
    #!bttom row + capillaries, +  pyramid entries + bcs +pbalance

    #!Initialise solution matrix arrays
    RHS = zeros(MatrixSize)
    update_resist = zeros(Int,(NonZeros,3))
    original_resist = zeros(NonZeros)
    #!Put matrix together
    nnz = 1  #count entries
    radfac = (int_rad_aout-int_rad_ain)/numgens
    for ng in 1:numgens
        #!Radius of branching villi at theis generation
        int_radius_gen = int_rad_ain + (ng-1)*radfac
        #!Resistance of a segment of branching villi between terminal convolute segments
        visc_factor = calc_visc_factor(1000*int_radius_gen, 0.45)
        R_seg = (8.0*mu*visc_factor*seg_length)/(pi*(int_radius_gen^4.0))     #! %resistance of each intermediate villous segment
        #Pa s mm^-3
        
        #!Add volume and surface area to this segment
        # total_art_volume = total_art_volume + (np.pi*int_radius_gen**2.0*seg_length*numconvolutes)*(2.0**ng)
        # total_art_surface_area = total_art_surface_area + np.pi*2.0*int_radius_gen*seg_length*numconvolutes*(2.0**ng)
        #!Add corresponding "capillary" or terminal conduit volume and surface area to this segment
        # if(capillary_model_type.eq.3)then
        # total_cap_volume = total_cap_volume + 0.65E-3*numconvolutes*numseries*numparallel*(2.0**ng)
        # total_cap_surface_area = total_cap_surface_area + 0.135*numconvolutes*numseries*numparallel*(2.0**ng)
        # else if(capillary_model_type.eq.2)then
        #    total_cap_volume = total_cap_volume +   PI*cap_rad**2.0*dble(cap_length)*dble(numconvolutes)*&
        #           dble(numseries)*dble(numparallel)*dble(numparallel_cap)*(2.0**(ng))
        #   total_cap_surface_area = 2.0*PI*cap_rad*cap_length*dble(numconvolutes)*dble(numseries)*dble(numparallel)* &
        #           dble(numparallel_cap)*(2.0**dble(ng))
        # endif
        #!Add matrix entries for each convolute within this generation
        for nc in 1:numconvolutes
            i = (ng-1)*numconvolutes + nc   #!row number
            #!outward branches (arteries)
            if ((nc == 1) && (ng == 1)) #!Inlet -pressure out -QR =-Pressure in
                RHS[i] = -press_in
                # SparseRow(nnz) = 1
                # SparseCol(nnz) = i
                # SparseVal(nnz) = -1.0
                irows[nnz] = i
                icols[nnz] = i
                entries[nnz] = -1.0
                nnz += 1
                # SparseCol(nnz) = MatrixSize
                # SparseVal(nnz) = -R_seg/2.0 !divide by 2 because bifurcating here
                irows[nnz] = i
                icols[nnz] = MatrixSize
                entries[nnz] = -R_seg/2.0 #!divide by 2 because bifurcating here
                if update_resistance
                    update_resist[nnz,1] = 1    #!indicating this one needs to be updated and is an artery
                    update_resist[nnz,2] = 1    #!power for resistance division
                    update_resist[nnz,3] = i    #!segment that this resistance is associated with
                    original_resist[nnz] = entries[nnz]
                end
                nnz += 1
            else
                # SparseCol(nnz) = i-1
                # SparseVal(nnz) = 1.0
                irows[nnz] = i
                icols[nnz] = i-1
                entries[nnz] = 1.0
                nnz += 1
                # SparseCol(nnz) = i
                # SparseVal(nnz) = -1.0
                irows[nnz] = i
                icols[nnz] = i
                entries[nnz] = -1.0
                nnz += 1
                divider = 2.0^(ng-1)   #!FLOW DIVIDES BY 2
                curgen = ng
                checkgen = 1
                for j in 1:(i-1)
                    # SparseCol(nnz) = p_unknowns + j
                    # SparseVal(nnz) = R_seg/divider
                    irows[nnz] = i
                    icols[nnz] = p_unknowns+j
                    entries[nnz] = R_seg/divider
                    if update_resistance
                        update_resist[nnz,1] = 1    #!indicating this one needs to be updated and is an artery
                        update_resist[nnz,2] = ng   #!current generation
                        update_resist[nnz,3] = j   #!segment that this resistance is associated with
                        original_resist[nnz] = entries[nnz]
                    end
                    nnz += 1
                    if j == checkgen*numconvolutes
                        checkgen = checkgen + 1
                        curgen = curgen - 1
                        divider = 2.0^(curgen - 1)
                    end
                end
                # SparseCol(nnz) = MatrixSize
                # SparseVal(nnz) = -R_seg/(2.0**ng)
                irows[nnz] = i
                icols[nnz] = MatrixSize
                entries[nnz] = -R_seg/(2.0^ng)
                if update_resistance
                    update_resist[nnz,1] = 1    #!indicating this one needs to be updated and is an artery
                    update_resist[nnz,2] = ng   #!current generation
                    update_resist[nnz,3] = i  #!segment that this resistance is associated with
                    original_resist[nnz] = entries[nnz]
                end
                nnz += 1
            end
        end
    end


    #!Repeat for veins
    for ng in 1:numgens
        int_radius_gen = int_rad_vin + (int_rad_vout-int_rad_vin)/numgens*(ng-1.0)
        visc_factor = calc_visc_factor(1000*int_radius_gen, 0.45)
        R_seg=(8.0*mu*visc_factor*seg_length)/(pi*(int_radius_gen^4.0))   #! %resistance of each intermediate villous segment
        # total_vein_volume = total_vein_volume + pi*(int_radius_gen^2.0)*seg_length*numconvolutes*(2.0^ng)
        # total_vein_surface_area = total_vein_surface_area + pi*2.0*int_radius_gen*seg_length*numconvolutes*(2.0^ng)
        for nc in 1:numconvolutes
            i = (ng-1)*numconvolutes + nc      #!pointer to arteries row number
            #!inward branches (veins)
            i2 = Int(p_unknowns/2 + (ng-1)*numconvolutes + nc)
            # SparseCol(nnz) = i2
            # SparseVal(nnz) = 1.0
            irows[nnz] = i2
            icols[nnz] = i2
            entries[nnz] = 1.0
            nnz += 1
            if i2 == p_unknowns
                RHS[i2] = press_out      #assume p is 0 at outlet
                # SparseCol(nnz) = MatrixSize
                # SparseVal(nnz) = -R_seg/2.0
                irows[nnz] = i2
                icols[nnz] = MatrixSize
                entries[nnz] =  -R_seg/2.0
                if update_resistance
                    update_resist[nnz,1] = 2    #!indicating this one needs to be updated and is a vein
                    update_resist[nnz,2] = numgens - ng + 1     #!current generation
                    update_resist[nnz,3] = i2   #!segment that this resistance is associated with
                    original_resist[nnz] = entries[nnz]
                end
                nnz += 1
            else
                #SparseCol(nnz) = i2+1
                #SparseVal(nnz) = -1.0
                irows[nnz] = i2
                icols[nnz] = i2+1
                entries[nnz] = -1.0
                nnz += 1
                checkgen = 1
                curgen = Int(ceil((p_unknowns - i2)/numconvolutes))
                divider = 2.0^(curgen-1)
                for j in 1:(p_unknowns - i2)
                    if j == (p_unknowns - i2)
                        #SparseCol(nnz) = p_unknowns+j
                        #SparseVal(nnz) = R_seg
                        irows[nnz] = i2
                        icols[nnz] = p_unknowns+j
                        entries[nnz] = R_seg
                        if update_resistance
                            update_resist[nnz,1] = 2    #!indicating this one needs to be updated and is a vein
                            update_resist[nnz,2] = numgens - ng + 1     #!current generation
                            update_resist[nnz,3] = i2   #!segment that this resistance is associated with
                            original_resist[nnz] = entries[nnz]
                        end
                        nnz += 1
                    else
                        #SparseCol(nnz) = p_unknowns + j
                        #SparseVal(nnz) = R_seg/divider
                        irows[nnz] = i2
                        icols[nnz] = p_unknowns+j
                        entries[nnz] = R_seg/divider
                        if update_resistance
                            update_resist[nnz,1] = 2    #!indicating this one needs to be updated and is a vein
                            update_resist[nnz,2] = numgens - ng + 1     #!current generation
                            update_resist[nnz,3] = i2  #!segment that this resistance is associated with
                            original_resist[nnz] = entries[nnz]
                        end
                        nnz += 1
                    end

                    if j == checkgen*numconvolutes
                        checkgen = checkgen + 1
                        curgen = curgen - 1
                        divider = 2.0^(curgen-1.0)
                    end
                end
                curgen = Int(ceil((p_unknowns - i2)/numconvolutes))
                # SparseCol(nnz) = MatrixSize
                # SparseVal(nnz) = -R_seg/(2.0**curgen)
                irows[nnz] = i2
                icols[nnz] = MatrixSize
                entries[nnz] = -R_seg/(2.0^curgen)
                if update_resistance
                    update_resist[nnz,1] = 2    #!indicating this one needs to be updated and is a vein
                    update_resist[nnz,2] = numgens - ng + 1     #!current generation
                    update_resist[nnz,3] = i2   #!segment that this resistance is associated with
                    original_resist[nnz] = entries[nnz]
                end
                nnz += 1
            end
        end
    end


    for ng in 1:numgens
        for nc in 1:numconvolutes
            #!Capillary segments
            i = (ng-1)*numconvolutes + nc   #!pointer to arteries row number
            # SparseCol(nnz) = i
            # SparseVal(nnz) = 1.0
            irows[nnz] = p_unknowns+i
            icols[nnz] = i
            entries[nnz] = 1.0
            nnz += 1
            # SparseCol(nnz) = p_unknowns + 1 - i
            # SparseVal(nnz) = -1.0
            irows[nnz] = p_unknowns+i
            icols[nnz] = p_unknowns-i + 1
            entries[nnz] = -1.0
            nnz += 1
            # SparseCol(nnz) = p_unknowns+i
            # SparseVal(nnz) = -R_cap!/2.0
            irows[nnz] = p_unknowns+i
            icols[nnz] = p_unknowns+i
            entries[nnz] = -R_cap  #/2.0
            if update_resistance
                update_resist[nnz,1] = 3    #!indicating this one needs to be updated and is a vein
                update_resist[nnz,2] = ng   #!current generation
                update_resist[nnz,3] = i    #!segment that this resistance is associated with
                original_resist[nnz] = entries[nnz]
            end
            nnz += 1
        end
    end

    for ng in 1:numgens
        for nc in 1:numconvolutes
            #!Conservation of flow
            i = (ng-1)*numconvolutes + nc    #!pointer to arteries row number
            # SparseCol(nnz) = p_unknowns+i
            # SparseVal(nnz) = 2.0**ng
            irows[nnz] = MatrixSize
            icols[nnz] = p_unknowns+i
            entries[nnz] = 2.0^ng
            nnz += 1
        end
    end

    # SparseCol(nnz) = MatrixSize
    # SparseVal(nnz) = -1.0
    irows[nnz] = MatrixSize
    icols[nnz] = MatrixSize
    entries[nnz] = -1.0
    nnz += 1

    #do solve
    sparse_smat = sparse(irows,icols,entries)
    Solution = sparse_smat\RHS

    if update_resistance
        h=0.8
        elastance = 5.0E5
        count_its = 0
        converged = false
        while converged == false
            for nnz in 1:NonZeros
                if update_resist[nnz,1] == 1   #artery
                    if update_resist[nnz,3] == 1
                        p_in = press_in - 11.0*133.0
                    else
                        p_in  = Solution[update_resist[nnz,3]-1] - 11.0*133.0
                    end
                    p_out = Solution[update_resist[nnz,3]] - 11.0*133.0
                    int_radius_gen = int_rad_ain + (update_resist[nnz,2]-1.0)*
                                     (int_rad_aout-int_rad_ain)/numgens

                    r_in = int_radius_gen + 3.0*(int_radius_gen^2)*p_in/(4.0*elastance*h)
                    r_out = int_radius_gen + 3.0*(int_radius_gen^2)*p_out/(4.0*elastance*h)
                    r_ave = (r_in + r_out)/2.0
                    entries[nnz] = original_resist[nnz]*(int_radius_gen^4.0)/(r_ave^4.0)
                else
                    if update_resist[nnz,1] == 2  #vein
                        if update_resist[nnz,3] == p_unknowns
                            p_out = press_out - 11.0*133.0
                        else
                            p_out  = Solution[update_resist[nnz,3]+1] - 11.0*133.0
                        end
                        p_in = Solution[update_resist[nnz,3]] - 11.0*133.0
                        int_radius_gen = int_rad_vin + (update_resist[nnz,1]-1.0)*
                                         (int_rad_vout-int_rad_vin)/numgens

                        r_in = int_radius_gen + 3.0*(int_radius_gen^2)*p_in/(4.0*elastance*h)
                        r_out = int_radius_gen + 3.0*(int_radius_gen^2)*p_out/(4.0*elastance*h)
                        r_ave = (r_in + r_out)/2.0
                        entries[nnz] = original_resist[nnz]*(int_radius_gen^4.0)/(r_ave^4.0)
                    else
                        if (update_resist[nnz,1] == 3)   #capillary
                            # if(capillary_model_type == 2):   #We have a simple branching model and can update capillary resistamc
                            #     p_in = Solution[update_resist[nnz,2]] - 11.0*133.0
                            #     p_out = Solution[p_unknowns-update_resist[nnz,2]+1] - 11.0*133.0
                            #     r_in = cap_rad + 3.0*(cap_rad**2)*p_in/(4.0*elastance*h)
                            #     r_out = cap_rad + 3.0*(cap_rad**2)*p_out/(4.0*elastance*h)
                            #     r_ave = (r_in + r_out)/2.0
                            #     solver_mat[] = original_resist[nnz]*(cap_rad**4.0)/(r_ave**4.0)
                            # else:     #!We have an anatomically derived Erlich Model
                            entries[nnz] = original_resist[nnz]
                        end
                    end
                end
            end

            #Solve system
            sparse_smat = sparse(irows,icols,entries)
            SolutionNew = sparse_smat\RHS
            # sparse_smat = scipy.sparse.coo_matrix((entries, (irows, icols)), shape=(MatrixSize, MatrixSize))
            # SolutionNew = sparse.linalg.spsolve(sparse_smat,RHS)
            err = 0.0
            for i in 1:MatrixSize
                err = err + (SolutionNew[i] - Solution[i])^2.0/Solution[i]^2.0
            end
            err = err/MatrixSize
            if err < 1.0e-6
                converged = true
            else
                if count_its > 20
                    converged = true
                end
            end
            Solution = SolutionNew
            count_its += 1
        end
    end
    #get solution parts
    psol = Solution[1:p_unknowns]
    qsol = Solution[(p_unknowns+1):(p_unknowns+q_unknowns)]       #flow through an individual bud (and final entry = total flow)
    pdrop = psol[1:numgens*numconvolutes] - 
            psol[reverse(collect((numgens*numconvolutes+1):2*numgens*numconvolutes))]

    resistance = (press_in-press_out)/(Solution[MatrixSize])
    #!!Effective length of capillary system
    terminal_length = resistance*(pi*(0.03)^4.0)/(8.0*mu)
    weights = zeros(length(qsol))
    for ng in 1:numgens
        weights[((ng-1)*numconvolutes + 1):(ng*numconvolutes)] .= 2^ng
    end

    #Confirmed that this is all correct

    return pdrop, qsol, resistance, weights
end

function calculate_and_save_allfluxes(OutputPath::String, TerminalProperties::Matrix{PyObject},
    Labels::Array{String,1}, FluxVals::Array{Float64,1},
    Nparallel::Int, Nseries::Int, Nparallel_cap::Int, Nconv::Int, Ngens::Int; 
    use_full_res::Bool=true)

    Ntrees = size(TerminalProperties)[2]
    Nfluxes = length(FluxVals)
    #ignore any(?) flux vals where there is a network missing
    N_meso = Ngens*Nconv+1
    for n in 1:Ntrees, i in 1:Nfluxes
        if TerminalProperties[i,n] != 0
            print("Tree no ", n," flux val: ", FluxVals[i],'\n')
            TermPdrop = TerminalProperties[i,n][1]
            TermFlux = TerminalProperties[i,n][2]
            TermRadIn = TerminalProperties[i,n][3]
            TermRadOut = TerminalProperties[i,n][4]
            TermPin = TerminalProperties[i,n][5]
            TermPout = TerminalProperties[i,n][6]
            
            Nterm = length(TermPdrop)
            Qh = zeros(N_meso*Nterm)
            Ntv = zeros(N_meso*Nterm)

            for kt in 1:Nterm
                P_meso, Q_meso, Res, N_weights = solve_mesoscale_tree(Ngens, Nconv, Nparallel,
                        Nparallel_cap, Nseries, TermRadIn[kt], TermRadOut[kt], TermPin[kt], 
                        TermPout[kt]; update_resistance=use_full_res)   #need to consider occluded case
                if use_full_res
                    Qh[((kt-1)*N_meso+1):(kt*N_meso)] = Q_meso 
                else
                    Qactual = NetworkDict["Flows"][eout+1]
                    Qh[((kt-1)*N_meso+1):(kt*N_meso)] = Qactual*(Q_meso/Q_meso[end])  #scale to actual flow
                end
                Ntv[((kt-1)*N_meso+1):(kt*N_meso)] = N_weights

                kt += 1
            end
            filename = string(OutputPath,"tv_fluxes_",Labels[n],"_flow",string(Int(FluxVals[i])))
            if use_full_res == false
                filename = string(filename,"_approx")
            end
            filename = string(filename,".pkl")
            Data = Dict([("TV_fluxes",Qh),("Nweights",Ntv)])
            pickle_dump(PyDict(Data),filename)
        end
    end
end

function update_min(List::Array{Float64,1}, MinNow::Float64)
    MinLocal = min(List...)
    if MinLocal < MinNow
        MinNow = MinLocal
    end
    return MinNow
end

function update_max(List::Array{Float64,1}, MaxNow::Float64)
    MaxLocal = max(List...)
    if MaxLocal > MaxNow
        MaxNow = MaxLocal
    end
    return MaxNow
end

function calculate_meso_hetero_sensitivity(Filename::String, Fluxes::Array{Float64,1}, 
    NetworksTermData::Matrix{PyObject}, Nconv::Int ,Ngens::Int)
    #meso heterogeneity cases
    #fit the meso tree parameters
    NPvals = [1,2,3,4,5]
    NSvals = [1,2,3,4,5]
    #sample over pmean and pdrop
    maxPav = -1E6
    minPav = 1E6
    maxPdrop = -1E6
    minPdrop = 1E6
    maxRin = 0.0
    maxRout = 0.0
    for (i, q) in enumerate(Fluxes)
        if NetworksTermData[i,3] != 0
            Pin = NetworksTermData[i,3][5]
            Pout = NetworksTermData[i,3][6]
            Pav = 0.5 .* (Pin .+ Pout)
            Pdrop = (Pin .- Pout)
            minPav = update_min(Pav, minPav)
            maxPav = update_max(Pav, maxPav)
            minPdrop = update_min(Pdrop, minPdrop)
            maxPdrop = update_max(Pdrop, maxPdrop)

            Rin = NetworksTermData[i,3][3]
            Rout = NetworksTermData[i,3][4]
            maxRin = update_max(Rin, maxRin)
            maxRout = update_max(Rout, maxRout)
        end
    end
    print(maxRin,' ', maxRout)
    
    Pmeanvals = minPav .+ (maxPav - minPav).*collect(0:0.025:1.001)
    Pdropvals = minPdrop .+ (maxPdrop - minPdrop).*collect(0:0.025:1.001)
    NQmeso = Ngens*Nconv + 1
    Qmeso_vals = zeros((length(NPvals),length(NSvals),length(Pmeanvals),length(Pdropvals),NQmeso))
    Qmeso_weights = zeros((length(NPvals),length(NSvals),length(Pmeanvals),length(Pdropvals),NQmeso))
    for (inp,npar) in enumerate(NPvals)
        for (ins,ns) in enumerate(NSvals)
            print(inp,' ',ins,'\n')
            for (ipm,pmean) in enumerate(Pmeanvals)
                for (ipd,pdrop) in enumerate(Pdropvals)
                    Pin = pmean + 0.5*pdrop
                    Pout = pmean - 0.5*pdrop
                    P_meso, Q_meso, Res, N_weights = solve_mesoscale_tree(Ngens, Nconv, npar,
                        npar, ns, maxRin, maxRout, Pin, Pout; update_resistance=true) 
                    Qmeso_vals[inp,ins,ipm,ipd,:] = Q_meso
                    Qmeso_weights[inp,ins,ipm,ipd,:] = N_weights
                end
            end
        end
    end
    Data = Dict([("Nparallel",NPvals),("Nserial",NSvals),("Pmean",Pmeanvals),
                 ("Pdrop",Pdropvals),("LoopOrder",["Nparallel","Nserial","Pmean","Pdrop"]),
                 ("Qmeso",Qmeso_vals),("Qweights",Qmeso_weights)])

    pickle_dump(PyDict(Data),Filename)
end

function get_vals_from_lognormal_cdf(mu::Float64, sigma::Float64, offset::Float64, N_meso::Int)
    Dx = (1.0/N_meso)
    LowerCDFVals = collect(0:Dx:(1.0 - Dx))
    UpperCDFVals = collect(Dx:Dx:1.0)
    LowerXVals = quantile.(Ref(LogNormal(mu,sigma)),LowerCDFVals)
    UpperXVals = quantile.(Ref(LogNormal(mu,sigma)),UpperCDFVals)
    alpha = 0.5*exp(mu + sigma^2/2)
    beta = (mu + sigma^2)/(sqrt(2)*sigma)
    EV = alpha.*(erf.(beta .- log.(LowerXVals) ./ (sqrt(2)*sigma)) .- erf.(beta .- log.(UpperXVals) ./ (sqrt(2)*sigma))) ./
            (UpperCDFVals .- LowerCDFVals)

    return (offset .+ EV)
end

function run_meso_lognormal_loop(SummaryOutputPath::String, FluxOutputPath::String, 
        Fluxes::Array{Float64,1}, Networks::Matrix{PyObject}, TermProps::Matrix{PyObject}, 
        std_vals::Array{Float64,1}, Ngens::Int, Nconv::Int)
    #run for a range of sigma values
    Ntrees = 1
    Nfluxes = length(Fluxes)
    N_meso = 2*Nconv*(2^Ngens-1)   #total number of terminal villi
    Offset_vals = 1/N_meso .- 0.06261.*std_vals
    Offset_vals[Offset_vals .< 0] .= 0.0
    NetworksOut = Array{PyObject,2}(undef, Nfluxes, length(std_vals))
    TermPropsOut = Array{PyObject,2}(undef, Nfluxes, length(std_vals))
    #ignore any(?) flux vals where there is a network missing
    for (j,std) in enumerate(std_vals)
        mean = 1.0/N_meso - Offset_vals[j]
        mu = log(mean^2/sqrt(std^2 + mean^2))
        sigma = sqrt(log(1 + std^2/mean^2))
        Q_meso_rel = get_vals_from_lognormal_cdf(mu,sigma,Offset_vals[j],N_meso)
        for n in 1:Ntrees, i in 1:Nfluxes
            if TermProps[i,n] != 0
                Tfluxes = TermProps[i,n][2]
                print("Flux check ", Fluxes[i], ' ', sum(Tfluxes)*60/1000, "\n")
                Qh = vec(Tfluxes' .*  Q_meso_rel)
                Ntv = ones(N_meso*length(Tfluxes))
                filename = string(FluxOutputPath,"tv_fluxes_meso_std_",j,
                                  "_flow",string(Int(Fluxes[i])))
                filename = string(filename,".pkl")
                Data = Dict([("TV_fluxes",Qh),("Nweights",Ntv)])
                pickle_dump(PyDict(Data),filename)
            end
            TermPropsOut[i,j] = TermProps[i,n]
            NetworksOut[i,j] = Networks[i,n]
        end
    end
    filename = string(SummaryOutputPath,"meso_hetero_sweep.pbz2")
    data = PyObject([Fluxes,NetworksOut])
    compress_pickle_dump(data,filename)
    data2 = PyObject([TermPropsOut])
    filename2 = string(SummaryOutputPath,"meso_hetero_sweep_terminal.pbz2")
    compress_pickle_dump(data2,filename2)
end



