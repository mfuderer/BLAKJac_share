## Script for generating optimized sequences 
#  Note: "manual_run" variable has to be cycled manually from 1 to 8

include("setup.jl")
recon_options = Dict() # erase all existing settings
nsweeps = 6
nky = 224
nTR = round(Int64, nsweeps*nky)
ky = 1.0 .*(repeat(1:nky, inner=1, outer=nsweeps)); 
kz = ones(nsweeps*nky)
trajectorySet = BLAKJac.TrajectorySet(ky,kz)
BLAKJac.BLAKJac_defaults!(trajectorySet, recon_options)
dataFolder = "/home/mfuderer/Documents/Julia/Capture"

manual_run = 7  # to be increased manually from 1 to 8
# (phase, startstate, cyclic, RMSflip, nametag)
tasklist = [(true, 1, false, 40, "anh"),
            (true, 1, true,  40, "aph"),
            (true,-1, false, 40, "aih"),
            (true, 1, false, 20, "anl"),
            (false,1, false, 40, "nnh"),
            (false,1, true,  40, "nph"),
            (false,-1,false, 40, "nih"),
            (false,1, false, 20, "nnl")]
(slow_phase,start_state,cyclic,RMSflip,nametag) = tasklist[manual_run]

recon_options["opt_slow_phase"] = slow_phase
recon_options["startstate"] = start_state
recon_options["considerCyclic"] = cyclic 
recon_options["sar_limit"]  = RMSflip^2/recon_options["TR"]

recon_options["useSymmetry"] = true     
recon_options["TR"]      = 0.01
recon_options["sigma_ref"] = 1.4 
recon_options["optpars"]   = Optim.Options(time_limit = 3000, iterations = 100000, f_tol=1.0e-10, g_tol = 1.0e-10)  
recon_options["opt_criterion"] = "sar-limited noise"  
recon_options["emphasize_low_freq"] = true 
recon_options["opt_focus"] = "max"      
recon_options["opt_initialize"] = "cRandom30" 
recon_options["opt_complex"] = false   
recon_options["opt_account_maxFlip"] = true 
recon_options["opt_keep_positive"] = false                           
recon_options["TW"] = 0.0             # 0  
recon_options["opt_emergeCriterion"] = 2000
ph = [] 
recon_options["opt_imposed_2nd_derivative_of_phase"] = ph
recon_options["opt_iterations_limit"] = 1
recon_options["sizeSteps"] = [6]   
nRealizations = 6 

fn_base = "20231024"*nametag
for i in 1:nRealizations
    stageText = ""
    portionRange = 0:0
    fn = dataFolder*"/$fn_base($i).jld2"
    RFdeg = BLAKJac.BLAKJac_optimize(trajectorySet, recon_options, i);
    FileIO.save(fn,"RFdeg",RFdeg)
end     

scores=zeros(nRealizations,2)
saved_H = Dict()
recon_options["rfName"]  = "from_file"
recon_options["rfFunction"] = rfDictionary[recon_options["rfName"]]
rfFunction = recon_options["rfFunction"]
(fig,ax)=(subplots(Int(ceil(nRealizations/3)),3,figsize=(9,3)))
for i in 1:nRealizations
    fn = "$fn_base($i)"
    recon_options["rfFile"]  = fn
    RFdeg = rfFunction(recon_options["nTR"], recon_options["nky"])
    ax[i].plot(abs.(RFdeg))
    anglesdd = zeros(length(RFdeg))
    for i in 1:length(RFdeg)-2
        anglesdd[i] = (rad2deg(angle(conj(RFdeg[i])*RFdeg[i+1]*RFdeg[i+1]*conj(RFdeg[i+2])))+270.0) % 180.0 -90.0
    end
    ax[i].plot(anglesdd)
    noises, ItotAll = BLAKJac.BLAKJac_analysis!(cpu, RFdeg, trajectorySet, recon_options, saved_H) 
    @show noises, mean(noises), ItotAll
    scores[i,1] = noises[3] # maximum(noises[2:3])
    scores[i,2] = ItotAll
end
