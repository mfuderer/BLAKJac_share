####
# File created 2023-08-08 to regenerate the set of graphs intended for resubmission of phase-benefits paper


## from 2023-01-17 Graphs of 6 different realizations of the optimization 
# taken as example: Amplitude+Phase no-pauze 
# "Fig 9" 
include("setup.jl")
PyPlot.rc("font", family="serif")
recon_options = Dict() # erase all existing settings
nsweeps = 6                                            
nky = 224
nTR = round(Int64, nsweeps*nky)
ky = 1.0 .*(repeat(1:nky, inner=1, outer=nsweeps)); 
kz = ones(nsweeps*nky)
trajectorySet = BLAKJac.TrajectorySet(ky,kz)
BLAKJac.BLAKJac_defaults!(trajectorySet, recon_options)

recon_options["useSymmetry"] = true     
recon_options["TR"]      = 0.01
recon_options["sigma_ref"] = 1.4 # See logbook 20220815
recon_options["emphasize_low_freq"] = false            # this is for optimization; for evaluation, use it without
recon_options["startstate"]     = 1
recon_options["considerCyclic"] = true
fn_base =  "20220920A" # "20231013Y" #  "20220920A" # change 2023-10-13 !!!!!!!!!!!!!!!!!!!!!!
nRealizations = 6 

scores=zeros(nRealizations,2)
saved_H = Dict()
recon_options["rfName"]  = "from_file"
recon_options["rfFunction"] = rfDictionary[recon_options["rfName"]]
rfFunction = recon_options["rfFunction"]
(fig,ax)=(subplots(Int(ceil(nRealizations/3)),3,figsize=(14,10)))
for i in 1:nRealizations
    fn = "$fn_base($i)"
    recon_options["rfFile"]  = fn
    RFdeg = rfFunction(recon_options["nTR"], recon_options["nky"])
    ax[i].plot(abs.(RFdeg), label="Amplitude")
    anglesdd = zeros(length(RFdeg))
    for i in 1:length(RFdeg)-2
        anglesdd[i] = (rad2deg(angle(conj(RFdeg[i])*RFdeg[i+1]*RFdeg[i+1]*conj(RFdeg[i+2])))+270.0) % 180.0 -90.0
    end
    ax[i].plot(anglesdd, label="ϕ\'\'")
    noises, ItotAll = BLAKJac.BLAKJac_analysis!(cpu, RFdeg, trajectorySet, recon_options, saved_H) 
    @show noises, mean(noises), ItotAll
    scores[i,1] = noises[3] # maximum(noises[2:3])
    scores[i,2] = ItotAll
    @show mean(anglesdd.^2)
end
ax[2].legend(fontsize=18)
for rrr in 1:2
    ax[rrr].set_ylabel("Absolute value (degrees) \n or phase'' (degrees per TR²)", fontsize=18)
end
for ccc in 2:2:6; 
    ax[ccc].set_xlabel("RF pulse number", fontsize=18); 
end;





## From 2022-10-11 ; fig.1
#  ("New set of BLAKJac noise level outputs for sequences used in phantoms 2022-09-23 and (envisaged) in volunteer 2022-10-14,
#  copied from 2022-08-24 script
#  re-used 2022-10-14 for display purposes
#  adapted 2022-10-28
#  messed with 2023-02-08 ")
include("setup.jl")
PyPlot.rc("font", family="serif")

recon_options = Dict() # erase all existing settings
nsweeps = 6                                            # !!!!!!!
nky = 224
nTR = round(Int64, nsweeps*nky)
ky = 1.0 .*(repeat(1:nky, inner=1, outer=nsweeps)); 
kz = ones(nsweeps*nky)
trajectorySet = BLAKJac.TrajectorySet(ky,kz)
BLAKJac.BLAKJac_defaults!(trajectorySet, recon_options)
recon_options["emphasize_low_freq"] = false # Statement (dummy), to emphasize difference to optimization results

recon_options["useSymmetry"] = true     
recon_options["TR"]      = 0.01
recon_options["handleB1"] =  "no" 
recon_options["rfName"]  = "from_file"
cases=["20220818W(5)", "20220815M(5)",  "20220815P(1)", "20220815Q(2)", "20220921B(3)", "20220920A(1)", "20220815S(3)", "20220816T(4)", ]
names=["nih6_818W5"  , "aih6_815M5",       "nnh6"        , "anh6"     , "nph6",        "aph6",        "nnl6_815S3",   "anl6_816T4",   ]
inversions=[ -1,            -1,             1,               1,               1,            1,            1,             1]
cyclic =[false,          false,              false,          false,       true,         true,          false,        false,        ] 
labels =["A-Def",      "B-Def",         "A1",               "B1",         " A2",         "B2",          "A3",         "B3"]

saved_H = Dict()
(fig,ax)=(subplots(Int(ceil(length(cases)/4)),4,figsize=(14,10)))

for ccc in 1:length(cases)
    recon_options["considerCyclic"] = cyclic[ccc]
    recon_options["startstate"] = inversions[ccc]
    recon_options["rfFile"]  = cases[ccc]; 
    name                     = names[ccc]
    recon_options["rfFunction"] = rfDictionary[recon_options["rfName"]]
    rfFunction = recon_options["rfFunction"]
    RFdeg = rfFunction(recon_options["nTR"], recon_options["nky"])
    ax[ccc].plot(abs.(RFdeg),label="Amplitude") #ax[ccc].plot([0.0]) #ax[ccc].plot(abs.(RFdeg))    
    anglesdd = zeros(length(RFdeg))
    for i in 1:length(RFdeg)-2
        anglesdd[i] = (rad2deg(angle(conj(RFdeg[i])*RFdeg[i+1]*RFdeg[i+1]*conj(RFdeg[i+2])))+270.0) % 180.0 -90.0
    end
    ax[ccc].plot(anglesdd, label="ϕ\'\'")
    ax[ccc].text(0.5,0.7,labels[ccc],fontsize=14, transform=ax[ccc].transAxes)
    noises, ItotAll, b1f = BLAKJac.BLAKJac_analysis!(cpu, RFdeg, trajectorySet, recon_options, saved_H) 
    @show noises, mean(noises), name, ItotAll, b1f
    @show abs(RFdeg[1])
end
ax[2].legend(fontsize=14)
ax[1].set_ylabel("Amplitude-only\n\n Absolute value (degrees) \n or phase'' (degrees per TR²)", fontsize=18)
ax[2].set_ylabel("Amplitude+Phase\n\n Absolute value (degrees) \n or phase'' (degrees per TR²)", fontsize=18)
for ccc in 2:2:8; 
    ax[ccc].set_xlabel("RF pulse number", fontsize=18); 
end;
ax[1].text(0,1.17,"Optimized for ...", fontsize=18,      transform=ax[1].transAxes);
ax[1].text(0,1.1,"Initial inversion (Def)", fontsize=18, transform=ax[1].transAxes); # Note: label change 2023-02-08
ax[3].text(0,1.1,"No inversion (1)", fontsize=18, transform=ax[3].transAxes);
ax[5].text(0,1.1,"No pause (2)", fontsize=18, transform=ax[5].transAxes);
ax[7].text(0,1.1,"Low SAR (3)", fontsize=18, transform=ax[7].transAxes);





## Fig 6  
# From 2023-02-08 (...) last modified 2023-04-09
using Statistics
PyPlot.rc("font", family="serif")

mapnames = ["std(T1)", "std(T2)", "rho"];
# names=["aih6","anh6","aph6","anl6","nih6","nnh6","nph6","nnl6"]
casenames =["Initial inversion (1)", "No inversion (2)", "No pause (3)", "Very low SAR (4)"]
pos  =[  1.0,   2.0,    3.0,   4.0,   1.3,  2.3,  3.3,  4.3]
kleur=["orange","orange","orange","orange", "blue","blue","blue","blue"]
refSequence = 5
collection = Dict()
ddd = Dict()

ddd["label"]      = "Phantom \n"
ddd["s_label"]    = "Phantom"
ddd["scatter_color"]="black"
dummy             = ["aih6","anh6","aph6","anl6","nih6","nnh6","nph6","nnl6"]
ddd["deviations"] = [[ 14.1,  17.1,   18.9,  39.3,  18.0, 62.1, 49.0,  129],
                       [2.3,   2.4,    2.7,   8.1,   3.6,  6.0,  4.3,  18.2]] # phantom
ddd["predictions"]= [[  2.4,   2.5,   3.0,    5.0,   3.3,  4.9,  5.2,   16.1],
                     [  2.4,   2.4,   2.7,    5.0,   3.3,  4.3,  4.5,   16.5]]   # 7-points BLAKJac
ddd["normalize"]  = false
ddd["values"]     = []
ddd["casenames"]  = casenames
ddd["description"]= "standard deviation \n in phantom vials (a.u.)"
collection["phantom_N"] = copy(ddd)

# -------------------------------------------------------------
ddd["label"]      = "White \n matter"
ddd["s_label"]      = "White matter"
ddd["scatter_color"]="green"
dummy             = ["aih6","anh6","aph6","anl6","nih6","nnh6","nph6","nnl6"]
ddd["deviations"] = [[   98,   127,  171,    195,    103,   130,  201,   409],
                     [  5.1,   5.2,  7.0,   18.9,   5.7,    6.1, 9.3,   31.1]]   # spatially measured ROI in WM, added 2022-11-03
ddd["predictions"]= [[3.2,  3.2,   4.0,    8.3,   3.2,  4.6,  5.6,   14.2],
                     [4.7,  4.7,   5.2,    9.6,   5.7,  6.2,  7.2,   15.6]]   # added 2022-11-29, 90% (936,58) and 10% (10,10) 
                                                                            # per publication on SPIJN, 2nd component: myelin water
ddd["normalize"]  = false # true
ddd["values"]     = [[1000, 1060,  1110,  790,    940,    1040,  1090,  970],
                     [   40,   44,    44,   50,     35,     40,     44,   48]]
ddd["casenames"]  = casenames
ddd["description"]= "relative standard deviation \n in white matter ROI (a.u.)"
collection["white_matter"] = copy(ddd)

# -------------------------------------------------------------
ddd["label"]      = "Marrow \n (tibia)"
ddd["s_label"]      = "Bone marrow"
ddd["scatter_color"]="orange"
dummy             = ["aih6","anh6","aph6","anl6","nih6","nnh6","nph6","nnl6"]
# ddd["deviations"] = [[   60,    47,   54,     63,    70,    89,   93,    143],
#                      [ 19.1,    35,   23,     26,    31,    33,   36,     39]]   # spatially measured ROI in marrow, added 2023-01-02
# ddd["deviations"] = [[   54,    42,   49,     57,    61,    81,   79,    136],
#                      [ 19.0,    36, 19.4,     25,    29,    31,   32,     42]]   # spatially measured in small ROI in marrow, see 2023-01-03
# ddd["deviations"] = [[   65,    47,   53,     62,    65,    86,   89,    148],
#                      [   20,    26,   29,     24,    30,    35,   35,     46]]   # spatially measured in more superior small ROI in marrow, see 2023-01-11
ddd["deviations"] = [[   68,    49,   54,     63,    70,    87,   96,    151],
                     [   23,    33,   33,     25,    32,    38,   40,     46]]   # on modified ROI position dd. 2023-01-19 per suggestion of team 
ddd["predictions"]= [[  2.3,   2.3,  2.6,    3.2,   3.7,   5.3,  4.9,   10.1],
                     [  1.8,   1.7,  2.0,    2.3,   3.3,   4.7,  4.5,    8.2]]   # predictions on (0.33,0.08)
ddd["normalize"]  = false # true
ddd["values"]     = [[  310,  290,   284,  317,    368,    404,   377,  557],
                     [  128,  147,   131,  104,    145,    151,   144,  154]]
ddd["casenames"]  = ["Initial inversion (1)", "No inversion (2)", "No pause (3)", "Low SAR (4)"]
ddd["description"]= "relative standard deviation \n in tibial bone-marrow ROI (a.u.)"
collection["marrow"] = copy(ddd)

# -------------------------------------------------------------
ddd["label"]      = "Muscle \n (calf)"
ddd["s_label"]      = "Muscle tissue"
ddd["scatter_color"]="red"
dummy             = ["aih6","anh6","aph6","anl6","nih6","nnh6","nph6","nnl6"]
# ddd["deviations"] = [[   59,    73,   95,     83,    59,   104,  120,    118],
#                      [  1.8,   1.7,  2.3,    2.2,   1.7,   1.8,  2.1,    2.2]]   # spatially measured ROI in marrow, added 2023-01-02
ddd["deviations"] = [[   55,    68,   89,     79,    56,    94,  116,    108],
                     [  1.6,   1.5,  2.1,    2.0,   1.5,   1.7,  2.0,    2.0]]   # on modified ROI position dd. 2023-01-19 per suggestion of team
ddd["predictions"]= [[  7.4,   6.4,  8.9,    6.4,   6.7,   8.6, 11.2,    9.1],
                     [  8.0,   6.9,  8.4,    6.6,   8.2,   7.6, 10.6,    8.9]]
ddd["normalize"]  = false # true
ddd["values"]     = [[ 1417, 1282,  1281, 1274,   1304,   1392,  1408, 1371],
                     [   27,   23,    23,   23,     22,     23,    23,   22]]
ddd["casenames"]  = ["Initial inversion (1)", "No inversion (2)", "No pause (3)", "Low SAR (4)"]
ddd["description"]= "relative standard deviation \n in calf muscle ROI (a.u.)"
collection["muscle"] = copy(ddd)
#
for (key,ddd) in collection
    dev    = ddd["deviations"]
    means  = ddd["values"]
    relDev = [(dev[i]./dev[i][refSequence] ) for i in 1:2]
    if ddd["normalize"]
        relDev = [(dev[i]./means[i]./dev[i][refSequence].*means[i][refSequence] ) for i in 1:2]
    end
    bljDev = ddd["predictions"]
    bljRel = [(bljDev[i]./bljDev[i][refSequence])   for i in 1:2]

    brief_names = ["(Def)","(1)","(2)","(3)"]            # added 2023-01-19
    (fig,ax)=(subplots(1,2,figsize=(12,5)))
    fig.subplots_adjust(bottom=0.1) # was 0.3
    fig.subplots_adjust(top=0.95) 
    for m in 1:2
        ax[m].set_xticks(1:length(ddd["casenames"]))
        #ax[m].set_xticklabels(ddd["casenames"], rotation=60, fontsize=14)
        ax[m].set_xticklabels(brief_names, fontsize=16) 
        p1 = ax[m].bar(pos, relDev[m], width=0.3 ,color=kleur) 
        ax[m].scatter(pos, bljRel[m], color="black") 
        ax[m].bar(4.6,7.0,width=0.05,color="white") 
        if ddd["s_label"]=="Phantom"
            ax[m].text(2.3,6.3,"T$m",fontsize=18)
        end
    end
    ax[1].set_ylabel("relative std.", fontsize=18, horizontalalignment="left", y=0.0)
    ax[2].text(-0.3,1.4,ddd["s_label"],  fontsize=18, backgroundcolor="grey")
end


## Fig 7
# (!) dependent on fig. 6 !!!
# from 2023-01-11 "Further analysis on measured noise data in phantom, knee and brain" (messed with dd. 2023-02-03)
using Statistics
comparisonSet = [[(1,2),(5,6)],[(1,3),(5,7)],[(1,5),(2,6),(3,7)]]
task          = [ "compare",     "compare",    "average"]
description   = [ "benefit of inversion pulse",
                                "benefit of pause",
                                            "average benefit of phase"]
barColor     = ["orange","blue"]
textColor    = ["black", "white"]
textValue    = ["Amplitude+Phase", "Amplitude only"]

barHeight = 0.25

figure()
for (setNumber,setElement) in enumerate(comparisonSet)
    sumRatios = zeros(2,length(setElement))
    for (i,refTestPair) in enumerate(setElement)
        (refSequence,testSequence) = refTestPair
        for (key,ddd) in collection
            dev    = ddd["deviations"]
            bljDev = ddd["predictions"]
            for m in 1:2
                sumRatios[1,i] += log(   dev[m][testSequence]/   dev[m][refSequence])
                sumRatios[2,i] += log(bljDev[m][testSequence]/bljDev[m][refSequence])
            end
        end
    end
    sumRatios ./= 2*length(collection)
    @show sumRatios

    sumRatios .*= 100.0

    vPos = 0.0
    if task[setNumber]=="compare"
        for i in 1:2
            verPos = setNumber+barHeight*(i-1.5)
            barh(   verPos, sumRatios[1,i],color=barColor[i], height=barHeight)
            # scatter(sumRatios[2,i], verPos, color="black")                    # removed 2023-02-03 per suggestion Aless
            horPos = (i==2) ? 2.0 : 2.0+sumRatios[1,1]
            text(horPos, verPos-0.03, textValue[i], color=textColor[i])
        end
        #tPos = maximum(sumRatios)
        vPos = setNumber + barHeight + 0.05
    elseif task[setNumber]=="average"
        sumRatios = mean(sumRatios,dims=2)
        barh(        setNumber,            sumRatios[1,1],color="gray", height=barHeight)
        # scatter(     sumRatios[2,1],         setNumber,  color="black")          # removed 2023-02-03 per suggestion Aless
        # tPos = maximum(sumRatios,dims=1)[1,1]
        vPos = setNumber + 0.5*barHeight + 0.05
    end
    text(2.0,vPos,description[setNumber])
end
plot([0,0],[1-barHeight-0.05,length(comparisonSet)+barHeight+0.05], color="black")
yticks([])
#yticklabels(description)
xlabel("SNR benefit (in %)")



## Fig 8
# (!) dependent on fig. 6 !!!
# from "added 2023-01-20: measurement-to-BLAKJac correlation"
markers = ["+","o"]
figure()
allLnM = []
allLnB = []
for (key,ddd) in collection
    dev    = ddd["deviations"]
    bljDev = ddd["predictions"]
    for m in 1:2
        lnm = log.(dev[m])
        lnb = log.(bljDev[m])
        lnm .-= mean(lnm)
        lnb .-= mean(lnb)
        append!(allLnM, lnm)
        append!(allLnB, lnb)
        scatter(lnb,lnm,marker=markers[m],color=ddd["scatter_color"],label=ddd["s_label"]*" T$m")
        slope = cov(lnm,lnb) / cov(lnb);
        spread = std(lnm .- slope.*lnb)
        @show m,ddd["s_label"],cor(lnm,lnb),slope, spread
    end 
end 
plot([-0.7,1.4],[-0.7,1.4],"--",color="grey")
xlabel("Log of BLAKJac-predicted noise level")
ylabel("Log of measured noise level")
legend()

@show cor(allLnM,allLnB)
slope = cov(allLnM,allLnB) / cov(allLnB); @show slope
interc = mean(allLnM) - slope*mean(allLnB); @show interc; 


############################## Prep for MRM Note

## 2023-10-13 Experiment with nonzero waiting time
# Taken from 2022-09-20 "Aanmaken van pauzenloze sequenties" 
include("setup.jl")
recon_options = Dict() # erase all existing settings
nsweeps = 6
nky = 224
nTR = round(Int64, nsweeps*nky)
ky = 1.0 .*(repeat(1:nky, inner=1, outer=nsweeps)); 
kz = ones(nsweeps*nky)
trajectorySet = BLAKJac.TrajectorySet(ky,kz)
BLAKJac.BLAKJac_defaults!(trajectorySet, recon_options)

recon_options["useSymmetry"] = true     
recon_options["TR"]      = 0.01
recon_options["startstate"] = 1  
recon_options["sigma_ref"] = 1.4 #  See logbook 20220815
recon_options["optpars"]   = Optim.Options(time_limit = 3000.0, iterations = 100000, f_tol=1.0e-10, g_tol = 1.0e-10)  
#recon_options["optpars"]   = Optim.Options(time_limit = 3000.0, iterations = 100000, f_tol=1.0e-8, g_tol = 1.0e-8)  
recon_options["opt_criterion"] = "sar-limited noise"  
recon_options["sar_limit"]  = 40^2/0.01   # 40 # 100 T
recon_options["emphasize_low_freq"] = true 
recon_options["opt_focus"] = "max"      
recon_options["opt_initialize"] = "cRandom30" 
recon_options["opt_complex"] = false      
recon_options["opt_account_maxFlip"] = true 
recon_options["opt_keep_positive"] = false                           
recon_options["opt_slow_phase"] =  true 
            # true Y,A,B,D,F,J,L,N,Q,S,U
            # false X,Z,C,E,G,H,I,K,M,P,R                   
recon_options["considerCyclic"] = true   
recon_options["TW"] = 0.0               # 0  # 0,1 for Z, A # 0.25 C,D # 0.5 E,F 
recon_options["opt_emergeCriterion"] = 2000
ph = [] 
#ph = zeros(nTR); ph .= 2.0;
recon_options["opt_imposed_2nd_derivative_of_phase"] = ph
recon_options["opt_iterations_limit"] = 1 # 3 B, H # 2 G,U
recon_options["sizeSteps"] = [6, 10, 15, 20] # [6] 
                # [6, 10, 15]  B,H # [6 10] G,U # 5 I,J # 4 K,L # 3 M,N # 2 P,Q # 1 R,S    
nRealizations = 6 

fn_base = "20231014U"
for i in 1:nRealizations
    stageText = ""
    portionRange = 0:0
    fn = "/home/mfuderer/Documents/Julia/Capture/$fn_base($i).jld2"
    RFdeg = BLAKJac.BLAKJac_optimize(trajectorySet, recon_options, i);
    FileIO.save(fn,"RFdeg",RFdeg)
end 

#fn_base = "20220817U"; nRealizations=6
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
figure(); scatter(scores[:,1], scores[:,2]);

## 2023-10-15 (Continuation of the above - but I ran out of alphabet)
include("setup.jl")
recon_options = Dict() # erase all existing settings
nsweeps = 6
nky = 224
nTR = round(Int64, nsweeps*nky)
ky = 1.0 .*(repeat(1:nky, inner=1, outer=nsweeps)); 
kz = ones(nsweeps*nky)
trajectorySet = BLAKJac.TrajectorySet(ky,kz)
BLAKJac.BLAKJac_defaults!(trajectorySet, recon_options)

recon_options["useSymmetry"] = true     
recon_options["TR"]      = 0.01
recon_options["startstate"] = 1  
recon_options["sigma_ref"] = 1.4 #  See logbook 20220815
time_limit                 = 3000 # 3000 # 5000 A, C ,W
recon_options["optpars"]   = Optim.Options(time_limit = time_limit, iterations = 100000, f_tol=1.0e-10, g_tol = 1.0e-10)  
recon_options["opt_criterion"] = "sar-limited noise"  
recon_options["sar_limit"]  = 40^2/0.01   # 40 
                # 100 T,E # 70 D,F # 30 G,L # 25 H.oopsG,M?,N # 20 I,P # 15 J,Q # 10 K,R
recon_options["emphasize_low_freq"] = true 
recon_options["opt_focus"] = "max"      
recon_options["opt_initialize"] = "cRandom30" 
recon_options["opt_complex"] = false   
recon_options["opt_account_maxFlip"] = true 
recon_options["opt_keep_positive"] = false                           
recon_options["opt_slow_phase"] =  true
            # true A,C,D,L,M,N,P,Q,R,X
            # false B,E,F,G,H,I,J,K,S,T,U,V,W
recon_options["considerCyclic"] = true   
recon_options["TW"] = 0.1             # 0  # 0.1 S # 0.1 nospoil T # 0.1 spoil U,X # 0.25 V
recon_options["opt_emergeCriterion"] = 2000
ph = [] 
#ph = zeros(nTR); ph .= 2.0;
recon_options["opt_imposed_2nd_derivative_of_phase"] = ph
recon_options["opt_iterations_limit"] = 1 #1 #  4 A,B # 2 C # 5 W
recon_options["sizeSteps"] = [6] # [6] # [6, 10, 15,20] A,B # [6,10] C # [6,10,15,20,30] W
                #    
nRealizations = 6 

fn_base = "20231016X"
for i in 1:nRealizations
    stageText = ""
    portionRange = 0:0
    fn = "/home/mfuderer/Documents/Julia/Capture/$fn_base($i).jld2"
    RFdeg = BLAKJac.BLAKJac_optimize(trajectorySet, recon_options, i);
    FileIO.save(fn,"RFdeg",RFdeg)
end     

#fn_base = "20220817U"; nRealizations=6
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

## 2023-10-17 Further continuation, ran out of alphabet again
include("setup.jl")
recon_options = Dict() # erase all existing settings
nsweeps = 6
nky = 224
nTR = round(Int64, nsweeps*nky)
ky = 1.0 .*(repeat(1:nky, inner=1, outer=nsweeps)); 
kz = ones(nsweeps*nky)
trajectorySet = BLAKJac.TrajectorySet(ky,kz)
BLAKJac.BLAKJac_defaults!(trajectorySet, recon_options)

recon_options["useSymmetry"] = true     
recon_options["TR"]      = 0.01
recon_options["startstate"] = 1  
recon_options["sigma_ref"] = 1.4 #  See logbook 20220815
time_limit                 = 3000 # 3000 # 5000 N  
recon_options["optpars"]   = Optim.Options(time_limit = time_limit, iterations = 100000, f_tol=1.0e-10, g_tol = 1.0e-10)  
recon_options["opt_criterion"] = "sar-limited noise"  
recon_options["sar_limit"]  = 20^2/0.01   # 40 # 25 P # 20 Q 
recon_options["emphasize_low_freq"] = true 
recon_options["opt_focus"] = "max"      
recon_options["opt_initialize"] = "cRandom30" 
recon_options["opt_complex"] = false   
recon_options["opt_account_maxFlip"] = true 
recon_options["opt_keep_positive"] = false                           
recon_options["opt_slow_phase"] =  false
            # true A,B,C,D,E,K,M
            # false F,G,H,I,J,L,N,P,Q 
recon_options["considerCyclic"] = true   
recon_options["TW"] = 0.0             # 0  
            # 0.25 A # 0.5 B,F # 1.0 C,G # 2.0 D,H # 3.0 E,I # 5.0 J,K # 10.0 L,M
recon_options["opt_emergeCriterion"] = 2000
ph = [] 
#ph = zeros(nTR); ph .= 2.0;
recon_options["opt_imposed_2nd_derivative_of_phase"] = ph
recon_options["opt_iterations_limit"] = 1 #1 # 6 N 
recon_options["sizeSteps"] = [6]  # [6] # [6,10,15,20,30,40] N 
                #    
nRealizations = 6 

fn_base = "20231017Q"
for i in 1:nRealizations
    stageText = ""
    portionRange = 0:0
    fn = "/home/mfuderer/Documents/Julia/Capture/$fn_base($i).jld2"
    RFdeg = BLAKJac.BLAKJac_optimize(trajectorySet, recon_options, i);
    FileIO.save(fn,"RFdeg",RFdeg)
end     

#fn_base = "20220817U"; nRealizations=6
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


## 2023-10-17 Re-entering the data of above simulations
using PyPlot
PyPlot.rc("font", family="serif")

resols =    [2,    3,    4,     5,    6,   10,   15,   20,   30,   40]
noises_r=[[3.58, 3.43, 3.46, 3.37, 3.23, 3.14, 3.10, 3.02],
          [8.62, 7.85, 6.76, 5.85, 5.66, 5.42, 4.51, 4.47, 4.08, 3.96]]

rms_flip =  [10,  15,   20,   25,   30,   40,  70,  100]
noises_p=[[5.62, 4.13, 3.80, 3.44, 3.38, 3.23, 3.26, 3.27],
          [15.5, 11.0, 7.86, 7.86, 6.46, 5.66, 5.21, 4.90]]

wait_time = [0,  0.1,   0.25, 0.5,  1.0,  2.0,  3.0,  5.0,  10]
noises_w=[[3.23, 3.21, 3.29, 3.21, 3.07, 2.99, 2.92, 2.77, 2.68],
          [5.66, 5.89, 6.03, 6.01, 6.06, 6.06, 5.57, 5.38, 5.24]]
segment_time = 224*6*0.01
noises_wn = [noises_w[i].*sqrt.((wait_time.+segment_time)./segment_time) for i in 1:2]

fig,(ax1,ax2) = subplots(1,2,sharey=true)
#fig.subplots_adjust(vspace=2.0)
ax1.scatter(resols[1:length(noises_r[1])],noises_r[1],color="orange",label="phase")
ax1.scatter(resols[1:length(noises_r[2])],noises_r[2],color="blue",label="no phase")
ax2.scatter([1344],[3.5],color="blue")
ax1.set_ylim(0.0,9.0)
ax2.set_ylim(0.0,9.0)
ax1.set_xlim(0,42)
ax2.set_xlim(1340,1350)
ax2.set_xticks([1344])
ax2.set_xticklabels(["1344"])
ax1.spines["right"].set_visible(false)
ax2.spines["left"].set_visible(false)
ax1.set_ylabel("Noise level [a.u.]", fontsize=18)
ax1.set_xlabel("Cubic-spline resolution", fontsize=18)
ax1.legend()
# ax1.plot([2,1344],[3.5,3.5],"--",color="blue")
# ax2.plot([2,1344],[3.5,3.5],"--",color="blue")
# ax1.plot([2,1344],[3.02,3.02],"--",color="orange")
# ax2.plot([2,1344],[3.02,3.02],"--",color="orange")

figure()
scatter(rms_flip,noises_p[1],color="orange",label="phase")
scatter(rms_flip,noises_p[2],color="blue",label="no phase")
scatter([0.1],[0.1],color="white")
ylabel("Noise level [a.u.]", fontsize=18)
xlabel("RMS excitation angle [deg]", fontsize=18)
legend()

figure()
scatter(wait_time,noises_w[1],color="orange",label="phase")
scatter(wait_time,noises_w[2],color="blue",label="no phase")
scatter(wait_time,noises_wn[1],color="orange",marker="x",label="phase, normalized")
scatter(wait_time,noises_wn[2],color="blue",marker="x",label="no phase, normalized")
scatter([0.1],[0.1],color="white")
ylabel("Noise level [a.u.]", fontsize=18)
xlabel("Wait time between segments [s]", fontsize=18)
legend()


## 2023-10-18 Repeat logging of last optimization run (because logging failed)
include("setup.jl")
recon_options = Dict() # erase all existing settings
nsweeps = 6
nky = 224
nTR = round(Int64, nsweeps*nky)
ky = 1.0 .*(repeat(1:nky, inner=1, outer=nsweeps)); 
kz = ones(nsweeps*nky)
trajectorySet = BLAKJac.TrajectorySet(ky,kz)
BLAKJac.BLAKJac_defaults!(trajectorySet, recon_options)

fn_base = "20231013Y"; nRealizations=6
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

## 2023-10-18 Attempts to make a no-phase full-resolution single-realsation run
include("setup.jl")
recon_options = Dict() # erase all existing settings
nsweeps = 6
nky = 224
nTR = round(Int64, nsweeps*nky)
ky = 1.0 .*(repeat(1:nky, inner=1, outer=nsweeps)); 
kz = ones(nsweeps*nky)
trajectorySet = BLAKJac.TrajectorySet(ky,kz)
BLAKJac.BLAKJac_defaults!(trajectorySet, recon_options)

recon_options["useSymmetry"] = true     
recon_options["TR"]      = 0.01
recon_options["startstate"] = 1  
recon_options["sigma_ref"] = 1.4 #  See logbook 20220815
time_limit                 = 2500 # 3000 # 5000 N # 2500 R 
recon_options["optpars"]   = Optim.Options(time_limit = time_limit, iterations = 100000, f_tol=1.0e-10, g_tol = 1.0e-10)  
recon_options["opt_criterion"] = "sar-limited noise"  
recon_options["sar_limit"]  = 40^2/0.01 
recon_options["emphasize_low_freq"] = true 
recon_options["opt_focus"] = "max"      
recon_options["opt_initialize"] = "cRandom30" 
recon_options["opt_complex"] = false   
recon_options["opt_account_maxFlip"] = true 
recon_options["opt_keep_positive"] = false                           
recon_options["opt_slow_phase"] =  false
recon_options["considerCyclic"] = true   
recon_options["TW"] = 0.0             # 0  
recon_options["opt_emergeCriterion"] = 2000
ph = [] 
#ph = zeros(nTR); ph .= 2.0;
recon_options["opt_imposed_2nd_derivative_of_phase"] = ph
recon_options["opt_iterations_limit"] = 50 #1 # 50 R 
recon_options["sizeSteps"] = [6,10,15,20,30]   # [6] # [6,10,15,20,30] R  
recon_options["portion_size"]=38
nRealizations = 1 # 6 # 1 R 

fn_base = "20231019dummy"
for i in 1:nRealizations
    stageText = ""
    portionRange = 0:0
    fn = "/home/mfuderer/Documents/Julia/Capture/$fn_base($i).jld2"
    RFdeg = BLAKJac.BLAKJac_optimize(trajectorySet, recon_options, i);
    FileIO.save(fn,"RFdeg",RFdeg)
    end     

#fn_base = "20220817U"; nRealizations=6
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


## 2023-10-19 Note: for a remake of the original figures, consider first 335 lines of this 
##  file as well as ScanAnalysisScript4.jl on line 109
## For generation of sequences: this is a bit tricky due to re-use (read: overwriting) of
##   some of the scripts.
##   Useful entry: script 2022-08-18 from BLAKJacScript3,
##   which is actually very similar to the scripts above. Needs double-check, e.g. generate
##   series S (no phase, inversion, 40^2) and compare to 20220818W 

## 2023-10-20 Note: script for former fig. 9 (first in this file) needs no refactoring
## 2023-10-21  ... nor is refactoring needed for script leading to fig. 8

## 2023-10-20 Refactored (not yet modified) script for fig. 1
#   (Adapted to new condition sequence and names 2023-10-27)
include("setup.jl")
PyPlot.rc("font", family="serif")

recon_options = Dict() # erase all existing settings
nsweeps = 6  
nky = 224
nTR = round(Int64, nsweeps*nky)
ky = 1.0 .*(repeat(1:nky, inner=1, outer=nsweeps)); 
kz = ones(nsweeps*nky)
trajectorySet = BLAKJac.TrajectorySet(ky,kz)
BLAKJac.BLAKJac_defaults!(trajectorySet, recon_options)
recon_options["emphasize_low_freq"] = false # Statement (dummy), to emphasize difference to optimization results

recon_options["useSymmetry"] = true     
recon_options["TR"]      = 0.01
recon_options["handleB1"] =  "no" 
recon_options["rfName"]  = "from_file"

# caseSet = [ ("20220818W(5)","nih6_818W5",-1,false,"A-def","Initial inversion (Def)"),
#             ("20220815M(5)","aih6_815M5",-1,false,"B-def","Initial inversion (Def)"),
#             ("20220815P(1)","nnh6" ,      1,false,"A1",   "No inversion (1)"),
#             ("20220815Q(2)","anh6",       1,false,"B1",   "No inversion (1)"),
#             ("20220921B(3)","nph6",       1,true, "A2",   "No pause (2)"),
#             ("20220920A(1)","aph6",       1,true ,"B2",   "No pause (2)"),
#             ("20220815S(3)","nnl6_815S3", 1,false,"A3",   "Low SAR (3)"),
#             ("20220816T(4)","anl6_816T4", 1,false,"B3",   "Low SAR (3)") ]
caseSet = [ ("20220815P(1)","nnh6" ,      1,false,"A-Base",   "No inversion (1)"),
            ("20220815Q(2)","anh6",       1,false,"B-Base",   "No inversion (1)"),
            ("20220921B(3)","nph6",       1,true, "A-NoPause",   "No pause (2)"),
            ("20220920A(1)","aph6",       1,true ,"B-NoPause",   "No pause (2)"),
            ("20220818W(5)","nih6_818W5",-1,false,"A-Invert","Initial inversion (Def)"),
            ("20220815M(5)","aih6_815M5",-1,false,"B-Invert","Initial inversion (Def)"),
            ("20220815S(3)","nnl6_815S3", 1,false,"A-LowSAR",   "Low SAR (3)"),
            ("20220816T(4)","anl6_816T4", 1,false,"B-LowSAR",   "Low SAR (3)") ]

saved_H = Dict()
(fig,ax)=(subplots(Int(ceil(length(caseSet)/4)),4,figsize=(14,10)))

for ccc in eachindex(caseSet)
    (case,name,inversion,cyclic,label,axTxt) = caseSet[ccc]
    recon_options["considerCyclic"] = cyclic
    recon_options["startstate"] = inversion
    recon_options["rfFile"]  = case
    name                     = name
    recon_options["rfFunction"] = rfDictionary[recon_options["rfName"]]
    rfFunction = recon_options["rfFunction"]
    RFdeg = rfFunction(recon_options["nTR"], recon_options["nky"])
    ax[ccc].plot(abs.(RFdeg),label="Amplitude") #ax[ccc].plot([0.0]) #ax[ccc].plot(abs.(RFdeg))    
    anglesdd = zeros(length(RFdeg))
    for i in 1:length(RFdeg)-2
        anglesdd[i] = (rad2deg(angle(conj(RFdeg[i])*RFdeg[i+1]*RFdeg[i+1]*conj(RFdeg[i+2])))+270.0) % 180.0 -90.0
    end
    ax[ccc].plot(anglesdd, label="ϕ\'\'")
    ax[ccc].text(0.5,0.7,label,fontsize=14, transform=ax[ccc].transAxes)
    if ccc%2==1
        ax[ccc].text(0,1.1,axTxt, fontsize=18, transform=ax[ccc].transAxes);
    else
        ax[ccc].set_xlabel("RF pulse number", fontsize=18); 
    end
    noises, ItotAll, b1f = BLAKJac.BLAKJac_analysis!(cpu, RFdeg, trajectorySet, recon_options, saved_H) 
    @show noises, mean(noises), name, ItotAll, b1f
    @show abs(RFdeg[1])
end
ax[2].legend(fontsize=14)
ax[1].set_ylabel("Amplitude-only\n\n Absolute value (degrees) \n or phase'' (degrees per TR²)", fontsize=18)
ax[2].set_ylabel("Amplitude+Phase\n\n Absolute value (degrees) \n or phase'' (degrees per TR²)", fontsize=18)
ax[1].text(0,1.17,"Optimized for ...", fontsize=18,      transform=ax[1].transAxes);


## 2023-10-20 Fig. 6 refactoring 
# From 2023-02-08 (...) last modified 2023-04-09
using Statistics
PyPlot.rc("font", family="serif")

mapnames = ["std(T1)", "std(T2)", "rho"];
pos  =[  1.0,   2.0,    3.0,   4.0,   1.3,  2.3,  3.3,  4.3]
kleur=["orange","orange","orange","orange", "blue","blue","blue","blue"]
brief_names = ["(Def)","(1)","(2)","(3)"]  
resequence  = ["anh6","aph6","aih6","anl6","nnh6","nph6","nih6","nnl6"]
refSequence = 5
collection = Dict()
ddd = Dict()

ddd["label"]      = "Phantom \n"
ddd["s_label"]    = "Phantom"
ddd["scatter_color"]="black"
ddd["key"]        = ["aih6","anh6","aph6","anl6","nih6","nnh6","nph6","nnl6"]
ddd["deviations"] = [[ 14.1,  17.1,   18.9,  39.3,  18.0, 62.1, 49.0,  129],
                       [2.3,   2.4,    2.7,   8.1,   3.6,  6.0,  4.3,  18.2]] # phantom
ddd["predictions"]= [[  2.4,   2.5,   3.0,    5.0,   3.3,  4.9,  5.2,   16.1],
                     [  2.4,   2.4,   2.7,    5.0,   3.3,  4.3,  4.5,   16.5]]   # 7-points BLAKJac
ddd["normalize"]  = false
ddd["values"]     = []
ddd["description"]= "standard deviation \n in phantom vials (a.u.)"
collection["phantom_N"] = copy(ddd)

# -------------------------------------------------------------
ddd["label"]      = "White \n matter"
ddd["s_label"]      = "White matter"
ddd["scatter_color"]="green"
ddd["key"]          = ["aih6","anh6","aph6","anl6","nih6","nnh6","nph6","nnl6"]
ddd["deviations"] = [[   98,   127,  171,    195,    103,   130,  201,   409],
                     [  5.1,   5.2,  7.0,   18.9,   5.7,    6.1, 9.3,   31.1]]   # spatially measured ROI in WM, added 2022-11-03
ddd["predictions"]= [[3.2,  3.2,   4.0,    8.3,   3.2,  4.6,  5.6,   14.2],
                     [4.7,  4.7,   5.2,    9.6,   5.7,  6.2,  7.2,   15.6]]   # added 2022-11-29, 90% (936,58) and 10% (10,10) 
                                                                            # per publication on SPIJN, 2nd component: myelin water
ddd["normalize"]  = true # false # true
ddd["values"]     = [[1000, 1060,  1110,  790,    940,    1040,  1090,  970],
                     [   40,   44,    44,   50,     35,     40,     44,   48]]
ddd["description"]= "relative standard deviation \n in white matter ROI (a.u.)"
collection["white_matter"] = copy(ddd)

# -------------------------------------------------------------
ddd["label"]      = "Marrow \n (tibia)"
ddd["s_label"]      = "Bone marrow"
ddd["scatter_color"]="orange"
ddd["key"]         = ["aih6","anh6","aph6","anl6","nih6","nnh6","nph6","nnl6"]
ddd["deviations"] = [[   68,    49,   54,     63,    70,    87,   96,    151],
                     [   23,    33,   33,     25,    32,    38,   40,     46]]   # on modified ROI position dd. 2023-01-19 per suggestion of team 
ddd["predictions"]= [[  2.3,   2.3,  2.6,    3.2,   3.7,   5.3,  4.9,   10.1],
                     [  1.8,   1.7,  2.0,    2.3,   3.3,   4.7,  4.5,    8.2]]   # predictions on (0.33,0.08)
ddd["normalize"]  = true # false # true
ddd["values"]     = [[  310,  290,   284,  317,    368,    404,   377,  557],
                     [  128,  147,   131,  104,    145,    151,   144,  154]]
ddd["description"]= "relative standard deviation \n in tibial bone-marrow ROI (a.u.)"
collection["marrow"] = copy(ddd)

# -------------------------------------------------------------
ddd["label"]      = "Muscle \n (calf)"
ddd["s_label"]      = "Muscle tissue"
ddd["scatter_color"]="red"
ddd["key"]        = ["aih6","anh6","aph6","anl6","nih6","nnh6","nph6","nnl6"]
ddd["deviations"] = [[   55,    68,   89,     79,    56,    94,  116,    108],
                     [  1.6,   1.5,  2.1,    2.0,   1.5,   1.7,  2.0,    2.0]]   # on modified ROI position dd. 2023-01-19 per suggestion of team
ddd["predictions"]= [[  7.4,   6.4,  8.9,    6.4,   6.7,   8.6, 11.2,    9.1],
                     [  8.0,   6.9,  8.4,    6.6,   8.2,   7.6, 10.6,    8.9]]
ddd["normalize"]  = true # false # true
ddd["values"]     = [[ 1417, 1282,  1281, 1274,   1304,   1392,  1408, 1371],
                     [   27,   23,    23,   23,     22,     23,    23,   22]]
ddd["description"]= "relative standard deviation \n in calf muscle ROI (a.u.)"
collection["muscle"] = copy(ddd)
#
for (key,ddd) in collection
    oriseq = ddd["key"]
    dev    = ddd["deviations"]
    bljDev = ddd["predictions"]

    (fig,ax)=(subplots(1,2,figsize=(12,5)))
    fig.subplots_adjust(bottom=0.1) # was 0.3
    fig.subplots_adjust(top=0.95) 
    for m in 1:2
        reorderDevs = [dev[m][findfirst(oriseq.==resequence[i])] for i in eachindex(resequence)]
        reorderBLJ = [bljDev[m][findfirst(oriseq.==resequence[i])] for i in eachindex(resequence)]

        relDev = reorderDevs./reorderDevs[refSequence]
        if ddd["normalize"]
            means  = ddd["values"]
            reorderMeans= [means[m][findfirst(oriseq.==resequence[i])] for i in eachindex(resequence)]
            relDev = reorderDevs./reorderMeans./reorderDevs[refSequence].*reorderMeans[refSequence] 
        end 
        bljRel = reorderBLJ./reorderBLJ[refSequence]

        ax[m].set_ylim(0,5.0)
        ax[m].set_xticks(1:4)
        ax[m].set_xticklabels(brief_names, fontsize=16) 
        p1 = ax[m].bar(pos, relDev, width=0.3 ,color=kleur) 
        ax[m].scatter(pos, bljRel, color="black") 
        if ddd["s_label"]=="Phantom"
            ax[m].text(2.3,4.3,"T$m",fontsize=18)
                end
    end
    ax[1].set_ylabel("relative std.", fontsize=18, horizontalalignment="left", y=0.0)
    ax[2].text(-0.3,1.4,ddd["s_label"],  fontsize=18, backgroundcolor="grey")
end


## 2023-10-20 Redesign of figure 7
# (!) dependent on fig. 6 !!!
# from 2023-01-11 "Further analysis on measured noise data in phantom, knee and brain" (messed with dd. 2023-02-03)
using Statistics

barCollection=[ ("(h)",[(1,2),(5,6)],"compCorr","SNR efficiency gain of pause"),
                ("(g)",[(1,2),(5,6)],"compare","SNR gain of pause"),
                ("(f)",[(3,1),(7,5)],"compare","SNR gain of inversion pulse"),
                ("(e)",[(4,8)],      "average","SNR gain of phase in SAR-restricted situation"),
                ("(d)",[(3,7)],      "average","SNR gain of phase with inversion"),
                ("(c)",[(2,6)],      "average","SNR gain of phase without inversion or pause"),
                ("(b)",[(1,5)],      "average","SNR gain of phase without inversion"),
                ("(a)",[(1,5),(2,6),(3,7),(4,8)],"average","average SNR gain of phase")
                ]

barColor     = ["orange","#7070ff"]
textColor    = ["black", "black"]
textValue    = ["Amplitude+Phase", "Amplitude only"]

barHeight = 0.25

figure(figsize=(6,10))
for (setNumber,(tag,setElement,task,description)) in enumerate(barCollection)
    sumRatios = zeros(2,length(setElement))
    for (i,refTestPair) in enumerate(setElement)
        (refSequence,testSequence) = refTestPair
        for (key,ddd) in collection
            oriseq = ddd["key"]
            dev    = ddd["deviations"]
            bljDev = ddd["predictions"]
            for m in 1:2
                reorderDevs = [dev[m][findfirst(oriseq.==resequence[i])] for i in eachindex(resequence)]
                reorderBLJ = [bljDev[m][findfirst(oriseq.==resequence[i])] for i in eachindex(resequence)]
                relDev = reorderDevs[testSequence]./reorderDevs[refSequence]
                if ddd["normalize"]
                    means  = ddd["values"]
                    reorderMeans= [means[m][findfirst(oriseq.==resequence[i])] for i in eachindex(resequence)]
                    relDev = reorderDevs[testSequence]./reorderMeans[testSequence]./reorderDevs[refSequence].*reorderMeans[refSequence] 
                end 
                sumRatios[1,i] += log(relDev)
                sumRatios[2,i] += log(reorderBLJ[testSequence]/reorderBLJ[refSequence])
            end
        end
    end
    sumRatios ./= 2*length(collection)
    @show sumRatios

    sumRatios .*= 100.0

    vPos = 0.0
    if task=="compare" || task=="compCorr"
        scanEff = (task=="compCorr") ? 100*log(sqrt(13.44/18.43)) : 0.0
        for i in 1:2
            verPos = setNumber+barHeight*(i-1.5)
            barh(   verPos, sumRatios[1,i]+scanEff,color=barColor[i], height=barHeight)
            # scatter(sumRatios[2,i], verPos, color="black")                    # removed 2023-02-03 per suggestion Aless
            horPos = (i==2) ? 2.0 : 2.0+sumRatios[1,1]+scanEff
            text(horPos, verPos-0.03, textValue[i], color=textColor[i])
        end
        #tPos = maximum(sumRatios)
        vPos = setNumber + barHeight + 0.05
    elseif task=="average"
        sumRatios = mean(sumRatios,dims=2)
        barh(        setNumber,            sumRatios[1,1],color="gray", height=barHeight)
        # scatter(     sumRatios[2,1],         setNumber,  color="black")          # removed 2023-02-03 per suggestion Aless
        # tPos = maximum(sumRatios,dims=1)[1,1]
        vPos = setNumber + 0.5*barHeight + 0.05
    end
    text(2.0,vPos,description)
    text(-4.0,vPos-0.2,tag)
end
plot([0,0],[1-barHeight-0.05,length(barCollection)+barHeight+0.05], color="black")
yticks([])
xlabel("SNR (or SNR efficiency) gain (in %)")
grid(axis="x")


## 2023-10-21 Refactoring of image display part 
# Taken from ScanAnalysisScript4.jl,
#      2023-02-03 Repeat of a number of figures for the phase-benefits paper 
# (taken from:) 2023-01-19 Repeat knee (change ROIs); see 2023-01-18 (messed with 2023-01-23 to show ROI-edges only)
# (Modified 2023-10-25 for the new ordering of conditions)
include("../startup.jl")
include("../scanAnalysis/ScanAnalysis.jl")
fn_common = "/smb/user/mfuderer/BLD_RT_RESEARCH_DATA/PROJECT/MRSTAT/Scandata_miha/MF_654407i/mrstat"
nRecons=10;
# !!! following line modified !!!!!
folderDetails       = ["_nnh", "_nph", "_nih", "_nnm", "_anh", "_aph", "_aih", "_anm"]   
#folderDetails       = ["_nih", "_nnh", "_nph", "_nnm", "_aih", "_anh", "_aph", "_anm"]   
sx = 224; sy = 224; maxIt = 21; 
nZoom = 1; zoomtype = 1; 
saParams = Dict()
saParams["nCases"]    = length(folderDetails)
saParams["nRecons"]   = nRecons
nSubtypes = length(folderDetails) ÷ 2;

# Read data into 7-d array
maps = ReadInto7DArray(folderDetails, nZoom, nRecons,3,maxIt,sx,sy,fn_common);

# means and standard deviations
meanoverdyn= mean(maps,dims=3);
devoverdyn = std(maps,dims=3); # case zoom [dummy] type iter x y

# Show and display
nTypesShown = 2
dispMax  = [2.0,0.28,1.0] # [2.0,0.4,1.0] 2022-11-01
#paper_description = ["A1", "A2", "A3", "A4", "B1", "B2", "B3", "B4"]
paper_description =["","","","","","","",""]  # modification 2023-02-03, removing red labels
DisplayImageSet(meanoverdyn, paper_description, nTypesShown, nSubtypes)    


## Brain 
include("../startup.jl")
include("../scanAnalysis/ScanAnalysis.jl")
fn_common = "/smb/user/mfuderer/BLD_RT_RESEARCH_DATA/PROJECT/MRSTAT/Scandata_miha/MR_631163/mrstat"
# !!! following line modified !!!!!
folderDetails       = [ "_nnh6", "_nph6", "_nih6", "_nnl6", "_anh6", "_aph6b", "_aih6", "_anl6"]
#folderDetails       = [ "_nih6", "_nnh6", "_nph6", "_nnl6", "_aih6", "_anh6", "_aph6b", "_anl6"]

maps = ReadInto7DArray(folderDetails, nZoom, nRecons,3,maxIt,sx,sy,fn_common);

meanoverdyn= mean(maps,dims=3);
devoverdyn = std(maps,dims=3); # case zoom [dummy] type iter x y
# Show and display
nTypesShown = 2
dispMax  = [2.0,0.28,1.0] # [2.0,0.4,1.0] 2022-11-01
paper_description =["","","","","","","",""]  # modification 2023-02-03, removing red labels
DisplayImageSet(meanoverdyn, paper_description, nTypesShown, nSubtypes)   

## Phantom
include("../startup.jl")
include("../scanAnalysis/ScanAnalysis.jl")
fn_common = "/smb/user/mfuderer/BLD_RT_RESEARCH_DATA/PROJECT/MRSTAT/Scandata_miha/MF_623226/mrstat"
# !!! following line modified !!!!!
folderDetails       = ["_nnh62", "_nph62",  "_nih62", "_nnl62", "_anh6b2", "_aph62", "_aih62", "_anl62" ]
#folderDetails       = ["_nih62", "_nnh62",  "_nph62", "_nnl62", "_aih62", "_anh6b2", "_aph62", "_anl62" ]

nRecons = 10; sx = 224; sy = 224; maxIt = 21; 
nRecons = 10; off=10
nZoom = 1; zoomtype = 1; 
saParams = Dict()
saParams["nCases"]    = length(folderDetails)
saParams["nRecons"]   = nRecons
nSubtypes = length(folderDetails) ÷ 2;

# Read data into 7-d array
maps = ReadInto7DArray(folderDetails, nZoom, nRecons,3,maxIt,sx,sy,fn_common);

# means and standard deviations
#maps = reverse(maps, dims=6)
meanoverdyn= mean(maps,dims=3);
devoverdyn = std(maps,dims=3); # case zoom [dummy] type iter x y
# display the sets and tube-average correspondence
paper_description =["","","","","","","",""]  # modification 2023-02-03, removing red labels
DisplayImageSet(meanoverdyn, paper_description, nTypesShown, nSubtypes)


## 2023-10-23 Checking whether 20220818W reproduces 
include("setup.jl")
recon_options = Dict() # erase all existing settings
nsweeps = 6
nky = 224
nTR = round(Int64, nsweeps*nky)
ky = 1.0 .*(repeat(1:nky, inner=1, outer=nsweeps)); 
kz = ones(nsweeps*nky)
trajectorySet = BLAKJac.TrajectorySet(ky,kz)
BLAKJac.BLAKJac_defaults!(trajectorySet, recon_options)

recon_options["useSymmetry"] = true     
recon_options["TR"]      = 0.01
recon_options["startstate"] = -1 # 1 # -1 S 
recon_options["sigma_ref"] = 1.4 #  See logbook 20220815
time_limit                 = 3000 # 3000 # 5000 N # 2500 R 
recon_options["optpars"]   = Optim.Options(time_limit = time_limit, iterations = 100000, f_tol=1.0e-10, g_tol = 1.0e-10)  
recon_options["opt_criterion"] = "sar-limited noise"  
recon_options["sar_limit"]  = 40^2/0.01 
recon_options["emphasize_low_freq"] = true 
recon_options["opt_focus"] = "max"      
recon_options["opt_initialize"] = "cRandom30" 
recon_options["opt_complex"] = false   
recon_options["opt_account_maxFlip"] = true 
recon_options["opt_keep_positive"] = false                           
recon_options["opt_slow_phase"] =  false
recon_options["considerCyclic"] = false # true # false S   
recon_options["TW"] = 0.0             # 0  
recon_options["opt_emergeCriterion"] = 2000
ph = [] 
#ph = zeros(nTR); ph .= 2.0;
recon_options["opt_imposed_2nd_derivative_of_phase"] = ph
recon_options["opt_iterations_limit"] = 1 #1 # 50 R 
recon_options["sizeSteps"] = [6]   # [6] # [6,10,15,20,30] R  
#recon_options["portion_size"]=38 # only R
nRealizations = 6 # 6 # 1 R 

fn_base = "20231023S"
for i in 1:nRealizations
    stageText = ""
    portionRange = 0:0
    fn = "/home/mfuderer/Documents/Julia/Capture/$fn_base($i).jld2"
    RFdeg = BLAKJac.BLAKJac_optimize(trajectorySet, recon_options, i);
    FileIO.save(fn,"RFdeg",RFdeg)
    end     

##
fn_base = "20220818W"; nRealizations=6
#fn_base = "20231023S"; nRealizations=6
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


## 2023-10-23 Observed that the above reproduces; streamlining the script 
include("setup.jl")
recon_options = Dict() # erase all existing settings
nsweeps = 6
nky = 224
nTR = round(Int64, nsweeps*nky)
ky = 1.0 .*(repeat(1:nky, inner=1, outer=nsweeps)); 
kz = ones(nsweeps*nky)
trajectorySet = BLAKJac.TrajectorySet(ky,kz)
BLAKJac.BLAKJac_defaults!(trajectorySet, recon_options)

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
    fn = "/home/mfuderer/Documents/Julia/Capture/$fn_base($i).jld2"
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



## 2023-10-24 Further refactoring of what was once Fig.6: plot SNR rather than noise 
#     (built upon 2023-10-20 Fig. 6 refactoring, which came 
#      from 2023-02-08 (...) last modified 2023-04-09)
using Statistics
PyPlot.rc("font", family="serif")

mapnames = ["std(T1)", "std(T2)", "rho"];
pos  =[  1.0,   2.0,    3.0,   4.0,   1.3,  2.3,  3.3,  4.3]
kleur=["orange","orange","orange","orange", "blue","blue","blue","blue"]
brief_names = ["Base","No pause","Inv.","Low SAR"]  
resequence  = ["anh6","aph6","aih6","anl6","nnh6","nph6","nih6","nnl6"]
refSequence = 5
collection = Dict()
ddd = Dict()

ddd["label"]      = "Phantom \n"
ddd["s_label"]    = "Phantom"
ddd["scatter_color"]="black"
ddd["key"]        = ["aih6","anh6","aph6","anl6","nih6","nnh6","nph6","nnl6"]
ddd["deviations"] = [[ 14.1,  17.1,   18.9,  39.3,  18.0, 62.1, 49.0,  129],
                       [2.3,   2.4,    2.7,   8.1,   3.6,  6.0,  4.3,  18.2]] # phantom
ddd["predictions"]= [[  2.4,   2.5,   3.0,    5.0,   3.3,  4.9,  5.2,   16.1],
                     [  2.4,   2.4,   2.7,    5.0,   3.3,  4.3,  4.5,   16.5]]   # 7-points BLAKJac
ddd["normalize"]  = false
ddd["values"]     = []
ddd["description"]= "standard deviation \n in phantom vials (a.u.)"
collection["phantom_N"] = copy(ddd)

# -------------------------------------------------------------
ddd["label"]      = "White \n matter"
ddd["s_label"]      = "White matter"
ddd["scatter_color"]="green"
ddd["key"]          = ["aih6","anh6","aph6","anl6","nih6","nnh6","nph6","nnl6"]
ddd["deviations"] = [[   98,   127,  171,    195,    103,   130,  201,   409],
                     [  5.1,   5.2,  7.0,   18.9,   5.7,    6.1, 9.3,   31.1]]   # spatially measured ROI in WM, added 2022-11-03
ddd["predictions"]= [[3.2,  3.2,   4.0,    8.3,   3.2,  4.6,  5.6,   14.2],
                     [4.7,  4.7,   5.2,    9.6,   5.7,  6.2,  7.2,   15.6]]   # added 2022-11-29, 90% (936,58) and 10% (10,10) 
                                                                            # per publication on SPIJN, 2nd component: myelin water
ddd["normalize"]  = true # false # true
ddd["values"]     = [[1000, 1060,  1110,  790,    940,    1040,  1090,  970],
                     [   40,   44,    44,   50,     35,     40,     44,   48]]
ddd["description"]= "relative standard deviation \n in white matter ROI (a.u.)"
collection["white_matter"] = copy(ddd)

# -------------------------------------------------------------
ddd["label"]      = "Marrow \n (tibia)"
ddd["s_label"]      = "Bone marrow"
ddd["scatter_color"]="orange"
ddd["key"]         = ["aih6","anh6","aph6","anl6","nih6","nnh6","nph6","nnl6"]
ddd["deviations"] = [[   68,    49,   54,     63,    70,    87,   96,    151],
                     [   23,    33,   33,     25,    32,    38,   40,     46]]   # on modified ROI position dd. 2023-01-19 per suggestion of team 
ddd["predictions"]= [[  2.3,   2.3,  2.6,    3.2,   3.7,   5.3,  4.9,   10.1],
                     [  1.8,   1.7,  2.0,    2.3,   3.3,   4.7,  4.5,    8.2]]   # predictions on (0.33,0.08)
ddd["normalize"]  = true # false # true
ddd["values"]     = [[  310,  290,   284,  317,    368,    404,   377,  557],
                     [  128,  147,   131,  104,    145,    151,   144,  154]]
ddd["description"]= "relative standard deviation \n in tibial bone-marrow ROI (a.u.)"
collection["marrow"] = copy(ddd)

# -------------------------------------------------------------
ddd["label"]      = "Muscle \n (calf)"
ddd["s_label"]      = "Muscle tissue"
ddd["scatter_color"]="red"
ddd["key"]        = ["aih6","anh6","aph6","anl6","nih6","nnh6","nph6","nnl6"]
ddd["deviations"] = [[   55,    68,   89,     79,    56,    94,  116,    108],
                     [  1.6,   1.5,  2.1,    2.0,   1.5,   1.7,  2.0,    2.0]]   # on modified ROI position dd. 2023-01-19 per suggestion of team
ddd["predictions"]= [[  7.4,   6.4,  8.9,    6.4,   6.7,   8.6, 11.2,    9.1],
                     [  8.0,   6.9,  8.4,    6.6,   8.2,   7.6, 10.6,    8.9]]
ddd["normalize"]  = true # false # true
ddd["values"]     = [[ 1417, 1282,  1281, 1274,   1304,   1392,  1408, 1371],
                     [   27,   23,    23,   23,     22,     23,    23,   22]]
ddd["description"]= "relative standard deviation \n in calf muscle ROI (a.u.)"
collection["muscle"] = copy(ddd)
#
for (key,ddd) in collection
    oriseq = ddd["key"]
    dev    = ddd["deviations"]
    bljDev = ddd["predictions"]

    (fig,ax)=(subplots(1,2,figsize=(12,5)))
    fig.subplots_adjust(bottom=0.1) # was 0.3
    fig.subplots_adjust(top=0.95) 
    for m in 1:2
        reorderDevs = [dev[m][findfirst(oriseq.==resequence[i])] for i in eachindex(resequence)]
        reorderBLJ = [bljDev[m][findfirst(oriseq.==resequence[i])] for i in eachindex(resequence)]

        relDev = reorderDevs./reorderDevs[refSequence]
        if ddd["normalize"]
            means  = ddd["values"]
            reorderMeans= [means[m][findfirst(oriseq.==resequence[i])] for i in eachindex(resequence)]
            relDev = reorderDevs./reorderMeans./reorderDevs[refSequence].*reorderMeans[refSequence] 
        end 
        bljRel = reorderBLJ./reorderBLJ[refSequence]

        ax[m].set_ylim(0,5.0)
        ax[m].set_xticks(1:4)
        ax[m].set_xticklabels(brief_names, fontsize=16) 
        p1 = ax[m].bar(pos, relDev.^-1, width=0.3 ,color=kleur)  #inverted
        ax[m].scatter(pos, bljRel.^-1, color="black")            # inverted
        if ddd["s_label"]=="Phantom"
            ax[m].text(2.3,4.3,"T$m",fontsize=18)
        end
    end
    ax[1].set_ylabel("relative SNR", fontsize=18, horizontalalignment="left", y=0.0)
    ax[2].text(-0.3,1.4,ddd["s_label"],  fontsize=18, backgroundcolor="grey")
end


## 2023-10-27 Re-display of some DOF-results for figure S-&
include("setup.jl")
recon_options = Dict() # erase all existing settings
nsweeps = 6
nky = 224
nTR = round(Int64, nsweeps*nky)
ky = 1.0 .*(repeat(1:nky, inner=1, outer=nsweeps)); 
kz = ones(nsweeps*nky)
trajectorySet = BLAKJac.TrajectorySet(ky,kz)
BLAKJac.BLAKJac_defaults!(trajectorySet, recon_options)

recon_options["useSymmetry"] = true     
recon_options["TR"]      = 0.01
recon_options["startstate"] = 1
recon_options["sigma_ref"] = 1.4 #  See logbook 20220815
recon_options["emphasize_low_freq"] = true 
recon_options["considerCyclic"] = true  
recon_options["TW"] = 0.0             # 0  

#
#fn_base = "20220920A"; nRealizations=1; skip = 0       # DOF 6 phase
#fn_base = "20231018R"; nRealizations=1; skip = 0      # DOF full no phase
fn_base = "20231013X"; nRealizations=1; skip=2      # DOF 6 no phase
scores=zeros(nRealizations,2)
saved_H = Dict()
recon_options["rfName"]  = "from_file"
recon_options["rfFunction"] = rfDictionary[recon_options["rfName"]]
rfFunction = recon_options["rfFunction"]
#(fig,ax)=(subplots(Int(ceil(nRealizations/3)),3,figsize=(9,3)))
figure(figsize=(3,3)) # modified 2023-10-27 !!!!!!!!!!
for i in 1:nRealizations
    fn = "$fn_base($(i+skip))"
    recon_options["rfFile"]  = fn
    RFdeg = rfFunction(recon_options["nTR"], recon_options["nky"])
    plot(abs.(RFdeg))
    anglesdd = zeros(length(RFdeg))
    for i in 1:length(RFdeg)-2
        anglesdd[i] = (rad2deg(angle(conj(RFdeg[i])*RFdeg[i+1]*RFdeg[i+1]*conj(RFdeg[i+2])))+270.0) % 180.0 -90.0
    end
    plot(anglesdd)
    noises, ItotAll = BLAKJac.BLAKJac_analysis!(cpu, RFdeg, trajectorySet, recon_options, saved_H) 
    @show noises, mean(noises), ItotAll
    scores[i,1] = noises[3] # maximum(noises[2:3])
    scores[i,2] = ItotAll
end










