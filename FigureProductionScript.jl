## Figures S-5, S-6 and S-7
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



## Figure S-7 inset subplots
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

recon_options["useSymmetry"] = true     
recon_options["TR"]      = 0.01
recon_options["startstate"] = 1
recon_options["sigma_ref"] = 1.4 
recon_options["emphasize_low_freq"] = true 
recon_options["considerCyclic"] = true  
recon_options["TW"] = 0.0             # 0  

#
#fn_base = "20220920A"; nRealizations=1; skip = 0       # DOF 6 phase
#fn_base = "20231018R"; nRealizations=1; skip = 0      # DOF full no phase
fn_base = "20231013X"; nRealizations=1; skip=2      # DOF 6 no phase
scores=zeros(nRealizations,2)
saved_H = Dict()

figure(figsize=(3,3)) 
for i in 1:nRealizations
    fn = "$fn_base($(i+skip))"
    fPath = dataFolder*"/"*fn*".jld2"
    vars = FileIO.load(fPath)
    RFdeg = vars["RFdeg"]
    RFdeg = vec(complex.(RFdeg))
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


## Figure S-8
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
recon_options["sigma_ref"] = 1.4 
recon_options["emphasize_low_freq"] = false            # this is for optimization; for evaluation, use it without
recon_options["startstate"]     = 1
recon_options["considerCyclic"] = true
fn_base =  "20231024aph"
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




## Elements of Figure 3
using Statistics
PyPlot.rc("font", family="serif")

mapnames = ["std(T1)", "std(T2)", "rho"];
pos  =[  1.3,   2.3,    3.3,   4.3,   1.0,  2.0,  3.0,  4.0]
kleur=["orange","orange","orange","orange", "blue","blue","blue","blue"]
labelText=["Ampl+Phase","Ampl+Phase","Ampl+Phase","Ampl+Phase","Ampl-only","Ampl-only","Ampl-only","Ampl-only"]
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
    fig.subplots_adjust(wspace=0.01,hspace=0.01,left=0.1,right=0.97,bottom=0.1,top=0.95)

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
            ax[2].text(0.86,1.5,"A. Amplitude-only",fontsize=18,rotation=90,color="blue")
            ax[2].text(1.16,2.8,"B. Ampl+Phase",fontsize=18,rotation=90,color="orange")
        end
    end
    ax[1].set_ylabel("relative SNR", fontsize=18, horizontalalignment="left", y=0.0)
    ax[2].text(-0.1,1.8,ddd["s_label"],  fontsize=18, backgroundcolor="grey")
end


## Figure 4 (depends on the above)
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



## Figure 5 (depends on script for fig 3 being run)
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





## Figure 1
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
dataFolder = "/home/mfuderer/Documents/Julia/Capture"
recon_options["emphasize_low_freq"] = false # Statement (dummy), to emphasize difference to optimization results

recon_options["useSymmetry"] = true     
recon_options["TR"]      = 0.01
recon_options["handleB1"] =  "no" 
recon_options["rfName"]  = "from_file"

caseSet = [ ("20231024nnh(1)","nnh6" ,      1,false,"A-Base",   "No inversion (1)"),
            ("20231024anh(2)","anh6",       1,false,"B-Base",   "No inversion (1)"),
            ("20231024nph(3)","nph6",       1,true, "A-NoPause",   "No pause (2)"),
            ("20231024aph(1)","aph6",       1,true ,"B-NoPause",   "No pause (2)"),
            ("20231024nih(5)","nih6_818W5",-1,false,"A-Invert","Initial inversion (Def)"),
            ("20231024aih(5)","aih6_815M5",-1,false,"B-Invert","Initial inversion (Def)"),
            ("20231024nnl(3)","nnl6_815S3", 1,false,"A-LowSAR",   "Low SAR (3)"),
            ("20231024anl(4)","anl6_816T4", 1,false,"B-LowSAR",   "Low SAR (3)") ]
saved_H = Dict()
(fig,ax)=(subplots(Int(ceil(length(caseSet)/4)),4,figsize=(14,10)))

for ccc in eachindex(caseSet)
    (case,name,inversion,cyclic,label,axTxt) = caseSet[ccc]
    recon_options["considerCyclic"] = cyclic
    recon_options["startstate"] = inversion
    fPath = dataFolder*"/"*case*".jld2"
    vars = FileIO.load(fPath)
    RFdeg = vars["RFdeg"]
    RFdeg = vec(complex.(RFdeg))
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



















