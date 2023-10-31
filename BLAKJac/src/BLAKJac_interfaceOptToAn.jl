# Interface functions between BLAKJac optimization and underlying BLAKJac analysis functions
# Added 2022-07-11


# Since the optimizer wants to interface to a function accepting only floats, but optimization of complex values may be needed,
# this function converts a complex vector into a real vector of doubled length
function ConvertCpxToConcatenatedFloat(RFC::Vector{ComplexF64})
    length = size(RFC,1)
    RF = zeros(length*2)
    for i=1:length
        RF[i]        = real(RFC[i])
        RF[i+length] = imag(RFC[i])
    end
    return RF
end

# If the vector-to-be-optimized is already non-complex float, the "Unwrap" is a dummy function
function ConvertCpxToConcatenatedFloat(v::Vector{Float64}) return v; end

# "Wrapping" is the reverse function: converting a double-length float vector into a complex vector
function ConvertConcatenatedFloatToCpx(RF::Vector{Float64})
    halfsize = size(RF,1)÷2
    RFC = zeros(ComplexF64,halfsize)
    for i=1:halfsize RFC[i] = complex(RF[i],RF[i+halfsize]) end
    return RFC
end

# Interface between (A) an optimizer that wants to optimize on a limited set FLOAT of values and (B) an analyzer that expects a full-length COMPLEX array
function ExpandCpxAndAnalyze(RF::Vector{Float64}, trajectorySet::Vector{Vector{TrajectoryElement}}, options::Dict)
    # convert array of floats into complex
    RFC = ConvertConcatenatedFloatToCpx(RF) # Combine float-array of double length into omplex array
    return ExpandAndAnalyze(RFC, trajectorySet, options)
end

# Interface between (A) an optimizer that wants to optimize on a limited set of values and (B) an analyzer that expects a full-length array
function ExpandAndAnalyze(RFshort::Array, trajectorySet::Vector{Vector{TrajectoryElement}}, options::Dict)
    # interpolate short RF array to full length
    nTR = length(trajectorySet)
    RFdeg = ExpandRF(RFshort, nTR, options)
    if ndims(RFshort)==1
        # The logic here is: if RFshort has two columns, then the second serves for a 2nd derivative of the phase;
        #   in that case it is pointless to expect complex values in the first column, which are the flip angles
        RFdeg = complex(RFdeg)  # convert to complex - if not so already
    else
        ; # no activity; the two columns of the array should be and remain real
    end
    penalty =  BLAKJac_criterion(RFdeg, trajectorySet, options)
    return penalty
end


# Interfaces between (A) an optimizer that only wants to optimize on a FLOAT portion of a vector and (B) an analyzer that considers the full COMPLEX length
function PlugInCpxAndAnalyze(RFpart::Vector{Float64}, portionRange, RFdeg::Vector{ComplexF64}, trajectorySet::Vector{Vector{TrajectoryElement}}, options::Dict)
    RFC = ConvertConcatenatedFloatToCpx(RFpart) # Combine float-array of double length into omplex array
    return PlugInAndAnalyze(RFC, portionRange, RFdeg, trajectorySet, options)
end

# Interfaces between (A) an optimizer that only wants to optimize on a portion of a vector and (B) an analyzer that considers the full length
function PlugInAndAnalyze(RFpart::Vector{T}, portionRange, RFdeg::Vector{T}, trajectorySet::Vector{Vector{TrajectoryElement}}, options::Dict) where T
    # plug the part into the full array
    RFdeg[portionRange] = RFpart
    RFdegC = complex(RFdeg)
    penalty =  BLAKJac_criterion(RFdegC, trajectorySet, options)
    return penalty
end

# Interpolates to full length.
# If there is a second column, this is interpreted as the 2nd derivative of the phase, which may be either to-be-optimized or preset
function ExpandRF(RFin::Array, nTR, options)
    ph  = options["opt_imposed_2nd_derivative_of_phase"]

    if (size(RFin,2)==1)
        RFall = IncreaseResolution(RFin, nTR)
        return RFall
    else
        # assumed to be one column for amplitude, other column for 2nd derivative of phase
        RFdegC= zeros(ComplexF64,nTR)
        RFall = IncreaseResolution(RFin, nTR)

        # if provided, overrule 'optimized' phase by predetermined (2nd derivative of) phase
        if (length(ph) > 0) RFall[:,2] .= ph; end;

        # double-integrate second column and add it as a phase
        dp=0.0; p = 0.0
        for i in 1:nTR
            dp+=deg2rad(RFall[i,2]);
            p+=dp;
            RFdegC[i] = RFall[i,1]*exp(im * p)
        end
        return RFdegC
    end
end

# For a given complex vector of RF angles and a trajectory, it outputs a criterion value, according to the selected option.
# It may also, conditionally, plot an output
function BLAKJac_criterion(RFdegC::Vector{ComplexF64}, trajectorySet::Vector{Vector{TrajectoryElement}}, options::Dict)

    # --------------------------------------------------
    # Analyze the sequence
    cpu  = ComputationalResources.CPU1()
    noises, infocon, b1s = BLAKJac_analysis!(cpu, RFdegC, trajectorySet, options)

    # --------------------------------------------------
    # Select the criterion value ("out") given the options
    noisePenalties = Dict();
    b1Penalties    = Dict();
    noisePenalties["rho"] = noises[1]
    noisePenalties["T1"]  = noises[2]
    noisePenalties["T2"]  = noises[3]
    noisePenalties["mean"]= mean(noises)
    noisePenalties["max"] = maximum(noises[2:3])
    noisePenalties["weighted"]=noises[2]/options["T1ref"]+noises[3]/options["T2ref"];
    b1Penalties["mean"] = mean(abs.(b1s))
    b1Penalties["max"]  = maximum(abs.(b1s[2:3]))
    b1Penalties["T1"]   = abs(b1s[2])
    b1Penalties["T2"]   = abs(b1s[3])    

    noisePenalty = noisePenalties[options["opt_focus"]]
    b1Penalty    = b1Penalties[options["opt_focus"]]

    # A factor can be taken into account expressing that we gain in sampling time if time needed for RF pulses is small
    maxFlip = maximum(abs.(RFdegC))
     lowFlipFactor = 1.0/sqrt(1.57 - 0.57 * min(maxFlip/180.0, 2.0)) # see logbook dd. 2021-05-25  # modified 2021-12-22
    nplff = lowFlipFactor*noisePenalty
    sarLevel = mean(abs2.(RFdegC))/options["TR"];
    sarPenalty = (sarLevel/options["sar_limit"])^10;

    if options["opt_criterion"] == "noise_level"
        out = noisePenalty
    elseif options["opt_criterion"] == "low-flip corrected noise"
        out = nplff
    elseif options["opt_criterion"] == "sar-limited noise"
        out = noisePenalty+sarPenalty
    elseif options["opt_criterion"] == "information content"
        out = -0.01*infocon  # negative, since the optimizer tries to minimize the output
                            # the multiplication by 0.01 serves to bring the figure to the same order of magnitude as the noise values,
                            # since the NelderMead optimizer applies absolute tolerance criteria
    elseif options["opt_criterion"] == "information content penalized"
        out = -0.01*infocon*sqrt(1.0/(1.3*lowFlipFactor)) # very houtje-touwtje; the 1.3 has been chosen to make a max of approx 80 deg 'neutral' 
    elseif options["opt_criterion"] == "B1 sensitivity sar-limited"
        out = b1Penalty+sarPenalty
    elseif options["opt_criterion"] == "B1 sensitivity/noise mix, sar-limited"
        out = 0.1*b1Penalty+sarPenalty+noisePenalty
    else
        throw("unknown optimization criterion")
    end
    negFlip = 0.0; imFlip = 0.0;
    if options["opt_keep_positive"]
        negFlip = maximum(max.(-real.(RFdegC),0.0))
        out += negFlip;
    end

    # --------------------------------------------------
    # Occasionally, display intermediate output
    options["optcount"] +=1
    i = options["optcount"]
    emergeCriterion = options["opt_emergeCriterion"] # (1*1000^3) ÷ (length(RFdegC)^2)
    emerge = (i % emergeCriterion == 1)
    if emerge
        options["plotfuncs"]["close"]()
        options["plotfuncs"]["intermediate"](RFdegC, options)
    end
    if (i % (emergeCriterion÷5+1) == 1)
        stageText = options["stage_text"]
        @show noisePenalty, nplff, infocon, out, stageText
        @show noisePenalty, sarPenalty, sarLevel, negFlip, out
    end

    return out
end

# Interfaces between (A) an outside world that wants to see a potentially complex array optimized using specific options
#                and (B) an optimizer that only wants to see a function accepting a float array and nothing else
function WrappedPortionOptimize(RFpart::Vector{ComplexF64}, portionRange, RFdeg::Vector{ComplexF64}, trajectorySet::Vector{Vector{TrajectoryElement}}, options::Dict)
    RFfloat = ConvertCpxToConcatenatedFloat(RFpart)
    opt_method = options["opt_method"]
    optpars    = options["optpars"]

    PlugInCpxAndAnalyze_(y) = PlugInCpxAndAnalyze(y, portionRange, RFdeg, trajectorySet, options)
    res = optimize(PlugInCpxAndAnalyze_, RFfloat, opt_method(), optpars)
    RFportion = ConvertConcatenatedFloatToCpx(Optim.minimizer(res))
    return res, RFportion
end

function WrappedPortionOptimize(RFpart::Vector{Float64}, portionRange, RFdeg::Vector{Float64}, trajectorySet::Vector{Vector{TrajectoryElement}}, options::Dict)
    opt_method = options["opt_method"]
    optpars    = options["optpars"]
    PlugInAndAnalyze_(y) = PlugInAndAnalyze(y, portionRange, RFdeg, trajectorySet, options)
    res = optimize(PlugInAndAnalyze_, RFpart, opt_method(), optpars)
    return res, Optim.minimizer(res)
end

function WrappedLowResOptimize(RFshort::Vector{ComplexF64}, trajectorySet::Vector{Vector{TrajectoryElement}}, options::Dict)
    RFfloat = ConvertCpxToConcatenatedFloat(RFshort) # Split complex array into float-array of double length
    opt_method = options["opt_method"]
    optpars    = options["optpars"]
    ExpandCpxAndAnalyze_(RF) = ExpandCpxAndAnalyze(RF, trajectorySet, options)
    res = optimize(ExpandCpxAndAnalyze_, RFfloat, opt_method(), optpars)
    RFold = ConvertConcatenatedFloatToCpx(Optim.minimizer(res))
    return res, RFold
end

function WrappedLowResOptimize(RFshort::Array{Float64}, trajectorySet::Vector{Vector{TrajectoryElement}}, options::Dict)
    opt_method = options["opt_method"]
    optpars    = options["optpars"]
    ExpandAndAnalyze_(y) = ExpandAndAnalyze(y, trajectorySet, options)
    res = optimize(ExpandAndAnalyze_, RFshort, opt_method(), optpars)
    return res, Optim.minimizer(res)
end



