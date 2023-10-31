module BLAKJac

    # Load external packages
    using StaticArrays
    using ComputationalResources
    using Distributed
    using DistributedArrays
    using StructArrays
    using DelimitedFiles
    using CUDA
    using Printf
    using LinearAlgebra # package for linear algebra operations (dot, )
    using Statistics
    using Optim
    using Interpolations
    using FileIO
    using Random       # Package for generating random numbers (comes with Julia)
    # using PyPlot

    # Load self-written packages
    #using BlochSimulators
    #using ..BlochDerivatives

    include("BLAKJac_analysis_structs.jl")
    include("BLAKJac_analysis.jl")
    include("BLAKJac_defaults.jl")
    include("BLAKJac_interfaceOptToAn.jl")
    include("BLAKJac_optimize.jl")
    include("PlotFunctions.jl")
    cpu  = ComputationalResources.CPU1()

    export BLAKJac_analysis!, BLAKJac_criterion, WrappedLowResOptimize, WrappedPortionOptimize, BLAKJac_optimize, BLAKJac_defaults!, TrajectorySet

end

