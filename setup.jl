# In Visual Studio Code, make sure the Julia extension is installed and
# start a Julia REPL with Ctrl + Shift + P > Julia: start REPL

## Define packages
using Revise # package to automatically load changes you make in your own modules
             # or in files that you "includet" (include & track)
using Pkg
Pkg.activate(".")

using LinearAlgebra # package for linear algebra operations (dot, )
using StaticArrays # package to create stack-allocated vectors/matrices/arrays

using Random       # Package for generating random numbers (comes with Julia)
using DelimitedFiles # For input of .txt files containing RF-patterns
using FFTW
using FileIO
using Interpolations

# Before adding PyPlot, do
#     ENV["PYTHON"] = "" in the Julia REPL.
#     Add PyPlot in Pkg mode
#     Pkg.build("PyCall")
#     In vscode, deactivate settings > Extensions > Julia > Use Plot Plane
#     Restart Julia
using PyPlot # same plotting backend as commonly used in python
using PyCall
#using ImageView
using ComputationalResources
using StructArrays
using Optim
using Interpolations
using Printf

cpu  = ComputationalResources.CPU1()
#using MRSTAT.BlochSimulators
#using MRSTAT
# using UMCUtils
using Colors
using Statistics

include("BLAKJac/src/recon_options.jl")
include("BLAKJac/src/RF_Shapes.jl")
include("BLAKJac/src/tissueparameters.jl")
include("BLAKJac/src/fisp3d.jl")
include("BLAKJac/src/BLAKJac.jl")
#include("../numerical_phantom/k_shapes.jl")
#if (!@isdefined(gelphantoms)) include("gelphantoms.jl") end
#if (!@isdefined(tissues)) include("tissues.jl") end
