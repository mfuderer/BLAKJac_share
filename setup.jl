# In Visual Studio Code, make sure the Julia extension is installed and
# start a Julia REPL with Ctrl + Shift + P > Julia: start REPL

## Define packages
using Pkg
Pkg.activate(".")

using LinearAlgebra # package for linear algebra operations (dot, )
using StaticArrays # package to create stack-allocated vectors/matrices/arrays

using Random       # Package for generating random numbers (comes with Julia)
using DelimitedFiles # For input of .txt files containing RF-patterns
using FFTW
using FileIO
using Interpolations

using PyPlot # same plotting backend as commonly used in python
using PyCall
using ComputationalResources
using StructArrays
using Optim
using Interpolations
using Printf
using BlochSimulators
using Colors
using Statistics

cpu  = ComputationalResources.CPU1()

include("BLAKJac/src/recon_options.jl")
include("BLAKJac/src/RF_Shapes.jl")
include("BLAKJac/src/BLAKJac.jl")

using .BLAKJac 