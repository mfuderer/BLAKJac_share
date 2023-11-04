# BLAKJac_share
Julia language script examples for the use of BLAKJac sequence-optimization.
More specifically, two scripts are provided:
* SequenceGenerationScript: for the generation of RF-pulse sequences as shown in figure 1 of a manuscript on the benefits of applying an optimized phase to RF pulses. This script uses BLAKJac_optimize, which in turn uses BLAKJac_analysis.
* FigureProductionScript (Not directly using BLAKJac) the formatting of figures used in aforementioned manuscript.

Required packages are:
BlochSimulators, Colors, ComputationalResources, CUDA, DelimitedFiles, Distributed, DistributedArrays, FFTW, FileIO, Interpolations, LinearAlgebra, Optim, Printf, PyCall, PyPlot, Random, StaticArrays, Statistics, StructArrays

