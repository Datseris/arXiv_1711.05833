#=
This file simply plots the potential and an example orbit on top of it,
using the model from our paper https://arxiv.org/abs/1711.05833
=#

using OrdinaryDiffEq, PyPlot, FileIO

cd(@__DIR__)
include("Library.jl")

# Potential type
pt = :br # :fkg, :br, :std

# Fundamental parameters:
B = 0.1
d0 = 0.3   #antidot diameter
c = 0.1  #cutoff distance, must be 0 < c < 0.5-d0/2 !!!

# FKG (old)
β = 4

# Boundary Roughness parameters (for ε = 0 no BR takes place)
ε = 0.1  #relative edge roughness
M = 16   #number of sine modes. complexity of edge roughness
matrix_size = 100  #size of the matrix used for the random coefficients
matrixnumber = 1

dt = 0.01; T = 100.0
t = 0:dt:T

parameters =
Dict(:pt => pt, :d0 => d0, :c=>c, :β=>β, :ε =>ε, :dt => dt, :T=>T, :N =>1)

# Set up ODEProblem
prob, potentialenergy = setup_ODEProblem(parameters)

# random initial condition:
var0 = initial(:ultra, potentialenergy)
prob.u0 .= var0

# solve:
sol = solve(prob, Tsit5(); abstol = 1e-9, reltol = 1e-9, saveat = t,
force_dtmin=true)

# Extract the variables from the solution
x = sol[1, :]
y = sol[2, :]

plot_orbit(x, y, potentialenergy)
