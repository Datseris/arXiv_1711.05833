#=
This code plots the resistivity curves
used in our paper https://arxiv.org/abs/1711.05833.

The user can optionally choose to include pinned orbits or not in the resistivity
calculations. `keep_pinned = false` does not keep the pinned orbits.
It integrates the chaotic correlations and then divides the resulting
resistivity curves with g_c. If one wants to include the pinned orbits, then with
`keep_pinned=true` correlation functions are transformed as follows:

C = C_chaos*g_c + 0.5*(1-g_c)*cos(2Bt)  [or sin()]

In this case the end result is (of course) not divided by g_c.
=#
using FileIO

include("Library.jl")
cd(@__DIR__)

## Potential type
pt = :PSB # :fkg, :br, :std, :PSB, :RPSB

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

dt = 0.01; T = 5000.0
t = 0:dt:T
N = 5000

par =
Dict(:pt => pt, :B => B, :d0 => d0, :c=>c,
:β=>β, :ε =>ε, :dt => dt, :T=>T, :N=>N)

# Resistivity Parameters
τ = [2.5, 5.0, 20.0]  #impurity scattering time in intrinsic units
keep_pinned = false #, false, "nogc"]

using PyPlot
figure()
for t in τ

    B, Rxx, Rxy, gc = generate_res(par, t, keep_pinned)

    plot(B, Rxy, label = "τ = $τ")
end

xlabel("\$B\$"); ylabel("\$R_{xy}\$")
