#=
This file creates the velocity correlation functions of the
billiard systems used
in our paper https://arxiv.org/abs/1711.05833


Notice that one has to modify the lines:
ttt, poss, vels, omeg = evolve!(randominside(bt, 2B), bt, T)
xt, yt, vxt, vyt, tim = construct(ttt, poss, vels, omeg, dt)

to:
ttt, poss, vels = evolve!(randominside(bt), bt, T)
xt, yt, vxt, vyt, tim = construct(ttt, poss, vels, dt)

to run the program with B=0 (zero magnetic field)
=#
using DynamicalBilliards, Interpolations, FileIO

cd(@__DIR__)
include("Library.jl")

B = isempty(args) ? 0.01 : parse(ARGS[1])

######################################
# Configure system                   #
######################################
function configure_system()
    # Potential type
    pt = :PSB # :PSB or :RPSB

    # Fundamental parameters:
    d0 = 0.3   #antidot diameter

    dt = 0.01; T = 2500.0
    t = 0:dt:T
    N = 2500

    par =
    Dict(:pt => pt, :d0 => d0, :dt => dt, :T=>T, :N=>N)

    return par, billiard_table(par)
end

######################################
# main                   #
######################################
function main(B)

par, bt = configure_system()
# Numerical Parameters
T = par[:T]
N = par[:N]
dt = par[:dt]
t = 0.0:dt:T; nt = length(t)

# Paths and names:
spath = "data/correlations/"
output = jldname(par, "Cor"; B=B)

println("\n----------------------------------")
println("Correlation Functions of Billiards.")
println("Parameters:")
for el in par
    println(el[1], " = ", el[2])
end
println("B = $B")

# Set up correlations:
print("\nPlanning FFT. . .")
# initialize fourier plans
rplan, iplan = prepare_fourier(nt)
# Initialize average cor
avrg_cxx = zeros(Float64, nt)
avrg_cxy = zeros(Float64, nt)
print("  Done.\n"); flush(STDOUT)

## COMPUTATION ##
pinned = 0
j = 0
perc = 0
println("Percentage:")
while j  < N
    println("$(j+1) / $(N)")

    ttt, poss, vels, omeg = evolve!(randominside(bt, 2B), bt, T)

    if ttt[end] == Inf
        pinned += 1
        continue
    end

    xt, yt, vxt, vyt, tim = construct(ttt, poss, vels, omeg, dt)
    cxx, cxy = velcor(t, tim, vxt, vyt, rplan, iplan)
    avrg_cxx .+= cxx; avrg_cxy .+= cxy

    j += 1
end
# Correlations are averaged only on chaotic phase-space.
# They DO NOT include the gc factor!
avrg_cxx ./= N; avrg_cxy ./= N;
g_c = 1.0 - (pinned)/(N + pinned)

ispath(spath) || mkpath(spath)
save(spath*output, Dict("cxx" => avrg_cxx,
"cxy" => avrg_cxy, "t" => t, "g_c" => g_c, "parameters" => par))

return t, avrg_cxx, avrg_cxx

end#main


################################
## Functions
################################

function gridit(t, T, vx, vy)
    # t is the time range (global)
    # T is the time vector of the specific particle
    VIX = interpolate( (T,), vx, Gridded(Linear()))
    VIY = interpolate( (T,), vy, Gridded(Linear()))
    VX = VIX[t]
    VY = VIY[t]
    return VX, VY
end

function velcor(t, T, VX, VY, rplan, iplan)
    # first create grided version
    vx, vy = gridit(t, T, VX, VY)
    nt = length(t)
    # zero-pad the velocity timeseries with as much zeroes
    append!(vx, zeros(Float64, nt))
    append!(vy, zeros(Float64, nt))
    # calculate fourier transforms using preplanned vectors
    ftvx = rplan*vx
    ftvy = rplan*vy
    ftvx_conj = conj(ftvx)
    # calculate correlation function based on product of fourier transforms
    cxx = (iplan*(ftvx_conj .* ftvx))[1:nt]/nt # devide by n again because of having 2 ffts
    cxy = (iplan*(ftvx_conj .* ftvy))[1:nt]/nt

    return cxx, cxy
end

function prepare_fourier(nt)
    # Define dummy that goes into rfft()
    dummy = zeros(Float64, 2*nt)
    normal_plan = plan_rfft(dummy)
    # Defne dummy that goes into irfft()
    dummy2 = complex(dummy)[1:(div(2*nt,2)+1)]
    inv_plan = plan_irfft(dummy2, 2*nt)
    return normal_plan, inv_plan
end

timer = time()
t, cxx, cxy = main(B)
timer = time() - timer

println("Process ended. Total time required was:")
println(round(timer, 3), " seconds, or ", round(timer/60, 3), " minutes")
