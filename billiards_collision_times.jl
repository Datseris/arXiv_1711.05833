#=
Calculates with great statistics the histograms of the collision times,
as well as the first collision times in the Periodic Sinai billiard / RPSB
used in our paper https://arxiv.org/abs/1711.05833

Requires package StatsBase (any version) and FileIO that can load
both .jld and .jld2 files
=#

using DynamicalBilliards, StatsBase, FileIO
cd(@__DIR__)
include("Library.jl")

bt = billiard_rectangle(setting = "periodic")
r = 0.15
# The following can be Disk or RandomDisk:
# (change the savename if you change this!!!)
push!(bt, RandomDisk([0.5,0.5], r))

B = 0.99
savedir = "billiards/coltimes/"
# Use the knowledge of maximum collision time for a given magnetic field
# (only needs to be computed once)
maxt = load(savedir*"mct.jld2")["mt"][round(Int, B*100)]

# Time vector used as the partition of the histogram
ε = 0.01
tvector = 0.0:ε:maxt+ε
savename = "Hist_RPSB_dt=$(ε)_r=$(r)_B=$(B).jld"

println("Starting calculation, r=$r, B=$B")
tim = time()

function colt_hist(B, tvector)

    hist = zeros(tvector)
    # First collision times
    fcthist = Float64[]
    tt=1000.0
    parnum = 100000
    pinned = 0
    pnum = 0; perc = 0

    while pnum < parnum
        perc = percentage(pnum, parnum, perc)
        #particle loop
        t, poss, vels = evolve!(randominside(bt, 2B), bt, tt)
        if t[end] == Inf
            #println("Pinned")
            pinned += 1
            continue
        else
            pnum +=1
        end

        append!(fcthist, t[2]) # t[1] is always zero by definition

        H = fit(Histogram, t[3:end], tvector, closed=:left)
        we =  H.weights/length(t[3:end])
        push!(we, zero(eltype(we)))
        hist .+= we

    end#particle loopS
    hist ./= parnum

    gc = 1 - (pinned)/(parnum+pinned)

    H = fit(Histogram, fcthist, tvector, closed = :left)
    we2 = H.weights/length(fcthist)
    push!(we2, zero(eltype(we2)))

    return tvector, hist, we2
end

ed, we, we2 = colt_hist(B, tvector)

tim = time() - tim
println("finished. Took a total of:")
println("$tim seconds, or $(round(tim/60, 3)) minutes")
println("\nWriting data in file now... (this should take at most 1 second)");
flush(STDOUT)
save(savedir*savename), "t", ed, "colt", we/ε, "fct", we2/ε)
println("Done!")
