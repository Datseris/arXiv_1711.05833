#=
Calculates the mean collision time in the Periodic Sinai billiard / RPSB
used in our paper https://arxiv.org/abs/1711.05833

Requires package StatsBase (any version) and FileIO that can load
both .jld and .jld2 files
=#

using DynamicalBilliards, FileIO

cd(@__DIR__)
cd("data/coltimes")

const r = 0.15

function mean_collision_time(B)

    PT = B == 0 ? Particle : MagneticParticle
    const bt = billiard_sinai(r, setting = "periodic")
    partnum = 20000
    n = 20000 # collisions
    mct = 0.0
    pinned = 0;  pnum = 0

    while pnum < partnum
        p::PT = B == 0 ? randominside(bt) : randominside(bt, 2B)

        t, poss, vels = evolve!(p, bt, n)
        if t[end] == Inf
            pinned += 1
            continue
        else
            pnum +=1
        end
        # calculate mean collision time of particle
        mct += sum(t)/n
    end
    gc = 1 - (pinned)/(partnum+pinned)
    return mct/partnum, gc
end

B = length(ARGS) == 0 ? 0.0 : parse(ARGS[1])

mct, g_c = mean_collision_time(B)
writedlm("mean_colt_r=$(r)._B=$(B).tsv", [mct g_c])
