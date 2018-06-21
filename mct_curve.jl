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
    const bt = billiard_sinai(r, setting = "periodic")
    partnum = 20000
    n = 20000 # collisions
    mct = 0.0
    pinned = 0;  pnum = 0

    while pnum < partnum
        p = randominside(bt, 2B)

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


function mean_collision_time()
    const bt = billiard_sinai(r, setting = "periodic")
    partnum = 20000
    n = 20000 # collisions
    mct = 0.0
    pinned = 0;  pnum = 0

    while pnum < partnum
        println(pnum)
        p = randominside(bt)

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


# function main()
#     Bs = 0:0.01:1.2
#
#     colts = Float64[3.100873430510771]
#     gcs = Float64[1.0]
#
#     for B in Bs[2:end]
#         println("B = $B")
#         mct, gc = mean_collision_time(B)
#         push!(colts, mct)
#         push!(gcs, gc)
#
#         save("mean_colt_r=$(r).jld2", Dict("mct" => colts, "gc" => gcs, "B" => Bs))
#     end
#     return Bs, colts, gcs
# end

# B = parse(ARGS[1])
B = 0.0

mct, g_c = mean_collision_time()
writedlm("mean_colt_r=$(r)._B=$(B).tsv", [mct g_c])
