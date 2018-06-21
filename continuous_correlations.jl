#=
This file creates the velocity correlation functions of the
continuous antidot system used
in our paper https://arxiv.org/abs/1711.05833
=#
using FileIO, OrdinaryDiffEq
cd(@__DIR__)
include("Library.jl")

B = 0.32              # use this when running interactively
# B = parse(ARGS[1])    # use this when running from terminal

######################################
# Configure system                   #
######################################
function configure_system(B)
    # Potential type
    pt = :std # :fkg, :br, :std

    # Fundamental parameters:
    d0 = 0.3   #antidot diameter
    c = 0.1  #cutoff distance, must be 0 < c < 0.5-d0/2 !!!

    # FKG (old)
    β = 4

    # Boundary Roughness parameters (for ε = 0 no BR takes place)
    ε = 0.0  #relative edge roughness
    M = 16   #number of sine modes. complexity of edge roughness
    matrix_size = 100  #size of the matrix used for the random coefficients
    matrixnumber = 1

    dt = 0.01; T = 5000.0
    t = 0:dt:T
    N = 5000

    par =
    Dict(:pt => pt, :d0 => d0, :c=>c,
    :β=>β, :ε =>ε, :dt => dt, :T=>T, :N=>N)

    prob, potentialenergy = setup_ODEProblem(par)
    return prob, potentialenergy, par
end


######################################
#                                    #
##------      MAIN: START     ------##
#                                    #
######################################
function main(B)

    prob, potentialenergy, par = configure_system(B)

    pt, d0, c, β, ε, dt, T, N = unzip(par)
    output = jldname(par, "Cor";  B=B)
    spath = "data/correlations/"
    ispath(spath) || mkpath(spath)

    if isfile(spath*output)
        println("File already exists!!!")
        return zeros(5), zeros(5), zeros(5)
    end

    if B>=0.17 # Value for n=21 peak. Don't care about higher orders
        cyclotron = 1/(2*B)    #cyclotron radius at given B
        chaotic_t = (π/B)*10  #time to check whether an orbit is chaotic or not.
        #(π/Β) is time for one rotation
    else
        cyclotron = 0.0  #cyclotron diameter for chaotic check
        chaotic_t = 5.0
    end


    start_time = time()
    println("\n----------------------------------")
    println("Trajectory Integrator and Correlation Calculator.")
    println("Integration time: T = $(par[:T]), total trajectories: $(par[:N]).")
    println("Output name: ")
    println(output); flush(STDOUT)

    approx_tick = dt;
    ### PREPARE FOURIER TRANSFORMS ###
    approx_nt = round(Int64, (T/approx_tick)) #approximate number of points
    nt = round(Int64,(2^(round(log(approx_nt)/log(2)))))   #number of points near pow(2)
    tick = T/nt
    t = collect(linspace(0.0, T, nt))
    id_chaotic = closest_point(t, chaotic_t)

    # Set up correlations:
    print("\nPlanning FFT. . .")
    avrg_cxx = zeros(Float64, nt); avrg_cxy = zeros(Float64, nt);
    cxx = zeros(Float64, nt); cxy = zeros(Float64, nt);
    vcf! = VCF(nt)
    print("  Done.\n"); flush(STDOUT)

    orbit_count = 0
    pinned_count = 0

    integ = init(prob, Tsit5(), saveat = t, abstol = 1e-9, reltol = 1e-9, dtmax = 0.05,
    maxiters = Inf)

    ###### Main Orbit loop
    println("!!!Starting main orbit loop!!!")
    while orbit_count < N
        println("$(orbit_count+1) / $N"); flush(STDOUT)

        # get proper initial conditions, with energy = 1:
        var_0 = initial(:ultra, potentialenergy)

        # Integrate up to chaotic check:
        reinit!(integ, SVector{4}(var_0))
        solve!(integ); sol = integ.sol

        # Calculate maximum distance from origin (within chaotic check):
        if B < 0.17
            max_dist = Inf
        else

            x_1 = sol[1, 1:id_chaotic]
            y_1 = sol[2, 1:id_chaotic]
            max_dist = sqrt(maximum( @. (x_1 - x_1[1])^2 + (y_1 - y_1[1])^2))
        end

        # Found orbit is pinned:
        if max_dist <= 2.0*cyclotron + 0.5
            pinned_count += 1
            continue
        # Found orbit is chaotic:
        else
            orbit_count += 1
        end

        # Calculate correlations
        vx = sol[3, :]
        vy = sol[4, :]
        vcf!(cxx, cxy, vx, vy)
        avrg_cxx .+= cxx; avrg_cxy .+= cxy

    end#end of N loop

    avrg_cxx ./= N; avrg_cxy ./= N;

    g_c = 1.0 - (pinned_count)/(N + pinned_count)

    # Create JLD (julia data) file with correlations (binary):
    save(spath*output, Dict("cxx" => avrg_cxx, "cxy" => avrg_cxy, "t" => t, "g_c" => g_c, "parameters" => par))
    end_time = time() - start_time
    println("\nEnded.")
    println("Process took $(round(end_time, 2)) seconds, or $(round(end_time/60, 3)) minutes
    equivalent to $(round(end_time/(T*N), 5)) seconds per unit of time."); flush(STDOUT)
    return avrg_cxx, avrg_cxy, t
end

avrg_cxx, avrg_cxy, t = main(B)
