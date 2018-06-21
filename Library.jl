using FileIO, StaticArrays

"""
Applies the relativistic equations of motion from the potential derivatives.
"""
function apply_eom(u, U, Ux, Uy, B)
    # Apply hyper-relativistic equations of motion:
    der = (u[3]*Uy - u[4]*Ux + 2*B)/(1-U)
    return SVector{4, Float64}(u[3], u[4], +u[4]*der, -u[3]*der)
end

"""
    antidot_eom(du, u, p, t)
Equations of motion for antidot model without boundary roughness.

The parameter container must be `p = (B, d0, c)`.
"""
function antidot_eom(u, p, t)
    B, d0, c = p
    x, y, vx, vy = u
    # Calculate quadrant of (x,y):
    U = Uy = Ux = 0.0
    xtilde = x - round(x);  ytilde = y - round(y)
    ρ = sqrt(xtilde^2 + ytilde^2)
    # Calculate derivatives and potential:
    if ρ < 0.5*d0 + c
        U = (1/(c^4))*(0.5d0 + c - ρ)^4 #potential
        sharedfactor = -(4*(c + d0/2 - ρ)^(3))/(c^4*ρ)
        Ux = sharedfactor*xtilde # derivatives
        Uy = sharedfactor*ytilde
    end
    return apply_eom(u, U, Ux, Uy, B)
end

function antidot_potential(x::Float64, y::Float64, p)
    B, d0, c = p
    δ = 4
    # Calculate quadrant of (x,y):
    xtilde = x - round(x);  ytilde = y - round(y)
    ρ = sqrt(xtilde^2 + ytilde^2)
    # Check distance:
    pot = ρ > 0.5*d0 + c ? 0.0 : (1/(c^δ))*(0.5*d0 + c - ρ)^δ
end

"""
    br_antidot_eom(du, u, p, t)
Equations of motion for antidot model with boundary roughness.

The parameter container must be `p = (B, d0, c, ε, R)` with R the
matrix of the random numbers (3D).
"""
function br_antidot_eom(u, p, t)
    B, d0, c, ε, R = p
    δ = 4
    x, y, vx, vy = u
    # Calculate quadrant of (x,y):
    U = Uy = Ux = 0.0
    Rsize, _dummy, M = size(R)
    xtilde = x - round(x);  ytilde = y - round(y)
    idx = mod(round(Int64, x), Rsize) + 1
    idy = mod(round(Int64, y), Rsize) + 1
    ρ = sqrt(xtilde^2 + ytilde^2)
    # Calculate effective distorted antidot diameter:
    φ = atan2(ytilde, xtilde)
    sinsum = 0.0
    cossum = 0.0
    for k in 1:M
        sinsum += sin(k*φ)*R[idx, idy, k]
        cossum += cos(k*φ)*k*R[idx, idy, k]
    end
    # current diameter:
    d = d0 + d0*ε*sinsum
    # Calculate derivatives and potential:
    if ρ < 0.5*d + c
        U = (1/(c^δ))*(0.5*d + c - ρ)^δ
        dddφ = d0*ε*cossum
        sharedfactor = -δ/((c^δ)*(ρ^2))*(0.5*d + c - ρ)^(δ-1)
        Ux = sharedfactor*(xtilde*ρ + 0.5*dddφ*ytilde)
        Uy = sharedfactor*(ytilde*ρ - 0.5*dddφ*xtilde)
    end
    return apply_eom(u, U, Ux, Uy, B)
end

function br_antidot_potential(x::Float64, y::Float64, p)
    B, d0, c, ε, R = p
    δ = 4
    # Calculate quadrant of (x,y):
    Rsize, _dummy, M = size(R)
    xtilde = x - round(x);  ytilde = y - round(y)
    idx = mod(round(Int64, x), Rsize) + 1
    idy = mod(round(Int64, y), Rsize) + 1
    ρ = sqrt(xtilde^2 + ytilde^2)
    # Calculate effective distorted antidot diameter:
    φ = atan2(ytilde, xtilde)
    sinsum = 0.0
    for k in 1:M
        sinsum += sin(k*φ)*R[idx, idy, k]
    end
    d = d0 + d0*ε*sinsum
    pot = ρ > 0.5*d + c ? 0.0 : (1/(c^δ))*(0.5*d + c - ρ)^δ
end

"""
    antidot_eom(du, u, p, t)
Equations of motion for antidot model used by Fleischmann (PRL 1992)

The parameter container must be `p = (B, d0, β)`.
"""
function fkg_antidot_eom(u, p, t)
    Β, d, β = p
    x, y = u
    pref = 1/(cos(d*π/2)^β)
    cosx = cos(π*x); cosy = cos(π*y)
    U = pref*(cosx*cosy)^β
    Ux = -π*β*U*tan(π*x)
    Uy = -π*β*U*tan(π*y)
    return apply_eom(u, U, Ux, Uy, B)
end

function fkg_antidot_potential(x::Float64, y::Float64, p)
    Β, d, β = p
    pref = 1/(cos(d*π/2)^β)
    cosx = cos(π*x); cosy = cos(π*y)
    U = pref*(cosx*cosy)^β
end

function unzip(p)
    pt = p[:pt]; d0 = p[:d0]; c = p[:c]; β = p[:β]; ε = p[:ε]
    dt = p[:dt]; T = p[:T]; N = p[:N]
    return pt, d0, c, β, ε, dt, T, N
end

function setup_ODEProblem(par)
    # Set up ODEProblem
    pt, d0, c, β, ε, dt, T, N = unzip(par)
    u0 = SVector{4}(rand(4)...)
    if pt == :std
        p0 = (B, d0, c)
        potentialenergy = (x, y) -> antidot_potential(x, y, p0)
        prob = ODEProblem{false}(antidot_eom, u0, (0.0, T), p0)
    elseif pt == :fkg
        p0 = (B, d0, β)
        potentialenergy = (x, y) -> fkg_antidot_potential(x, y, p0)
        prob = ODEProblem(fkg_antidot_eom, u0, (0.0, T), p0)
    elseif pt == :br
        # Change `matrixname` to use other random matrices
        matrixname = "Uniform_M=16_size=100_1.jld2"
        R = load("data/randommatrix/"*matrixname)["matrix"]
        p0 = (B, d0, c, ε, R)
        potentialenergy = (x, y) -> br_antidot_potential(x, y, p0)
        prob = ODEProblem(br_antidot_eom, u0, (0.0, T), p0)
    end
    return prob, potentialenergy
end


"""
    initial(system_type::Symbol, potentialenergy::Function)
Return a 4-element initial condition vector for the relativistic system that
has random initial conditions (for the momentum case).
Satisfies: `kineticenergy(vector) + potentialenergy(vector) = 1.0` .
# Arguments:
* system_type : `:relativistic` or `:ultra` (use `:ultra` for paper).
* potentialenergy : function that returns potential energy given `(x, y)` only!
"""
function initial(system_type::Symbol, potentialenergy::Function)
    var0 = zeros(Float64, 4)
    var0[1] = rand()
    var0[2] = rand()
    # Check if the initial condition is valid:
    pot0 = potentialenergy(var0[1], var0[2])
    while pot0 >= 0.99
        var0[1] = rand()
        var0[2] = rand()
        pot0 = potentialenergy(var0[1], var0[2])
    end

    θ0 = -π + 2*π*rand()

    if system_type == :relativistic
        var0[3] = (1-pot0)*cos(θ0) - var0[2]*B
        var0[4] = (1-pot0)*sin(θ0) + var0[1]*B
    elseif system_type == :ultra
        var0[3] = cos(θ0)
        var0[4] = sin(θ0)
    else
        error("Unspecified system.")
    end
    return var0
end

"""
    plot_orbit(x::Array{Float64}, y::Array{Float64}, potentialenergy;
    color = (0.0, 0.0, 1.0))
Plot given orbit `(x, y)` on top of the potential landscape,
defined by the function `potentialenergy`.
Plot's title will be the initial conditions `var_0_plot = (x(0), y(0), θ(0))`.

Optionally specify color of the orbit's curve in (r,g,b) format.

Optionally add a textbox containing the values of the parameters used.
Warning: the textbox (currently) needs all the parameters to be in global scope.
"""
function plot_orbit(x::Array{Float64}, y::Array{Float64}, potentialenergy;
    color = (0.0, 0.0, 1.0), textbox = true, n=4*512)

    x0 = x[1]; y0 = y[1]; θ0 = atan2(y0, x0)
    limits = [minimum(x), maximum(x), minimum(y), maximum(y)]


    plot_potential(potentialenergy, limits)

    # Plot orbit:
    PyPlot.plot(x,y, c = color, linewidth = 1.0)
    PyPlot.plot(x[1],y[1], c = color, marker = "o")
    PyPlot.plot(x[end],y[end], c = "black", marker = "o")
    # Add title with initial conditions:
    PyPlot.title("i.c.: \$ x_0 = $(round(x0,3)), y_0 = $(round(y0,3))" *
               ", \\theta_0 = $(round(θ0,3))\$",
               size = 20, family="serif")
    tight_layout()
end

function plot_potential(potentialenergy, limits=[-1,2,-1,2], n::Int64 = 4*512;
    cb = false)

    #point density for potential plot
    xlin = linspace(limits[1], limits[2], n)
    ylin = linspace(limits[3], limits[4], n)
    # define arrays necesarry for plotting
    xgrid = repmat(xlin', n, 1)
    ygrid = repmat(ylin ,1,n)
    pot = zeros(Float64, n, n)

    # calculate potential for plotting
    for col in 1:n, row in 1:n
        pot[row, col] = potentialenergy(xlin[col], ylin[row])
    end

    # Create figure and axes objects
    fig, ax = PyPlot.subplots()
    # make contour plot with color-filled areas:
    cpfilled = ax[:contourf](xgrid, ygrid, pot, 8, alpha=1.0, vmin = 0.1,
    extend="max", cmap="Greys" #=, levels=[1.0, 2.0],  =#)

    # add red contour line at potential = 1.0
    c1line = ax[:contour](xgrid, ygrid, pot, colors="red", linewidth=1.0,
    levels = [1.0])
                        # make contour plot with lines:
    cplines = ax[:contour](xgrid, ygrid, pot, 8, colors="black", linewidth=0.5,
    linestyle = "dashed", levels=[0.0])

    for c in cplines[:collections]
        c[:set_linestyle]("dashed")
    end

    gca()[:set_aspect]("equal")
    cb && (cbar = colorbar(cpfilled))
    cb && cbar[:add_lines](c1line)
    cb && cbar[:ax][:set_ylabel]("\$\\mathcal{U}\$")


    # add labels, adjust fonts and sizes:
    PyPlot.xlabel("\$x\$", size=28)
    PyPlot.ylabel("\$y\$", size=28)
    PyPlot.xticks(family="serif", size = 18)
    PyPlot.yticks(family="serif", size = 18)
end




"""
    generate_res(params, tau, keep_pinned, plot_res = true)
Load correlations (or load resistivities directly
if they exist) and then generate resistivities versus magnetic field. If
the resistivities already exist, it justs loads them. Otherwise, it also saves
them.

Both `tau` and `keep_pinned` can be arrays with possible values.

`keep_pinned` can only take three values: `true, false, "nogc"` which
correspond to using entire phasespace, using only chaotic orbits or using only
chaotic orbits without involving at all the g_c factor (like setting g_c to 1
for all magnetic fields)

Returns the lastly computed resistivities as:
    B, Rxx, Rxy, g_c
"""
function generate_res(parameters, τ, keep_pinned;
    cor_load_path = "data/correlations/",
    res_save_path = "data/resistivities/")

    Blist = Float64[]; Rxxlist = Float64[]
    Rxylist = Float64[]; gclist = Float64[]
    ### MAIN ###

    # Closure that loads correlations
    loadcors = (B) -> begin
        cname = jldname(parameters, "Cor"; B=B)
        loadpath = joinpath(cor_load_path, cname)
        f = load(loadpath)
        t = f["t"]; cxx = f["cxx"]; cxy=f["cxy"]; g_c = f["g_c"]

        return t, cxx, cxy, g_c
    end

    for tau in τ
        for kp in keep_pinned
            println("τ = $tau")
            # Save names:
            resname = jldname(parameters, "Res"; τ=tau, kp = kp)
            savepath = joinpath(res_save_path, resname)
            # Resistivities:
            if isfile(savepath)
                ff = load(savepath)
                Blist=ff["B"]; Rxxlist=ff["Rxx"]; Rxylist=ff["Rxy"]
                gclist = ff["g_c"]
            else
                Blist, Rxxlist, Rxylist, gclist = resistivities(tau, kp, loadcors)
                save(savepath, "B", Blist, "Rxx", Rxxlist, "Rxy", Rxylist, "g_c", gclist)
            end
        end
    end
    return Blist, Rxxlist, Rxylist, gclist
end

"""
    resistivities(τ, keep_pinned, loadcors)
Calculate resistivity curves with impurity scattering τ and with
handling of the g_c factor using `keep_pinned`.
Uses the closure `loadcors` to load correlations.

`keep_pinned` can only take three values: `true, false, "nogc"` which
correspond to using entire phasespace, using only chaotic orbits or using only
chaotic orbits without involving at all the g_c factor (like setting g_c to 1
for all magnetic fields).
"""
function resistivities(τ, keep_pinned, loadcors, Bs = 0.01:0.01:1.2)

    Blist = zeros(Bs)
    Rxxlist = zeros(Bs)
    Rxylist = zeros(Bs)
    gclist = zeros(Bs)

    for (i, B) in enumerate(Bs)
        t, cxx, cxy, g_c = loadcors(B)

        if keep_pinned == true
            Cxx = @. cxx*g_c + 0.5*(1-g_c)*cos(2B*t)
            Cxy = @. cxy*g_c + 0.5*(1-g_c)*sin(2B*t)
        elseif keep_pinned == false
            Cxx = @. cxx*g_c
            Cxy = @. cxy*g_c
        elseif keep_pinned == "nogc"
            Cxx = cxx
            Cxy = cxy
        else
            error("Incorrect value for `keep_pinned`")
        end

        sxx = trapezint(t[2]-t[1], Cxx .* exp.(-t ./ τ))
        sxy = trapezint(t[2]-t[1], Cxy .* exp.(-t ./ τ))
        ssq = sxx^2 + sxy^2
        Rxx = sxx/ssq; Rxy = sxy/ssq

        Rxxlist[i] = Rxx
        Rxylist[i] = Rxy
        gclist[i] = g_c
    end
    return Bs, Rxxlist, Rxylist, gclist
end

"""
    jldname(p, prefix, ending = ".jld2"; kwargs...)
"""
function jldname(p, prefix, ending = ".jld2"; kwargs...)
    kw = Dict(kwargs)

    # Core text:
    filename = "_pt=$(p[:pt])_T=$(p[:T])_N=$(p[:N])"
    # Potential/System parameters
    if p[:pt] == :std
        filename *= "_d0=$(p[:d0])_c=$(p[:c])_e=$(p[:ε])"
    elseif p[:pt] == :std
        filename *= "_d0=$(p[:d0])_c=$(p[:c])"
    elseif p[:pt] == :fkg
        filename *= "_d0=$(p[:d0])_b=$(p[:β])"
    elseif p[:pt] ∈ [:PSB, :RPSB]
        filename *= "_d0=$(p[:d0])"
    end

    # Adding prefix and extra variables:
    s = ""
    for k in keys(kw)
        s*="_"*string(k)*"=$(kw[k])"
    end

    return prefix*filename*s*ending
end


"""
    trapezint(x::AbstractVector, y::AbstractVector)
    trapezint(step::Real, y::AbstractVector)
Calculate the integral of `y` versus `x` using the trapezoidal rule.
The second method assumes equally-spaced `y`.
"""
function trapezint(x::AbstractVector{X}, y::AbstractVector{Y}) where {X<:Real, Y<:Real}
  integral::Float64 = 0.0
  @inbounds for i in 1:(length(x)-1)
    integral += (x[i+1] - x[i])*(y[i+1] + y[i])
  end
  integral *= 0.5
end
function trapezint(step::Real, y::AbstractVector{Y}) where {Y<:Real}
  integral::Float64 = 0.0
  @inbounds for i in 1:(length(y)-1)
    integral += step*(y[i+1] + y[i])
  end
  integral *= 0.5
end

# Velocity correlation functions
struct VCF{S, P}
    vxpadded::Vector{Float64}
    vypadded::Vector{Float64}
    normal_plan::S
    inv_plan::P
    nt::Int
end

function VCF(nt::Int, k::Int = 1)
    dummy = zeros(Float64, (k+1)*nt)
    normal_plan = plan_rfft(dummy)
    ivx = normal_plan*dummy
    inv_plan = plan_irfft(ivx, (k+1)*nt)
    VCF(dummy, deepcopy(dummy), normal_plan, inv_plan, nt)
end

function (vcf::VCF)(cx, cy, vx, vy)
    nt = vcf.nt
    vcf.vxpadded[1:nt] .= vx
    vcf.vypadded[1:nt] .= vy

    ftvx = vcf.normal_plan*vcf.vxpadded
    ftvy = vcf.normal_plan*vcf.vypadded
    ftvx_conj = conj(ftvx)

    # calculate correlation function based on product of fourier transforms
    cx .= +(vcf.inv_plan*(ftvx_conj .* ftvx))[1:nt]/nt # devide by n again because of having 2 ffts
    cy .= -(vcf.inv_plan*(ftvx_conj .* ftvy))[1:nt]/nt  #negative sign coz of XY instead of YX
    return
end

function billiard_table(p)
    if p[:pt] == :RPSB
        bt = billiard_rectangle(setting = "periodic")
        a = RandomDisk([0.5, 0.5], p[:d0]/2)
        push!(bt, a)
    elseif p[:pt] == :PSB
        bt = billiard_sinai(p[:d0]/2, setting = "periodic")
    end
    return bt
end

function closest_point(grid, x)
    best = 1
    dxbest = abs(x - grid[1])
    for i in 2:length(grid)
        dx = abs(x - grid[i])
        if dx < dxbest
            dxbest = dx
            best = i
        end
    end
    return best
end
