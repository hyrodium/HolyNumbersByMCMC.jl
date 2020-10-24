import Luxor
using Colors
using BenchmarkTools

include("_holynumbersbymcmc.jl")

##  Initial holynumbers
w01 = 40
w02 = 40
w03 = 40
w04 = 40
w05 = 60
w06 = 40
w07 = 40
w08 = 60
w09 = 50
w10 = 50
w11 = 60
w12 = 5
w13 = 15

holynumbers = HolyNumbers(w01,w02,w03,w04,w05,w06,w07,w08,w09,w10,w11,w12,w13)
draw("initial.gif", holynumbers)

## Energy functions

function energy(holynumbers::HolyNumbers; n=90)
    try
        return energy(Locus(holynumbers, n=n))
    catch
        return Inf
    end
end

function energy(locus::Locus)
    E = 0.0
    N = length(locus)
    L,R,D,U = -10,200,-100,-50
    l,r,d,u = lrdu(locus)
    l < L && return Inf
    r > R && return Inf
    d < D && return Inf
    u > U && return Inf

    # Gravity near ground
    s = 0.0
    for q in locus
        s += (q.y - D)^(2/3)
    end
    E += s/N

    # area
    E += (area(locus)-600)^2/10000

    # D - position
    # E +=(D+90)^2

    return E
end

energy(holynumbers)

holynumbers_log = iterate_mcmc(holynumbers, energy, 10000)
minimum(energy.(holynumbers_log))
best_holynumbers = bestsolution(holynumbers_log, energy)
holynumbers_log = iterate_mcmc(best_holynumbers, energy, 10000, β=4.0)
minimum(energy.(holynumbers_log))
best_holynumbers = bestsolution(holynumbers_log, energy)
holynumbers_log = iterate_mcmc(best_holynumbers, energy, 10000, β=16.0)
minimum(energy.(holynumbers_log))
best_holynumbers = bestsolution(holynumbers_log, energy)

draw("output/tmp_021.gif", best_holynumbers)
