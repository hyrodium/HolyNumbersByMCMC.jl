# struct
struct Point
    x::Float64
    y::Float64
end

struct Affine
    xx::Float64
    yx::Float64
    xy::Float64
    yy::Float64
    Δx::Float64
    Δy::Float64
end

function norm(a::Point)
    return sqrt(a.x^2+a.y^2)
end

function Base.:+(a::Point,b::Point)
    return Point(a.x+b.x, a.y+b.y)
end

function Base.:-(a::Point,b::Point)
    return Point(a.x-b.x, a.y-b.y)
end

function Base.:/(a::Point,b::Real)
    return Point(a.x/b, a.y/b)
end

function Base.:*(k::Real,p::Point)
    return Point(k*p.x, k*p.x)
end

function Base.:*(A::Affine, a::Point)
    return Point(A.xx*a.x+A.xy*a.y+A.Δx, A.yx*a.x+A.yy*a.y+A.Δy)
end

function cci(a::Point,b::Point,r,s)
    L=norm(b-a)
    m=(r+s)/L
    n=(r-s)/L
    A=Affine(m*n, sqrt((m^2-1)*(1-n^2)), -sqrt((m^2-1)*(1-n^2)), m*n, 0, 0)
    return (a+b)/2 + A*(b-a)/2
end

const Locus = Vector{Point}

function area(locus::Locus)
    n = length(locus)
    s = locus[1].x*locus[2].y-locus[1].y*locus[2].x
    for i in 2:n
        s += locus[i-1].x*locus[i].y-locus[i-1].y*locus[i].x
    end
    return s/2
end

struct HolyNumbers
    w01::Float64
    w02::Float64
    w03::Float64
    w04::Float64
    w05::Float64
    w06::Float64
    w07::Float64
    w08::Float64
    w09::Float64
    w10::Float64
    w11::Float64
    w12::Float64
    w13::Float64
end

Base.isless(h1::HolyNumbers, h2::HolyNumbers) = energy(h1)<energy(h2)

struct Points
    q1::Point
    q2::Point
    q3::Point
    q4::Point
    q5::Point
    q6::Point
    q7::Point
    q8::Point
    function Points(holynumbers::HolyNumbers, θ::Real)
        q1 = Point(0.0,0.0)
        q2 = q1 + Point(w13*cos(θ),w13*sin(θ))
        q3 = Point(holynumbers.w01, holynumbers.w12)
        q4 = cci(q2,q3,holynumbers.w10,holynumbers.w02)
        q5 = cci(q3,q2,holynumbers.w03,holynumbers.w11)
        q6 = cci(q4,q3,holynumbers.w05,holynumbers.w04)
        q7 = cci(q6,q5,holynumbers.w06,holynumbers.w07)
        q8 = cci(q7,q5,holynumbers.w08,holynumbers.w09)
        return new(q1,q2,q3,q4,q5,q6,q7,q8)
    end
end

function Locus(holynumbers::HolyNumbers; n::Int=90)
    ptss = [Points(holynumbers, 2π*i/n) for i in 0:n-1]
    locus=[ptss[i].q8 for i in 1:n]
    return locus
end

## draw part
function lrdu(pts::Points)
    q1 = pts.q1
    q2 = pts.q2
    q3 = pts.q3
    q4 = pts.q4
    q5 = pts.q5
    q6 = pts.q6
    q7 = pts.q7
    q8 = pts.q8

    l = min(q1.x,q2.x,q3.x,q4.x,q5.x,q6.x,q7.x,q8.x)
    r = max(q1.x,q2.x,q3.x,q4.x,q5.x,q6.x,q7.x,q8.x)
    d = min(q1.y,q2.y,q3.y,q4.y,q5.y,q6.y,q7.y,q8.y)
    u = max(q1.y,q2.y,q3.y,q4.y,q5.y,q6.y,q7.y,q8.y)

    return [l,r,d,u]
end

function lrdu(locus::Locus)
    xs = [q.x for q in locus]
    ys = [q.y for q in locus]

    l = minimum(xs)
    r = maximum(xs)
    d = minimum(ys)
    u = maximum(ys)

    return [l,r,d,u]
end

function lrdu(ptss::AbstractVector{Points})
    tmp = hcat(lrdu.(ptss)...)

    l = minimum(tmp[1,:])
    r = maximum(tmp[2,:])
    d = minimum(tmp[3,:])
    u = maximum(tmp[4,:])

    return [l,r,d,u]
end

function lxrpt(q;unitlength = 5)
    Luxor.Point(unitlength*q.x,-unitlength*q.y)
end

function _linkagefig(pts::Points, locus::Locus; unitlength=5, showlength=false)
    q1 = pts.q1
    q2 = pts.q2
    q3 = pts.q3
    q4 = pts.q4
    q5 = pts.q5
    q6 = pts.q6
    q7 = pts.q7
    q8 = pts.q8
    q0 = Point(q1.x, q3.y)

    segments = [(q1,q2), (q2,q4), (q3,q4), (q2,q5), (q3,q5), (q4,q6), (q3,q6), (q6,q7), (q5,q7), (q5,q8), (q7,q8)]
    if showlength
        append!(segments, [(q0,q1), (q0,q3)])
    end

    Luxor.background(RGB(0.2,0.2,0.2))

    Luxor.sethue("red")
    for i in 1:length(locus)-1
        Luxor.line(lxrpt(locus[i],unitlength=unitlength),lxrpt(locus[i+1],unitlength=unitlength),:stroke)
    end
    Luxor.line(lxrpt(locus[end],unitlength=unitlength),lxrpt(locus[1],unitlength=unitlength),:stroke)

    Luxor.sethue("gray")
    for (q_s,q_e) in segments
        Luxor.line(lxrpt(q_s,unitlength=unitlength),lxrpt(q_e,unitlength=unitlength),:stroke)
    end

    Luxor.sethue("black")
    for q in [q1,q2,q3,q4,q5,q6,q7,q8]
        Luxor.circle(lxrpt(q,unitlength=unitlength),5,:fill)
    end

    if showlength
        Luxor.sethue("white")
        Luxor.fontsize(20)
        for (q_s,q_e) in segments
            Luxor.text(string(round(norm(q_e-q_s), digits=4)),  lxrpt((q_s+q_e)/2), halign=:center, valign=:center)
        end
    end
end

function draw(filename, pts::Points, locus::Locus; up=5, down=-5, right=5, left=-5, unitlength=5, showlength=false)
    Luxor.Drawing(unitlength*(right-left),unitlength*(up-down),filename)
    Luxor.origin(-unitlength*left,unitlength*up)

    _linkagefig(pts, locus, unitlength=unitlength, showlength=showlength)

    Luxor.finish()
    Luxor.preview()
end

function draw(filename, holynumbers::HolyNumbers, θ::Real; n::Int=90, unitlength=5, showlength=false)
    ptss = [Points(holynumbers, 2π*i/n) for i in 0:n-1]
    pts = Points(holynumbers, θ)
    locus = Locus(holynumbers, n=90)
    l,r,d,u = Int.(round.(lrdu(ptss)))
    up = u+5
    down = d-5
    right = r+5
    left = l-5

    Luxor.Drawing(unitlength*(right-left),unitlength*(up-down),filename)
    Luxor.origin(-unitlength*left,unitlength*up)

    _linkagefig(pts, locus, unitlength=unitlength, showlength=showlength)

    Luxor.finish()
    Luxor.preview()
end

function draw(filename, holynumbers::HolyNumbers; n::Int=90, unitlength=5, showlength=false)
    ptss = [Points(holynumbers, 2π*i/n) for i in 0:n-1]
    locus = Locus(holynumbers, n=n)
    l,r,d,u = Int.(round.(lrdu(ptss)))
    up = u+5
    down = d-5
    right = r+5
    left = l-5

    function frame(scene, framenumber)
        Luxor.origin(-unitlength*left,unitlength*up)
        _linkagefig(ptss[framenumber], locus, unitlength=unitlength, showlength=showlength)
    end

    w = unitlength*(right-left)
    h = unitlength*(up-down)
    movie = Luxor.Movie(w,h,"linkage_motion", 1:n)
    Luxor.animate(movie, [Luxor.Scene(movie, frame, 1:n)], creategif=true, pathname=filename)
end

# mcmc part

function randnext(holynumbers::HolyNumbers)
    w01 = holynumbers.w01 + randn()/4
    w02 = holynumbers.w02 + randn()/4
    w03 = holynumbers.w03 + randn()/4
    w04 = holynumbers.w04 + randn()/4
    w05 = holynumbers.w05 + randn()/4
    w06 = holynumbers.w06 + randn()/4
    w07 = holynumbers.w07 + randn()/4
    w08 = holynumbers.w08 + randn()/4
    w09 = holynumbers.w09 + randn()/4
    w10 = holynumbers.w10 + randn()/4
    w11 = holynumbers.w11 + randn()/4
    w12 = holynumbers.w12 + randn()/4
    w13 = 15.0
    return HolyNumbers(w01,w02,w03,w04,w05,w06,w07,w08,w09,w10,w11,w12,w13)
end

function mcmcnext(X, energy::Function; β=1.0)
    X′ = randnext(X)
    E = energy(X)
    E′ = energy(X′)
    exp(-β*E′)/exp(-β*E) < rand() && return X
    return X′
end

function iterate_mcmc(X_init, energy::Function, N::Int; β=1.0)
    X_tmp = X_init
    Xs = Array{typeof(X_init),1}(undef, N)
    for i in 1:N
        X_tmp = mcmcnext(X_tmp, energy, β=β)
        Xs[i] = X_tmp
    end
    return Xs
end

function bestsolution(holynumbers_log::Vector{HolyNumbers}, energy::Function)
    energys = energy.(holynumbers_log)
    i = argmin(energys)
    return holynumbers_log[i]
end
