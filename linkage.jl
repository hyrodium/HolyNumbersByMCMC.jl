import Luxor
using Colors
using BenchmarkTools

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
    function Points(t,holynumbers::HolyNumbers)
        q1 = Point(0.0,0.0)
        q2 = q1 + Point(w13*cos(t),w13*sin(t))
        q3 = Point(holynumbers.w01,-holynumbers.w12)
        q4 = cci(q2,q3,holynumbers.w10,holynumbers.w02)
        q5 = cci(q3,q2,holynumbers.w03,holynumbers.w11)
        q6 = cci(q4,q3,holynumbers.w05,holynumbers.w04)
        q7 = cci(q6,q5,holynumbers.w06,holynumbers.w07)
        q8 = cci(q7,q5,holynumbers.w08,holynumbers.w09)
        return new(q1,q2,q3,q4,q5,q6,q7,q8)
    end
end

function Locus(holynumbers::HolyNumbers; n::Int=90)
    ptss = [Points(2π*i/n,holynumbers) for i in 0:n-1]
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

function draw(filename, pts::Points, locus::Locus; up=5, down=-5, right=5, left=-5, zoom=1, unitlength=5)
    q1 = pts.q1
    q2 = pts.q2
    q3 = pts.q3
    q4 = pts.q4
    q5 = pts.q5
    q6 = pts.q6
    q7 = pts.q7
    q8 = pts.q8

    # Luxor.Drawing(2000, 2000, filename)
    Luxor.Drawing(unitlength*(right-left),unitlength*(up-down),filename)
    Luxor.origin(-unitlength*left,unitlength*up)
    # Luxor.origin()
    Luxor.background(RGB(0.2,0.2,0.2))

    Luxor.sethue("red")
    for i in 1:length(locus)-1
        Luxor.line(lxrpt(locus[i],unitlength=unitlength),lxrpt(locus[i+1],unitlength=unitlength),:stroke)
    end
    Luxor.line(lxrpt(locus[n],unitlength=unitlength),lxrpt(locus[1],unitlength=unitlength),:stroke)

    Luxor.sethue("gray")
    for (q_s,q_e) in [(q1,q2), (q2,q4), (q3,q4), (q2,q5), (q3,q5), (q4,q6), (q3,q6), (q6,q7), (q5,q7), (q5,q8), (q7,q8)]
        Luxor.line(lxrpt(q_s,unitlength=unitlength),lxrpt(q_e,unitlength=unitlength),:stroke)
    end

    Luxor.sethue("black")
    for q in [q1,q2,q3,q4,q5,q6,q7,q8]
        Luxor.circle(lxrpt(q,unitlength=unitlength),5,:fill)
    end

    Luxor.finish()
    Luxor.preview()
end


##  Original holynumbers
w01 = 38.0
w02 = 41.5
w03 = 39.3
w04 = 40.1
w05 = 55.8
w06 = 39.4
w07 = 36.7
w08 = 65.7
w09 = 49.0
w10 = 50.0
w11 = 61.9
w12 = 7.8
w13 = 15.0

holynumbers = HolyNumbers(w01,w02,w03,w04,w05,w06,w07,w08,w09,w10,w11,w12,w13)

n=90
ptss = [Points(2π*i/n,holynumbers) for i in 0:n-1]
locus=Locus(holynumbers, n=90)

l,r,d,u = Int.(round.(lrdu(ptss)))
for i in 1:n
    draw("original_$(lpad(i,3,"0")).png", ptss[i], locus, up=u+5, down=d-5, right=r+5, left=l-5, unitlength=3)
end

## Energy functions
function miny(locus::Locus)
    return minimum(p.y for p in locus)
end

function maxy(locus::Locus)
    return maximum(p.y for p in locus)
end

function energy(holynumbers::HolyNumbers; n=90)
    energy(Locus(holynumbers, n=n))
end

## mcmc part
function randnext(holynumbers::HolyNumbers)
    w01 = holynumbers.w01 + randn()
    w02 = holynumbers.w02 + randn()
    w03 = holynumbers.w03 + randn()
    w04 = holynumbers.w04 + randn()
    w05 = holynumbers.w05 + randn()
    w06 = holynumbers.w06 + randn()
    w07 = holynumbers.w07 + randn()
    w08 = holynumbers.w08 + randn()
    w09 = holynumbers.w09 + randn()
    w10 = holynumbers.w10 + randn()
    w11 = holynumbers.w11 + randn()
    w12 = holynumbers.w12 + randn()
    w13 = 15.0
    return HolyNumbers(w01,w02,w03,w04,w05,w06,w07,w08,w09,w10,w11,w12,w13)
end

function mcmcnext(X;β=1.0)
    X′ = randnext(X)
    E = energy(X)
    try
        E′ = energy(X′)
        exp(-β*E′)/exp(-β*E) < rand() && return X
    catch
        E′ = Inf
        exp(-β*E′)/exp(-β*E) < rand() && return X
    end
    return X′
end

function iterate_mcmc(X_init, N; β=1.0)
    X_tmp = X_init
    Xs = Array{typeof(X_init),1}(undef, N)
    for i in 1:N
        X_tmp = mcmcnext(X_tmp, β=β)
        Xs[i] = X_tmp
    end
    return Xs
end

## try
function energy(locus::Locus)
    E = 0.0
    Y = miny(locus)
    N = length(locus)
    # A: thinner
    # return maxy(locus)-miny(locus)

    # B: check locus height
    # return (miny(locus)+92)^2

    # C: electric field
    s = 0.0
    for q in locus
        s += (q.y - Y)^(1/2)
    end
    E += s/N

    # D: electric field
    # s = 0.0
    # Y = miny(locus)
    # for q in locus
    #     s += (q.y - Y)^(1/2)
    # end
    # E += s/N

    L,R,D,U = 0,200,-200,-50
    l,r,d,u = lrdu(locus)
    l < L && return Inf
    r > R && return Inf
    d < D && return Inf
    u > U && return Inf

    # E: gravity near ground
    # s = 0.0
    # Y = miny(locus)
    # for q in locus
    #     s += (tanh(q.y - Y))^2
    # end
    # E += s/N

    # area
    E += (area(locus)+250)^2/10000

    return E
end

area(locus)

area(Locus(holynumbers))

energy(holynumbers)
holynumbers_log = iterate_mcmc(holynumbers,10000)
energy(minimum(holynumbers_log))
@time holynumbers_log = iterate_mcmc(minimum(holynumbers_log),10000, β=2.0)
energy(minimum(holynumbers_log))

## hoge
holynumbers′ = minimum(holynumbers_log)
# holynumbers′ = holynumbers_log[720]
energy(holynumbers′)

n=90
ptss = [Points(2π*i/n,holynumbers′) for i in 0:n-1]
locus=Locus(holynumbers′)

l,r,d,u = Int.(round.(lrdu(ptss)))

for i in 1:n
    draw("linkage_$(lpad(i,3,"0")).png", ptss[i], locus, up=u+5, down=d-5, right=r+5, left=l-5, unitlength=3)
end
