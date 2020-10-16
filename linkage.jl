import Luxor
using LinearAlgebra
using Colors
using BenchmarkTools

function cci(a,b,r,s)
    L=norm(b-a)
    m=(r+s)/L
    n=(r-s)/L
    M=[m*n -sqrt((m^2-1)*(1-n^2));sqrt((m^2-1)*(1-n^2)) m*n]
    return (a+b)/2 + M*(b-a)/2
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

struct PointPositions
    q1::Array{Float64,1}
    q2::Array{Float64,1}
    q3::Array{Float64,1}
    q4::Array{Float64,1}
    q5::Array{Float64,1}
    q6::Array{Float64,1}
    q7::Array{Float64,1}
    q8::Array{Float64,1}
    function PointPositions(t,holynumbers::HolyNumbers)
        q1 = [0.0,0.0]
        q2 = q1 + w13*[cos(t),sin(t)]
        q3 = [holynumbers.w01,-holynumbers.w12]
        q4 = cci(q2,q3,holynumbers.w10,holynumbers.w02)
        q5 = cci(q3,q2,holynumbers.w03,holynumbers.w11)
        q6 = cci(q4,q3,holynumbers.w05,holynumbers.w04)
        q7 = cci(q6,q5,holynumbers.w06,holynumbers.w07)
        q8 = cci(q7,q5,holynumbers.w08,holynumbers.w09)
        return new(q1,q2,q3,q4,q5,q6,q7,q8)
    end
end

function lxrpt(x;unitlength = 5)
    Luxor.Point(unitlength*x[1],-unitlength*x[2])
end


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
@benchmark ptss = [PointPositions(2π*i/n,holynumbers) for i in 0:n-1]
locus=[ptss[i].q8 for i in 1:n]


function draw(filename, pts::PointPositions, locus)
    q1 = pts.q1
    q2 = pts.q2
    q3 = pts.q3
    q4 = pts.q4
    q5 = pts.q5
    q6 = pts.q6
    q7 = pts.q7
    q8 = pts.q8

    Luxor.Drawing(800, 800, filename)
    Luxor.origin()
    Luxor.translate(-200,-150)
    Luxor.background(RGB(0.2,0.2,0.2))

    Luxor.sethue("red")
    for i in 1:length(locus)-1
        Luxor.line(lxrpt(locus[i]),lxrpt(locus[i+1]),:stroke)
    end
    Luxor.line(lxrpt(locus[n]),lxrpt(locus[1]),:stroke)

    Luxor.sethue("gray")
    for (q_s,q_e) in [(q1,q2), (q2,q4), (q3,q4), (q2,q5), (q3,q5), (q4,q6), (q3,q6), (q6,q7), (q5,q7), (q5,q8), (q7,q8)]
        Luxor.line(lxrpt(q_s),lxrpt(q_e),:stroke)
    end

    Luxor.sethue("black")
    for q in [q1,q2,q3,q4,q5,q6,q7,q8]
        Luxor.circle(lxrpt(q),5,:fill)
    end

    Luxor.finish()
    Luxor.preview()
end


for i in 1:n
    draw("linkage$(lpad(i,3,"0")).png", ptss[i], locus)
end
