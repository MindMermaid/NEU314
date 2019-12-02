include("standard_start.jl")
using LinearAlgebra
using Printf
##===== 5a
M = [0 1;-1 0];
    M_eig = eigvals(M);
    println("M eigenvalues =  [", M_eig[1], "]\n                 [", M_eig[2], "]");
# this doesn't work for some reason:  diag(eigvals(M))
 lambda = [M_eig[1] 0; 0  M_eig[2]];
 V = eigvecs(M);
 V_i = inv(V)
 println("is   M == V*lambda*V_i  true?   ------>  ", M == V*lambda*V_i)
#

##===== 5b
M1_2 = V*sqrt(lambda)*V_i
println("M1_2  =   [", M1_2[1,1], "   ", M1_2[1,2], "]\n         [", M1_2[2,1], "   ", M1_2[2,2], "]")
#

#=== 5c
    theta in radians for M = pi/2
    theta in radians for M1_2 = pi/3
=#

##====== 5d
include("animate_matrix.jl")
animate_matrix(M) # 90 degree rotation
animate_matrix(M1_2) # 45 degree rotation

#= M & M1_2 are rotational matrixes - M rotates points by 90 degrees
clockwise, M1_2 rotates them by 45 degrees clockwise =#


##===== 4
d = cos(pi/12)+1.0im*sin(pi/12)
d2 = d*d
d3 = d2*d

d_x = [real(d) real(d2)  real(d3)]
d_y = [imag(d) imag(d2) imag(d3)]

plot(d_x,d_y,"o");
    title("d, d*d, d*d*d")
    xlabel("real")
    ylabel("imaginary")
    grid("true")
    ylim(0,1)
    xlim(0,1)
#

e = MathConstants.e;
e1 = e^(1.0im*pi/12);
e2 = e1*e1
e3 = e2*e1

e_x = [real(e1) real(e2)  real(e3)]
e_y = [imag(e1) imag(e2) imag(e3)]

plot(e_x,e_y,"o");
    title("e, e*e, e*e*e")
    xlabel("real")
    ylabel("imaginary")
    grid("true")
    ylim(0,1)
    xlim(0,1)
#

println("is e3 = d3  ?   ---->  ", e3==d3)

println("atan((sin(pi/4)/cos(pi/4))) =  ", atan((sin(pi/4)/cos(pi/4))))
