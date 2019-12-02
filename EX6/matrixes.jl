include("standard_start.jl");
using LinearAlgebra;

v1 = [2;1]; #v1 = v1/norm(v1)
v2 = [0.5;2.5]; #v2 = v2/norm(v2)

V = [v1 v2]

x1 = 2*v1
x = [x1 v2]

x_new = inv(V)*x
Lambda = collect(Diagonal(x_new))


M = V*Lambda*inv(V)

# STRETCH & SQUEEZE
stretch = 1.3
squeeze = 0.8

ss = collect(Diagonal([stretch,squeeze]))
M_ss = M*ss

M_orig = M_ss*Lambda*inv(M_ss);

eig = eigen(M_orig)

M_orig*[0.967075; 0.254493]
2*[0.967075; 0.254493]

M_orig*[-0.242536;0.970143]
1*[-0.242536;0.970143]
