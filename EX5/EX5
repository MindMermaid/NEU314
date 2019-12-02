include("standard_start.jl");
include("animate_matrix.jl");
##=====1a
x = [x1 x2 x3]; # each column of x will be one of x1, x2, or x3
    x1 = [0.5; 1];
    x2 = [0.8; 0.2];
    x3 = [0.7; 0.3];


    M = [2 0; 0 2];

    y1 = M*x1;
    y2 = M*x2;
    y3 = M*x3;

    println(y1, y2, y3);

figure(1)
    animate_matrix(M);
    clf()
    animate_matrix(M, seed = x); # seed=x tells it which specific points to show
##=====

##===== 1b
M2 = [-3 1; 0 1];

    z1 = M2*x1
    z2 = M2*x2
    z3 = M2*x3

    println(z1, z2, z3);
##==== 1b_figure
animate_matrix(M2);
animate_matrix(M2,seed = x);
##=====

##====Part II _ 1c
theta = 2*pi/3;
    R = [cos(theta) sin(theta); -sin(theta) cos(theta)];

    y1 = R*x1
    r2 = R*x2
    r3 = R*x3
    r = [r1 r2 r3];

    R_i = inv(R);

    RR = R*R;

println(r1, r2, r3);
println("R = ", R, " & R_i = ", R_i, " RR = ", RR);

##===== 1c_ animate rotations by 120 degrees
    animate_matrix(R);
    animate_matrix(R, seed=x);
##=====
##===== 1d _ undo rotations
animate_matrix((inv(R),R));
animate_matrix((R,R,R));
##=====

##==== 1e
A = M2*R
    B = R*M2
    println(A,B);

animate_matrix((M2,R));
animate_matrix(A, seed= x);
animate_matrix((R,M2));
animate_matrix(B, seed = x);
##====

##====1f
M3 = [0.5 0; 0 2];
M4 = [1.25 0.75; 0.75 1.25];
animate_matrix(M3);
animate_matrix(M4);
#= EXPLANATION
M3 squashes the vectors along the x-axis
M4 stretches the vectors along the diagonal=#

##===== 1g
V3 = [1/sqrt(2) -1/sqrt(2); -1/sqrt(2) -1/sqrt(2)];
animate_matrix((V3,M3,V3));
animate_matrix(M4);
##=====


##=====2b
A = [4 -2; 1 1]
    x = [2;1]

    using LinearAlgebra

    println(eigvecs(A) * norm(x));
##======

##====== 3a
    A = [1 0.25; 0.5 1.5];
    B = [1 -.5; .5 -.25];

    animate_matrix(A)
    animate_matrix(B)

    println("det(A) = ",det(A));
    println("det(B) = ",det(B));

# B is invertible because det(B) = 0
##=====

##===== 3b
A = [1 2; -1 .5];
B = [-2 -4; -3 1];
println("det(AB) = ",det(A*B));
println("det(BA) = ",det(B*A));
println("det(AB) = det(A) * det(B) true or false?   ==", det(A)*det(B) == det(A*B));

 #= volume does not change along multiple transformations
 as we have shown that det(AB) = det(BA). Determinants are
 COMMUTATIVE =#
##======
