# PREP

using PyPlot
using LinearAlgebra
include("gsample.jl")
include("animate_matrix.jl")

v1 = [0.4;1]; v1 = v1/norm(v1)
v2 = [-0.4;1]; v2 = v2/norm(v2)

V = [v1 v2]
Lambda = collect(Diagonal([0.5,2]))

M = V*Lambda*inv(V)

b1 = [2;1]
b2 = [-1;1]
B = [b1 b2]


##  B moves cartesian axes to its columns:
animate_matrix(B, seed=zeros(2,0), vectors=[1 0;0 1], titleString="B")

##  inv(B) moves its columns back to the cartesian axes:
animate_matrix(inv(B), seed=zeros(2,0), vectors=B, titleString="inv(B)")

## So we casn use inv(B) to express x in basis set given by B's columns:

single_point = [0.5 ; 2]
D1 = animate_matrix(inv(B), plot_tails=false, vectors=B, seed=[0.5 ; 2],
    titleString="inv(B)")


##  Let's look at the transform induced by a matrix:

D = animate_matrix(M, titleString="M", titleColor="green")


##  Add some cartesian axes to it -- we're going to rotate them next

animate_matrix([1 0 ; 0 1], seed=D["data"], plot_tails=false, vectors=[1 0 ; 0 1],
    linesX=D["tailsX"], linesY=D["tailsY"], titleString="M", titleColor="green")

##  If we change bases by rotating, transform M still exists-- but we now write
# it down with new numbers

theta = 70*pi/180
R = [cos(theta) -sin(theta) ; sin(theta) cos(theta)]

animate_matrix(R, seed=D["data"], plot_tails=false, vectors=[1 0 ; 0 1],
    linesX=D["tailsX"], linesY=D["tailsY"], titleString="inv(R) M R", titleColor="green")


## Another example of changing basis, this time not a rotation

D = animate_matrix(M, titleString="M", titleColor="green")

D2 = animate_matrix([1 0 ; 0 1], seed=D["data"], plot_tails=false, vectors=B,
    linesX=D["tailsX"], linesY=D["tailsY"], titleString="M", titleColor="green")

## Change basis to columns of B

animate_matrix(inv(B), seed=D["data"], plot_tails=false, vectors=B,
    linesX=D["tailsX"], linesY=D["tailsY"], titleString="inv(B) M B", titleColor="green")


##################
# BACK TO IPAD: how does a change of basis affect a matrix?
##################

##  Let's show Mnew = inv(B) M B directly

figure(2)
Mnew = inv(B)*M*B

seed = inv(B)*inv(M)*D["data"]

animate_matrix(Mnew, seed=seed)
figure(1)

##  Could do change to other bases -- start by showing M again and new basis vectors
B2 = [[1;-2] [1;1]]

animate_matrix([1 0 ; 0 1], seed=D["data"], plot_tails=false, vectors=B2,
    linesX=D["tailsX"], linesY=D["tailsY"], titleString="M", titleColor="green")

##

animate_matrix(inv(B2), seed=D["data"], plot_tails=false, vectors=B2,
    linesX=D["tailsX"], linesY=D["tailsY"], titleString="inv(B2) M B2", titleColor="green")



## ========  Now, how about we change to M's eigenvector basis?

D = animate_matrix(M, titleString="M", titleColor="green", npoints=300)


animate_matrix([1 0 ; 0 1], seed=D["data"], plot_tails=false, vectors=V,
    linesX=D["tailsX"], linesY=D["tailsY"], titleString="M", titleColor="green")

##

animate_matrix(inv(V), seed=D["data"], plot_tails=false, vectors=V,
    linesX=D["tailsX"], linesY=D["tailsY"], titleString="inv(V) M V", titleColor="green")


########################
# BACK TO IPAD: changing to eigenvector basis diagonalizes
########################


# =============================================
#
#        PCA
#
# =============================================


C = [2 1 ; 1 1]
npoints = 300

data = gsample(C, npoints=npoints)

projection(data)


##

theta = 10*pi/180
v = [cos(theta) ; sin(theta)]

projection(data, v, plotProjDots=false, plotProjConn=false)
##

theta = 130*pi/180
v = [cos(theta) ; sin(theta)]
z = projection(data, v)
varz = sum(z.^2)/length(z)
title("variance is " * string(round(varz, digits=4)) )

##

######################################################
#
# BACK TO IPAD FOR WHICH DIRECTION HAS MOST VARIANCE
#
######################################################

PC = 1

E = eigen(data*data'/size(data,2))
z = projection(data, E.vectors[:,end-(PC-1)])
varz = sum(z.^2)/length(z)
title("variance is " * string(round(varz, digits=4)) * " eval is " * string(round(E.values[3-PC], digits=4)) )


##

PC = 1

E = eigen(data*data'/size(data,2))
Rdata = inv(E.vectors[:,end:-1:1])*data

EE = eigen(Rdata*Rdata'/size(Rdata,2))
z = projection(data, EE.vectors[:,end-(PC-1)])
varz = sum(z.^2)/length(z)

projection(Rdata, EE.vectors[:,end-(PC-1)])
