include("standard_start.jl")

##===== Problem 1
function dxdt(x)
    return -2x
end

function f(x)
    return x.^(2)
end

dt = 0.05
t = 0:dt:5
x = zeros(size(t))
x[1] = 5   # Need to define a starting point

for i=1:length(t)-1
    x[i+1] = x[i] + dt*dxdt(x[i])
end

figure(2); clf();
    plot(t, x, "o-")
    plot(f(t), x)
    xlabel("t"); ylabel("x") ; title("x(t)")
##=====



##===== Problem 3
using Roots
function dxdt(x)
    return (x.^4) - (2x.^3) - (5x.^2) + 6*x
end

figure(1); clf()

    x=-4:0.1:4
        plot(x, dxdt(x), label = "xdot", color="turquoise")
        hlines(0, xlim()[1], xlim()[2], color="black", linewidth=1)
        vlines(0, ylim()[1], ylim()[2], color="black", linewidth=1)
        grid("on")
        xlabel("x"); ylabel("xdot"); title("xdot(x)")
    #
    plot(x, x_dot(x), label = "xdot ~ @x=1", color="orange")
    dxdt_0 = find_zeros(dxdt, -4, 4)
        stable = [dxdt_0[1] dxdt_0[3]]
        unstable = [dxdt_0[2] dxdt_0[4]]
        plot(stable, dxdt(stable), ".", label = "Stable FP", color="red")
        plot(unstable, dxdt(unstable), ".", label = "Unstable FP", color="blue")
    #
    legend()
##

function x_dot(x)
    return -6 .* x .+ 6
end
    x=.5:0.1:1.5
        plot(x, x_dot(x), color="orange")
    #
#

##===== Problem 4
using Colors
using LinearAlgebra
include("animate_matrix.jl")

M = [.5 -1; 1 -0.5]

M_eig = eigen(M)

V = eigvals(M)

dt = 0.01; t = 0:dt:200
x = zeros(2,length(t))
x[:,1] = [-.1,.1]

for i=1:length(t)-1
    x[:,i+1] = x[:,i] + dt*M*x[:,i]  # Euler recipe
end

figure(1); clf(); axlim = 3
    xlim([-axlim, axlim]); ylim([-axlim, axlim])
    hlines(0, -axlim, axlim, color="black")
    vlines(0, -axlim, axlim, color="black")
    axis("scaled"); grid("on")

    plot(x[1,:], x[2,:], color="darkturquoise")
    # plot(x[1,:], x[2,:], "d")
    plot(x[1,1], x[2,1], "x", color="slategrey", label="start")
    plot(x[1,end], x[2,end], "o", color="mediumturquoise", label="stop")
    xlabel("rE"); ylabel("rI")
    legend()
# ======================================
#
#   Quiver plot
#
# ======================================
include("plot_arrow.jl")

gridpoints = collect(-3:0.25:3)
dt = 0.1

X = ones(length(gridpoints),1)*gridpoints[:]'
Y = copy(X')

pts = [X[:]' ; Y[:]']

pts2 = pts + dt*M*pts

figure(1); clf();
    plot_arrow(pts', pts2', color=[0.75, 0.75, 0.75], linewidth=1)
    axis("scaled")
    xlim([-axlim, axlim]); ylim([-axlim, axlim])
    hlines(0, -axlim, axlim, color="black")
    vlines(0, -axlim, axlim, color="black")
    xlabel("rE"); ylabel("rI")
