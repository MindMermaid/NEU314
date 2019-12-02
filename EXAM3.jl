##===== Files
include("standard_start.jl")

##===== Packages
using LinearAlgebra
using Colors
Colorrs = Colors.color_names;


##===== Functions
"""
quiver_plot(dxdt, axstart, axstop; nsteps=20, fignum=1)

Plots a 2d quiver_plot given the dynamics

"""

function quiver_plot(dxdt, axstart, axstop; nsteps=20, fignum=1)
    gridpoints = collect(axstart .+(axstop-axstart)*(0:(nsteps)-1)/(nsteps-1))

    X = ones(length(gridpoints),1)*gridpoints[:]'
    Y = copy(X')

    pts = [X[:]' ; Y[:]']


    figure(fignum); clf();
    h = quiver(pts[1,:], pts[2,:], dxdt(pts)[1,:], dxdt(pts)[2,:],
        color=(0.5, 0.5, 0.5))
    axis("scaled")

    xlim([axstart, axstop]); ylim([axstart, axstop])
    xlabel("x"); ylabel("y")
end

##===== Problem 1
M = [3+5/6 -5/3; 5/3 -1/3];
M_i = inv(M);
M_eig= eigen(M);

##===== Problem 1c
plot_center = [0, 0] # setting the center of the quiver plot
    gridpoints = collect(-10:1:20) # setting the plot siVe and the density of arrows
    grid_points_2D = ones(length(gridpoints),1)*gridpoints[:]'
    X = grid_points_2D  .+ plot_center[1]
    Y = grid_points_2D' .+ plot_center[2]
    X = X[:]; Y = Y[:]; # make X and Y into one-dimensional arrays

    U = (3+5/6) .* X + (-5/3).* Y .-2; # compute x_dot for each point
    V = (5/3) .*X +(-1/3) .* Y .-3; # compute y_dot for each point

    figure("EXAM3_1c");
        quiver(X, Y, U, V);plot(26/9, 49/9, "ro");
        xlim(-10,20);ylim(-10,20);
        hlines(0, -10, 20, color="grey")
        vlines(0, -10, 20, color="grey")
        title("Quiver Plot"); xlabel("x"); ylabel("y");
    #
#
##===== Problem 2a
function f(V)
    return (-V ./ (1 .+ 0.15.*ℯ.^(-0.08 .* V))) .- (0.2 .* V) .- 16
end

function dVdt(V)
    dV1dt = f(V[2,:]) .- V[1,:]
    dV2dt = f(V[1,:]) .- V[2,:]
    return [ dV1dt' ; dV2dt' ]
end

axstart = -100
axstop = 0
V=axstart:0.01:axstop
figure("EXAM3_2a")
    plot(V, f(V), label="V_dot(V)", color="darkturquoise");
    plot(0, -73.1, color="r");
    title("V_dot(V)")
    xlim(-100, 0); ylim(-20, 10);grid(.1);
    xlabel("V (mV)");ylabel("V_dot"); legend();
#
##=====

##===== Problem 2c
dt = 0.01; t = 0:dt:50;
    V = zeros(2,length(t))
    V[:,1] = [0,-44];

    for i=1:length(t)-1
        V[:,i+1] = V[:,i] + (-16.8 / (1 + (-0.08*0.15)*ℯ^(0.08*-8.1)) .* (V[:,i] .+ 8.1));
    end
#

dt = 0.01; t = 0:dt:50
vol = zeros(2,length(t));
    vol[2,1] = -45
    vol[1,1] = -44 # Need to define a starting point

for i=1:length(t)-1
    vol[1,i+1] = vol[1,i] + dt*dVdt(vol[1,i])
    vol[2,i+1] = vol[2,i] + dt*dVdt(vol[2,i])
end

figure("EXAM3_2c");
    plot(t, vol[1, :], color="g", label="V0=-45V");
    plot(t, vol[2,:], color="b", label="V0=-44V");
    plot(t[end], vol[1,end],"o", color="g", label="V_end = -21.53V")
    plot(t[end], vol[2,end],"o", color="b", label="V_end = -70.13V")
    xlabel("t");ylabel("V");title("V(t)");legend()
    xlim(0, 51); ylim(-72, -20);
#

##=====Problem 3a

function dVdt(V, W, a, b, I)
    return V .- V.^3 .- W +  I;
end
function dWdt(V, W, a, b)
    return b .*V .+ a .- W;
end

a = 2; b = 2; I = 0;

dt = 0.0001; t = 0:dt:50;
    x = zeros(2,length(t))
    x[1,1] = -1;
    x[2,1] = -1.9;

    for i=1:length(t)-1
        x[1,i+1] = x[1,i] + dt*dVdt(x[1,i], x[2,i], a, b, I)
        x[2,i+1] = x[2,i] + dt*dWdt(x[1,i], x[2,i], a, b)
    end

figure("EXAM3_3")
    plot(t, x[1,:], color = "lightgreen",label="V(t)")
    grid();title("Membrane Potential with I = 0")
    xlabel("t"); ylabel("V");legend();
#

figure("EXAM3_3")
    plot(x[1,:],x[2,:], color="darkblue")
    plot(x0, v0, "o", color="darkturquoise", label="F.P --> [V,W]=[-1.2879,-0.57582]")
    grid()
    xlabel("V"); ylabel("W");title("Fitzhugh-Nagumo System");legend()
#

V_0 = -1.2879 - ((-1.2879)^3/3) + 0.57582;
W_0 = 2*(-1.2879) + 2 + 0.57582;

V_0 = 0 - ((0)^3/3) - 2 + 2;
W_0 = 2*(0) + 2 - 2;


x0=-1.2879;
v0=-0.57582;

function jacobian(V,W,a,b,I)
    return [(1 - V^2 - W + I) (V-(V^3)/3 -1 + I); (b + a - W) (b*V + a -1)];
end

J1 = jacobian(-1.2879,-0.57582, 2, 2, 0)
J2 = jacobian(0,2,2,2,2)

eigenvalues =[eigen(J1).values eigen(J2).values];
