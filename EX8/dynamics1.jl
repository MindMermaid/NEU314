using PyPlot

##  =========================
#
#     EXAMPLE OF EULER INTEGRATION in 1-d NONLINEAR DYNAMICS
#
#   =========================

function dxdt(x)
    return 8*cos.(x).^2 + 2*x .- 3.5
end

figure(1); clf()

x=-3*pi:0.1:3*pi
plot(x, dxdt(x))
hlines(0, xlim()[1], xlim()[2], color="black", linewidth=1)
vlines(0, ylim()[1], ylim()[2], color="black", linewidth=1)
grid("on")
xlabel("x"); ylabel("xdot"); title("xdot(x)")
##

dt = 0.05
t = 0:dt:2
x = zeros(size(t))
x[1] = -0.65   # Need to define a starting point

for i=1:length(t)-1
    x[i+1] = x[i] + dt*dxdt(x[i])
end

figure(2); clf();
plot(t, x, "o-")
xlabel("t"); ylabel("x") ; title("x(t)")


# ======================================
#
#   BACK TO BOARD
#
# ======================================


##

# ======================================
#
#   EULER INTEGRATION AND TWO DIMENSIONAL DYNAMICS  -- 2-neuron network
#
# ======================================

M = [-1 1; -1 0.5]

dt = 0.05; t = 0:dt:20
x = zeros(2,length(t))
x[:,1] = [-1,1]

for i=1:length(t)-1
    x[:,i+1] = x[:,i] + dt*M*x[:,i]  # Euler recipe
end

figure(1); clf(); axlim = 2
xlim([-axlim, axlim]); ylim([-axlim, axlim])
hlines(0, -axlim, axlim, color="black")
vlines(0, -axlim, axlim, color="black")
axis("scaled"); grid("on")

plot(x[1,:], x[2,:], color="green")
# plot(x[1,:], x[2,:], "d")
plot(x[1,1], x[2,1], "o", color="green")
plot(x[1,end], x[2,end], "o", color="red")
xlabel("x1"); ylabel("x2")

# saveNmove("twod_euler_1")

# ======================================
#
#   BACK TO BOARD
#
# ======================================


##

# ======================================
#
#   Quiver plot
#
# ======================================

# M = [-1.6 2.5 ; -2.5 2.1]

include("plot_arrow.jl")

gridpoints = collect(-2:0.25:2)
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
xlabel("x1"); ylabel("x2")



# saveNmove("twod_quiver_2")
