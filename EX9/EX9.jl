##===== Files
include("standard_start.jl")
include("animate_matrix.jl")

##===== Packages
using LinearAlgebra
using Colors
Colorrs = Colors.color_names;
##===== Problem 1 Part 2

##===== 2-Neuron System
# Numerical
M = [0.3 -1.2; 1 -0.8];
    dt = 0.001; t = 0:dt:20;
        x = zeros(2,length(t))
        x[:,1] = [0.001,-0.001]

        for i=1:length(t)-1
            x[:,i+1] = x[:,i] + dt*M*x[:,i]  # Euler recipe
        end

        figure("EX9_1.2d1"); clf();
            plot(t, x[1,:], color="darkturquoise", label="x")
            plot(t, x[2,:], color="mediumvioletred", label="y")

            hlines(0, 0, 20, linestyle = "--", color="black")
            xlabel("t"); ylabel("values");grid("on");xlim(0,20)
            title("2 Neuron System (Numerical)");legend()
        #
    #

# Analytical
eigenspace = eigen(M).vectors
    eigenvalues = eigen(M).values
    initial_xy = [0.001,-0.001];
    uv_o = inv(eigenspace) * initial_xy
        # convert the initial values of x and y into the eigenspace of M
        # this will change the coordinates of the initial values from the standard basis set to the new basis set (the two eigenvectors of M)

        uv_t = zeros(length(initial_xy), length(t)) # make a matrix to store the values of the analytical solution in the eigenspace of M
        uv_t = complex(uv_t) # make the matrix able to accept complex numbers

        for i = 1:length(initial_xy)
            uv_t[i, :] = uv_o[i] * exp.(t * eigenvalues[i]) # compute the analytical solution for each dimension in the eigenspace of M
        end
        xy_t_analytical =  eigenspace * uv_t; # convert the resulting analytical solution in the eigenspace of M into the original Cartesian space

        figure("EX9_1.2d2"); clf();
            plot(t, xy_t_analytical[1,:], color="darkblue", label = "x")
            plot(t, xy_t_analytical[2,:], color="plum", label = "y")
                hlines(0, 0, 20, linestyle = "--", color="black")
                xlabel("t"); ylabel("values");grid("on");xlim(0,20)
                title("2 Neuron System (Analytical)");legend()
            #
        #
    #

fig = figure("EX9_1.2d3",figsize=(10,10));
    subplot(211)
    plot(t, x[1,:], color="darkturquoise", label="x")
    plot(t, x[2,:], color="mediumvioletred", label="y")

    hlines(0, 0, 20, linestyle = "--", color="black")
    xlabel("t"); ylabel("values");grid("on");xlim(0,20)
    title("Numerical");legend()

    subplot(212)
        plot(t, xy_t_analytical[1,:], color="darkblue", label = "x")
        plot(t, xy_t_analytical[2,:], color="plum", label = "y")
            hlines(0, 0, 20, linestyle = "--", color="black")
            xlabel("t"); ylabel("values");grid("on");xlim(0,20)
            title("Analytical");legend()
        #
    #
figure("EX9_1.2e",figsize=(8.5,7)); clf();
        xlim([-0.001, 0.002]); ylim([-0.0011, 0.0012]);
        hlines(0, -0.001, 0.002, color="black")
        vlines(0, -0.0011, 0.0012, color="black")
         grid("on")

        plot(x[1,:], x[2,:], color="darkturquoise", label= "[x(t),y(t)]")
        # plot(x[1,:], x[2,:], "d")
        plot(x[1,1], x[2,1], "x", color="slategrey", label="start")
        plot(x[1,end], x[2,end], "o", color="red", label="stop")
        xlabel("x"); ylabel("y"); title("[x,y] time course")
        legend()
##=====

##===== 3-Neuron System
M = [-0.5 -0.5 0; -0.5 0.1 -1; -1 0 -1]

dt = 0.001; t = 0:dt:20
    x = zeros(3,length(t))
    x[:,1,1] = [0.001,-0.001,0.002]

    for i=1:length(t)-1
        x[:,i+1] = x[:,i] + dt*M*x[:,i]  # Euler recipe
    end

    figure("EX9_1.2d4"); clf();
        plot(t, x[1,:], color="darkturquoise", label="x")
        plot(t, x[2,:], color="mediumvioletred", label="y")
        plot(t, x[3,:], color="lightgreen", label="z")

        hlines(0, 0, 20, linestyle = "--", color="black")
        xlabel("t"); ylabel("values");grid("on");xlim(0,20)
        title("3 Neuron System (Numerical)");legend()
    #


## Analytical
eigenspace = eigen(M).vectors
eigenvalues = eigen(M).values
initial_xy = [0.001,-0.001,0.002];
    uv_o = inv(eigenspace) * initial_xy

        uv_t = zeros(length(initial_xy), length(t))
        uv_t = complex(uv_t)

        for i = 1:length(initial_xy)
            uv_t[i, :] = uv_o[i] * exp.(t * eigenvalues[i])
        end
        xy_t_analytical =  eigenspace * uv_t;

        figure("EX9_1.2d5"); clf();
            plot(t, xy_t_analytical[1,:], color="darkblue", label = "x")
            plot(t, xy_t_analytical[2,:], color="plum", label = "y")
            plot(t, xy_t_analytical[3,:], color="chartreuse", label = "y")
            hlines(0, 0, 20, linestyle = "--", color="black")
            xlabel("t"); ylabel("values");grid("on");xlim(0,20)
            title("2 Neuron System (Analytical)");legend()
        #
    #
figure("EX9_1.2d6",figsize=(10,10));
    subplot(211)
        plot(t, x[1,:], color="darkturquoise", label="x")
        plot(t, x[2,:], color="mediumvioletred", label="y")
        plot(t, x[3,:], color="lightgreen", label="z")
            hlines(0, 0, 20, linestyle = "--", color="black")
            xlabel("t"); ylabel("values");grid("on");xlim(0,20)
            title("Numerical");legend()
    subplot(212)
        plot(t, xy_t_analytical[1,:], color="darkblue", label = "x")
        plot(t, xy_t_analytical[2,:], color="plum", label = "y")
        plot(t, xy_t_analytical[3,:], color="chartreuse", label = "y")
            hlines(0, 0, 20, linestyle = "--", color="black")
            xlabel("t"); ylabel("values");grid("on");xlim(0,20)
            title("Analytical");legend()
        #
    #
##=====
