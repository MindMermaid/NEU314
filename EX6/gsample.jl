"""
gsample(C; npoints=1)

returns samples from a 2-dimensional Gaussian with mean zero
and given covariance matrix C

= PARAMETERS:

- C    a 2-by-2 covariance matrix

= OPTIONAL PARAMETERS:

- npoints    number of points to sample

= RETURNS:

- data       a 2-by-npoints matrix of the sampled points

"""
function gsample(C; npoints=1)
    D = randn(2, npoints)
    return sqrt(C)*D
end


"""
projection(data, v; fig=1, dotColor="green",
    paxColor=[0.7, 0.7, 0.7],
    projLineColor="blue", projLineWidth=2,
    projDotColor = [0.6, 0, 0],
    connColor=[0.5, 0.5, 0.5], connLineStyle="--", connLineWidth=0.5,
    plotProjLine=true, plotProjDots=true, plotProjConn=true)

plots and projects a set of data points onto a given line

= PARAMETERS:

- data    a 2-by-npoints matrix, indicating points to be plotted

- v       a vector of 2 elements, indicating the direction of the projection line.
            The geometric length of this vector is ignored.

= OPTIONAL PARAMETERS:

- fig     The figure number to plot in

- dotColor   the color of the dots plotted

- paxColor   the color of the principal axes

- projLineColor   the color of the projection line (v)

- projLineWidth   the linewidth of the projection line

- projDotColor    the color with which dots projected onto the line will be shown

- connColor       the color of the lines connecting original dots to their projection

- connLineStyle   the line style of the lines connecting original dots to their projection

- connLineWidth   the line width of the lines connecting original dots to their projection

- plotProjLine    Boolean, determines whether the projection line (v) is shown or not

- plotProjDots    Boolean, determines whether the projected dots are shown or not

- plotProjConn    Boolean, determines whether connecting lines are shown or not

= RETURNS:

- projs   a 1-by-npoints array of the projection distances along v

"""
function projection(data, v; fig=1, dotColor="green",
    paxColor=[0.7, 0.7, 0.7],
    projLineColor="blue", projLineWidth=2,
    projDotColor = [0.6, 0, 0],
    connColor=[0.5, 0.5, 0.5], connLineStyle="--", connLineWidth=0.5,
    plotProjLine=true, plotProjDots=true, plotProjConn=true)

    v = copy(v)/norm(v); v = v[:]

    figure(fig); clf();
    plot(data[1,:], data[2,:], ".", color=dotColor)
    axis("equal")
    PyPlot.show; gcf().canvas.flush_events()

    vlines([0], ylim()[1], ylim()[2], color=paxColor, linestyle="--")
    hlines([0], xlim()[1], xlim()[2], color=paxColor, linestyle="--")

    PyPlot.show; gcf().canvas.flush_events()
    xl = xlim(); yl = ylim();

    if plotProjLine
        myv = 2 .* maximum([diff(collect(xlim())), diff(collect(ylim()))]) .* v;

        plot([-myv[1], myv[1]], [-myv[2], myv[2]], color=projLineColor,
            linewidth = projLineWidth,
            scalex=false, scaley=false)
    end

    pdata = v*(data'*v)'
    if plotProjDots
        plot(pdata[1,:], pdata[2,:], ".", color=projDotColor)
        xlim(xl); ylim(yl)
    end

    if plotProjConn
        Xconn = zeros(3*size(data,2))
        Yconn = zeros(3*size(data,2))
        for i=1:size(data,2)
            pos = (i-1)*3+1
            Xconn[pos]   = data[1,i]
            Yconn[pos]   = data[2,i]
            Xconn[pos+1] = pdata[1,i]
            Yconn[pos+1] = pdata[2,i]
            Xconn[pos+2] = NaN
            Yconn[pos+2] = NaN
        end
        plot(Xconn, Yconn, color=connColor, linestyle=connLineStyle,
            linewidth=connLineWidth)
            xlim(xl); ylim(yl)
    end

    return (data'*v)'
end

"""
projection(data)

Like projection(data, v), but no projection line, projected dots, or connecting
lines are shown. Other optional arguments to projection are allowed.
"""
function projection(data; kwargs...)
    return projection(data, [1;1],
        plotProjLine=false, plotProjDots=false, plotProjConn=false)
end

# projection(X[:,1:300], [2;1])
