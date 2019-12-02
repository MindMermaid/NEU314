
using LinearAlgebra


"""
plot_arrow(start, stop; hlen=0.2, theta=30, old_arrow_handle=nothing,
    coords_only=false, kwargs...)

plots an arrow. The length of the arrow head is a scaled version of
the length of the body; the angle is theta. kwargs are arguments that
are passed on to PyPlot (e.g., color="red" etc.)

= PARAMETERS:

- start   a 2-dimensional vector or an n-by-2 matrix (for n arrows)

- stop    a 2-dimensional vector or an n-by-2 matrix (for n arrows)

= OPTIONAL PARAMETERs:

- hlen   the scaling factor for the arrowhead

- theta   angle, in degrees of arrowhead with respect to body

- old_arrow_handle if this is passed, it must be a handle to a PyPlot Line2D
          object, whose xdata and ydata will be updated to the arrow coords
          in start and stop. New kwargs are ignored, since only the xdata
          and ydata will be updated.

- coords_only   If this is true, function does not plot, but returns LX, LY
        two vectors with complete X coords and complete Y coords for plotting the
        set of arrows

= RETURNS:

- h     A Line2D object that has all the arrows in a single line
        Each stroke is separated from the next one by a NaN. Thus
        each arrow, with three strokes, has a start and a stop and
        a NaN for each stroke and its xdata will have length 9*n where
        n is the number of arrows.

if optional parameter coords_only is passed as true, then instead returns
LX, LY two vectors with complete X coords and complete Y coords for plotting the
set of arrows.


= EXAMPLE CALL:

clf()
plot_arrow([0,0], [1,2], linewidth=4)
axis("scaled")


"""
function plot_arrow(start, stop; hlen=0.2, theta=30, old_arrow_handle=nothing,
    coords_only=false, kwargs...)

    if !all(size(start) .== size(stop))
        error("start and stop must have the same size")
    end
    if length(start)>0
        if (length(size(start))==1 && length(start) != 2) || (length(size(start))==2 && size(start)[2] != 2)
            error("start and stop must have two columns or be 2-long vectors")
        end
    end

    theta = theta * pi/180

    if length(start)==0
        start = zeros(0,2)
        stop  = zeros(0,2)
    elseif length(size(start))==1
        start = (Z = zeros(1,2); Z[:] = start; Z)
        stop  = (Z = zeros(1,2); Z[:] = stop;  Z)
    end

    LX = zeros(9*size(start)[1])
    LY = zeros(9*size(start)[1])

    function rotate(vec, theta)
        R = [cos(theta) sin(theta); -sin(theta) cos(theta)]
        return R*vec
    end

    for i=1:size(start)[1]
        me = (i-1)*9+1
        LX[me:me+2] = [start[i,1], stop[i,1], NaN]
        LY[me:me+2] = [start[i,2], stop[i,2], NaN]

        delta = stop[i,:] - start[i,:]
        headlen = norm(delta)*hlen
        if abs(delta[1]) <= 10*eps()
            if delta[2] > 0; rangle = pi/2
            else             rangle = -pi/2
            end
        else
            rangle = atan(delta[2]/delta[1])
            if delta[1]<0
                rangle = rangle+pi;
            end
        end

        head = rotate([headlen;0], -rangle+theta+pi)
        LX[me+3:me+5] = [stop[i,1], stop[i,1]+head[1], NaN]
        LY[me+3:me+5] = [stop[i,2], stop[i,2]+head[2], NaN]

        head = rotate([headlen;0], -rangle-theta+pi)
        LX[me+6:me+8] = [stop[i,1], stop[i,1]+head[1], NaN]
        LY[me+6:me+8] = [stop[i,2], stop[i,2]+head[2], NaN]
    end

    if coords_only
        return LX, LY
    end

    if old_arrow_handle == nothing
        return plot(LX, LY; kwargs...)[1]
    else
        old_arrow_handle.set_xdata(LX)
        old_arrow_handle.set_ydata(LY)
        return old_arrow_handle
    end

end
