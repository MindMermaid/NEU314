using JLD
using LinearAlgebra
using PyPlot
using Random


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


"""
animate_matrix(M::Array ; nsteps=30,
    seed=:random, mx=3, npoints=100, principal_axes=false,
    pax_linewidth = 0.5, pax_linecolor = (0.7, 0.7, 0.7),
    plot_tails = true, tailscolor=:green,
    markercolor=:mediumorchid, markersize=5,
    filename="", vectors=zeros(2,0), highpts=zeros(2,0),
    linesX = zeros(0,0), linesY = zeros(0,0), linescolor=:green, lineswidth=1,
    userpause=false, pause_steps=1, grid=true,
    pausetime=0.025, first_pause_time=1.0,
    anim=nothing, fps=10, srand_seed=nothing, plot_kwargs...)

Given a 2-by-2 matrix M, illustrates the effects of M by
taking many random points on the 2d plane and showing how M
transforms them. In addition, can also show and transform
highlighted points (shown in red); vectors (shown as black arrows
from the origin); and lines (arbitrary X,Y sequences, connected by
straight lines)

# PARAMETERS:

- M    a 2-by-2 matrix whose transformation will be illustrated; OR
       a tuple of 2-by-2 matrices, whose transformation will be shown
       in sequence, the last one first, and the first one last
       (this is so animate_matrix(M1*M2) produces the same as
       animate_matrix((M1, M2))).

# OPTIONAL PARAMETERS:

- nsteps    The number of steps to take in the animation.
            The more steps, the slower.

- seed      Currently can only take the value :random, or be a
            2-by-n matrix, or be a Dict. If it is :random, this means
            points will appear at random position on the
            plane. If it is a 2-by-n matrix, points will appear at each
            column of the matrix. If it is a Dict, it is expected to have
            the same keys and data format as the output of animate_matrix().
            This is useful for manually chaining calls to animate_matrix().

- mx        Only points within this radius will be shown. The
            limits of the plot are also set by mx.

- npoints   The number of points to show

- principal_axes  Boolean. If true, shows grey vertical and
            horizontal lines at zero.

- pax_linewidth   Linewidth of the principal axes, if they are shown

- pax_linecolor   Color of the prinicpal axes, if they are shown

- plot_tails.    Boolean. If true shows tails of the points,
            otherwise shows only points.

- tailscolor    The color of the tails (if they are plotted)

- markercolor   The color of the points.

- filename      String. By default empty string. If not
           empty, then a .mov file will be generated
           with that filename, showing a movie of the
           animation.

- vectors   A 2-by-n matrix, where n=0 is a legal value. Each column
            will be shown as a black vector arrow, base at the origin,
            with the tips following the transformation along.

- highpts   Same format and functionality as vectors, except here
            no arrow from the origin is shown, only a single
            red marker.

- linesX    linelength-by-nlines matrix, describing the X coordinates
            of nlines to be plotted and to be transfomed with space
            together with the matrix

- linesY    linelength-by-nlines matrix, describing the Y coordinates
            of nlines to be plotted and to be transfomed with space
            together with the matrix

- linescolor the color of these lines

- lineswidth the width of these lines

- userpause  if true, pauses for RETURN on keyboard between
            steps

- pause_steps  only relevant if userpause=true. Pause only before
             the very first pause_steps steps get taken. After
             that run smoothly, without pausing. To pause after
             every step, use pause_steps=Inf

- pausetime   Number of seconds to pause after every step

- first_pause_time  Number of seconds to pause before very first step

- grid   Boolean. If true, constant grey line grids are shown

- anim   If not nothing, assumed to be an Animation() object and
         each frame will be stored in it

- fps    Supposedly, frames per second for .mov animation. However,
         doesn't seem to work, the external movie build call is not
         affected by this parameter.

- srand_seed    An integer, a seed for the random number generator.
        A given seed always produces the same random number sequence, making
        that sequence reproducible.

- plot_kwargs  A Dict(), representing extra parameters to be sent
              to the plot() commands.


# RETURNS

- A Dict() with the following keys:

- "data"   A 2-by-n matrix with the final position of all the n points on the plot

- "tailsX" an ntailpoints-by-n matrix containing the horizontal position of each
    of the n tails. Will be 0-by-n if no tails were plotted

- "tailsY" an ntailpoints-by-n matrix containing the vertical position of each
        of the n tails. Will be 0-by-n if no tails were plotted

- "vectors"  A 2-by-n matrix of black arrow vector positions on the plot. The
        origin of each vector is the (0,0) point.

- "highpts"  A 2-by-n matrix of red marker positions on the plot.

- "linesX"  an nlinepoints-by-nlines matrix containing the final horizontal
    positions of each of the nlines drawn and transformed with M.
    Will be 0-by-0 if no lines were plotted.

- "linesY"  an nlinepoints-by-nlines matrix containing the final vertical
    positions of each of the nlines drawn and transformed with M.
    Will be 0-by-0 if no lines were plotted.


# EXAMPLES
```jldoctest
    animate_matrix([1.62 0.65 ; 0.65 0.88])

    animate_matrix([1.62 0.65 ; 0.65 0.88], plot_tails=false, npoints=300)

    animate_matrix([1.62 0.65 ; 0.65 0.88], filename="my_movie")

    M3 = [-0.125 0.22 ; -0.22 -0.125]
    theta = 60*pi/180; R = [cos(theta) sin(theta) ; -sin(theta) cos(theta)]
    animate_matrix((2*M3, R), srand_seed=123,
        principal_axes=true, plot_tails=false, vectors=[1 0 ; 0 1], highpts=[0.5, 0.5])
    animate_matrix(2*M3*R, srand_seed=123,
        principal_axes=true, plot_tails=false, vectors=[1 0 ; 0 1], highpts=[0.5, 0.5])

    out = animate_matrix(R);
    animate_matrix(M3, seed=out)
```

"""
function animate_matrix(M::Array ; nsteps=30,
    seed=:random, mx=3, npoints=100, principal_axes=false,
    pax_linewidth = 0.5, pax_linecolor = (0.7, 0.7, 0.7),
    plot_tails = true, tailscolor=:green,
    markercolor=:mediumorchid, markersize=5,
    filename="", vectors=zeros(2,0), highpts=zeros(2,0),
    linesX = zeros(0,0), linesY = zeros(0,0), linescolor=:green, lineswidth=1,
    userpause=false, pause_steps=1, grid=true,
    pausetime=0.025, first_pause_time=1.0,
    anim=nothing, fps=10, srand_seed=nothing, plot_kwargs...)

    vectors = copy(vectors); highpts = copy(highpts)

    if seed==:random
        max_radius = mx
        if srand_seed != nothing; Random.seed!(srand_seed); end

        X = (2*rand(npoints).-1)*max_radius
        Y = (2*rand(npoints).-1)*max_radius
        data = [X' ; Y']
        data = data[:,findall(data[1,:].^2 + data[2,:].^2 .< max_radius.^2)]

        tailsX = zeros(0, size(data)[2]) ;
        tailsY = zeros(0, size(data)[2]) ;
    elseif typeof(seed)<:Array
        seed = reshape(seed, 2, Int64(length(seed)/2))
        data = seed;
        tailsX = zeros(0, size(data)[2]) ;
        tailsY = zeros(0, size(data)[2]) ;
    elseif typeof(seed)<:Dict
        data = seed["data"]
        if !isempty(seed["tailsX"])
            tailsX = seed["tailsX"]
            tailsY = seed["tailsY"]
        else
            tailsX = zeros(0, size(data)[2]) ;
            tailsY = zeros(0, size(data)[2]) ;
        end
        if !isempty(seed["vectors"])
            vectors = [vectors seed["vectors"]]
        end
        if !isempty(seed["highpts"])
            highpts = [highpts seed["highpts"]]
        end
        if !isempty(seed["linesX"]) && !isempty(seed["linesY"])
            linesX = seed["linesX"]
            linesY = seed["linesY"]
        end
    else
        error("Sorry, I don't know how to work with this seed")
    end

    function standard_background()
        xlim(-(mx+0.1),mx+0.1); ylim(-(mx+0.1),mx+0.1)
        axis("scaled")
        if principal_axes
            vlines(0, ylim()[1]-1, ylim()[2]+1, color=pax_linecolor, linewidth=pax_linewidth; plot_kwargs...)
            hlines(0, xlim()[1]-1, xlim()[2]+1, color=pax_linecolor, linewidth=pax_linewidth; plot_kwargs...)
        end
        if typeof(mx) != Int64
            mxi = round(Int64, mx)+1
        else
            mxi = mx+1
        end
        if grid
            vlines(-mxi:mxi, ylim()[1]-1, ylim()[2]+1, color=(0.7, 0.7, 0.7), linewidth=0.5; plot_kwargs...)
            hlines(-mxi:mxi, xlim()[1]-1, xlim()[2]+1, color=(0.7, 0.7, 0.7), linewidth=0.5; plot_kwargs...)
        end
        axis([-(mx+0.1), mx+0.1, -(mx+0.1), mx+0.1], "scaled")
        gca().set_xticks([])
        gca().set_yticks([])
        gca().axis("off")
    end

    clf(); standard_background()
    tailsH   = plot([], [], linewidth=1, color=tailscolor; plot_kwargs...)[1]
    linesH   = plot([], [], linewidth=lineswidth, color=linescolor; plot_kwargs...)[1]
    vectorsH = plot_arrow([], [], linewidth=3; plot_kwargs...)
    dataH    = plot([], [], "o", markerfacecolor=markercolor, markersize=markersize,
        markeredgecolor="black"; plot_kwargs...)[1]
    highptsH = plot([], [], "D", markerfacecolor="red", markersize=6; plot_kwargs...)[1]

    function plot_it_all()
        tailsH.set_xdata(vcat(tailsX, NaN*ones(1,size(tailsX,2)))[:]);
        tailsH.set_ydata(vcat(tailsY, NaN*ones(1,size(tailsY,2)))[:]);

        linesH.set_xdata(vcat(real(linesX), NaN*ones(1,size(linesX,2)))[:]);
        linesH.set_ydata(vcat(real(linesY), NaN*ones(1,size(linesY,2)))[:]);
        starts = zeros(size(vectors,2), 2);
        stops  = vectors';
        LX, LY = plot_arrow(real(starts), real(stops), coords_only=true)
        vectorsH.set_xdata(LX); vectorsH.set_ydata(LY)
        dataH.set_xdata(data[1,:]);       dataH.set_ydata(data[2,:])
        highptsH.set_xdata(highpts[1,:]); highptsH.set_ydata(highpts[2,:])
    end

    if !isempty(filename); anim = Animation(); end

    # if typeof(anim)<:Plots.Animation; frame(anim); end
    PyPlot.show()

    E = eigen(M)
    D = E.values
    V = E.vectors
    # D = eigvals(M);
    # V = eigvecs(M);
    D1 = Diagonal(complex(D).^(1/nsteps))
    M1 = V*D1*inv(V)

    if plot_tails
        tailsX = [tailsX ; real(data[1,:])']
        tailsY = [tailsY ; real(data[2,:])']
    end
    for i=1:nsteps+1
        if i>1
            data = M1*data
            vectors = M1*vectors
            highpts = M1*highpts
            nlines = M1*[linesX[:]' ; linesY[:]']
        else
            nlines = [linesX[:]' ; linesY[:]']
        end
        if plot_tails
            tailsX = [tailsX ; real(data[1,:])']
            tailsY = [tailsY ; real(data[2,:])']
        end
        linesX = reshape(nlines[1,:], size(linesX)[1], size(linesX)[2])
        linesY = reshape(nlines[2,:], size(linesY)[1], size(linesY)[2])

        plot_it_all()
        if userpause && i<=pause_steps
            print("<RETURN> for next step")
            readline()
        end
        PyPlot.show; gcf().canvas.flush_events()
        if i==1 && first_pause_time > 0
            pause(first_pause_time)
        elseif pausetime > 0
            pause(pausetime)
        end
    end

    if !isempty(filename)
        mov(anim, filename * ".mov", fps=fps)
    end

    return Dict("data"=>real(data), "tailsX"=>tailsX,
        "tailsY"=>tailsY, "vectors"=>real(vectors),
        "highpts"=>real(highpts),
        "linesX"=>real(linesX), "linesY"=>real(linesY))
end



function animate_matrix(M::Tuple ; seed=:random, sleeptime=1,
    filename="", anim=nothing, fps=10,
    vectors = zeros(2,0), highpts = zeros(2,0), kwargs...)

    if !isempty(filename); anim=Animation(); end
    out = animate_matrix(M[end], seed=seed, anim=anim,
        vectors = vectors, highpts=highpts; kwargs...)

    for i=length(M)-1:-1:1
        sleep(sleeptime)
        out = animate_matrix(M[i], seed=copy(out), anim=anim; kwargs...)
    end

    if !isempty(filename)
        mov(anim, filename * ".mov", fps=fps)
    end
    return out
end
