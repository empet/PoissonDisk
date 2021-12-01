import StaticArrays: SVector, @SVector

function distsq(A::SVector{2, T}, B::SVector{2, T}) where T <: Real
    return (A[1]-B[1])^2 + (A[2]-B[2])^2
end

function grid_idxs(point:: SVector{2, T}, cellsize::T) where T <: Real 
    """ Returns cell indices, (i,j)=(row, col) for  `point` """
    return ceil(Int, point[2] / cellsize) + 1, ceil(Int, point[1] / cellsize) + 1
end

function ispoint_inrectangle(point::SVector{2, T}, rlims::Tuple{T, T}) where T <: Real
    """ 
    Check whether a point is within the rectangle 
    [0, rlims[1]] x [0, rlims[2]];
    here rlims =(width, height)
    """
    return (0 <= point[1] <= rlims[1]) && (0 <= point[2] <= rlims[2])
end    

function unif_sampling_annulus(center::SVector{2, T}, inner_rad::T, outer_rad::T; 
                               k=30) where T <: Real
    """ 
    Samples  uniformly k points from the annulus of center, `center`, 
    and radii, inner_radius, outer_radius
    """
    rad = sqrt.(inner_rad^2  .+ (outer_rad^2 - inner_rad^2)*rand(k)) 
    angle = 2*pi * rand(k)  
    x, y = center[1] .+ rad .* cos.(angle), center[2] .+ rad .*sin.(angle) 
    return [@SVector([xi, yi]) for (xi, yi) in zip(x,y)]
end

function Poisson_disk_sampling(;width =1.0, height=1.0, radius=0.05, k=30)
    """ 
    Bridson Algorithm to generate  Poisson disk samples 
    within a 2D rectangle, [0, width] x [0, height]
    Ref: R Bridson, Fast Poisson disk sampling in arbitrary dimensions,
    in  SIGGRAPH '07: ACM SIGGRAPH 2007 sketches
    https://www.cs.ubc.ca/~rbridson/docs/bridson-siggraph07-poissondisk.pdf
    """
    function  is_accepted(point:: SVector{2, T})  where T <: Real
        """ Check whether the point, `point` can be added to the active_list"""

        i, j = grid_idxs(point, cellsize)
        m, n = maximum([i-2, 1]), maximum([j-2, 1])
        p, q = minimum([i+2, grid_rows]), minimum([j+2, grid_cols])
        for k in n:q
            for l in m:p
                index = (k-1) * grid_rows +l
                gpoint = grid[index]
                if ismissing(gpoint)
                    continue
                end
                if distsq(point, gpoint) <= r2
                    return false
                end
            end
        end              
        return true
    end
    
    cellsize = radius/sqrt(2)
    grid_cols = ceil(Int, width / cellsize) + 1
    grid_rows = ceil(Int, height / cellsize) + 1
  
    grid = Vector{Any}(missing, grid_rows*grid_cols)
    r2 = radius*radius
    active_list = SVector{2, Float64}[]
    
    #generate a random point within the rectangle of interest
    # and  insert it into the active_list and grid
    point = @SVector([width*rand(1)[1], height*rand(1)[1]])
    push!(active_list, point)
    i, j = grid_idxs(point, cellsize)
    grid[(j-1)*grid_rows+i] = point
    
    while length(active_list) > 0
        i = rand(1:length(active_list), 1)[1]
        pt = active_list[i] 
        active_list  =  active_list[1:end .!= i] #removes the point in pos i
        rpoints = unif_sampling_annulus(pt, radius, 2*radius; k=k)
        for rp in rpoints
            if !ispoint_inrectangle(rp, (width, height)) || !is_accepted(rp)
                continue
            end    
            push!(active_list, rp)    
            i, j = grid_idxs(rp, cellsize)
            grid[(j-1)*grid_rows+i] = rp
        end
    end
    return filter(x->!ismissing(x), grid);
end    
