##-----     INTERPOLATIONS FUNCTIONS     -----##
"""
    interpolate1D(xq, xgrid)

    Linear interpolate in 1 dimension.

   Inputs:
   -------
   ``x``        :       Number,     number for placement on the grid\n
   ``xgrid``    :       array (n),  grid on which x is to be placed\n
   -------
   ```
"""
function interpolate1D(x,xgrid)
    nx = length(xgrid);
    if x<xgrid[1]
        xloc_min = 1;
        wx = 1;
    elseif x>xgrid[end]
        xloc_min = nx-1;
        wx = 0;
    else
        xloc_min = findlast(x.>=xgrid);
        wx = 1 - (x - xgrid[xloc_min])/(xgrid[xloc_min + 1] - xgrid[xloc_min]);
    end
    return xloc_min, wx;
end

"""
    interpolate_coord(x, x1, xqi, xqia ,xqpi)

    Linear interpolate:  `xq = xqpi * x[xqi] + (1-xqpi)*x[xqia]`.

    Code converted from Using the Sequence-Space Jacobian to Solve and Estimate Heterogeneous-Agent Models

   Inputs:
   -------
   ``x``        :       array (n), ascending data points\n
   ``xq``       :       array (nq), query points\n
   ``xqi``      :       array (nq), empty (to be filled with indices of lower bracketing gridpoints)\n
   ``xqia``     :       array (nq), empty (to be filled with indices of upper bracketing gridpoints)\n
   ``xqpi``     :       array (nq), empty (to be filled with weights of lower bracketing gridpoints)\n
   -------
   ```
"""
function interpolate_coord(x,xq,xqi,xqia,xqpi)
    #size of arrays
    nxq, nx = size(xq,1), size(x,1);

    #sort and keep track of initial order
    ind_new = sortperm(xq);
    ind_init = sortperm(ind_new);
    xq = xq[ind_new];

    #take care of value below and above minimum
    id = findall((x[1] .<= xq) .& (xq .< x[end]));
    xqi[(xq .< x[1])] .= 1;
    xqpi[(xq .< x[1])] .= 1;
    xqi[(xq .> x[nx])] .= nx;
    xqpi[(xq .> x[nx])] .= 1;

    #interpolation
    xi = 1;
    x_low = x[1];
    x_high = x[2];

    for xqi_cur in id
        xq_cur = xq[xqi_cur];
        while xi < (nx - 1)
           if x_high>=xq_cur
                break
            end
            xi += 1
            x_low = x_high;
            x_high = x[xi + 1];
        end
        xqpi[xqi_cur] = (x_high - xq_cur)/(x_high - x_low);
        xqi[xqi_cur] = xi;
    end

    # revert back to initial order
    xqpi[:] = xqpi[ind_init];
    xqi[:] = xqi[ind_init];

    # Compute index of point above, or same if last on the list
    xqia[:] = xqi[:] .+ 1
    xqia[(xqia .>= nx)] .= nx;
    xqia[(xq .< x[1])] .= xqi[(xq .< x[1])];

    return xqi, xqia, xqpi;
end
