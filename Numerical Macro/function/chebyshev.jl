"""
Computes n many Chebyshev nodes
"""
function chebynodes(n)
    ns = map(x->cos((2x-1)*π/2n),1:n);
    return sort(ns);
end


"""
Evaluates a Chebyshev polynomial of degree k at each entry of array y
"""
function Tcheb(y,k)
    Tj0 = ones(length(y));
    Tj1 = y;
    if k==0
        res =  Tj0;
    elseif k==1
        res = Tj1;
    else
        for j=2:k
            Tj = 2*y.*Tj1 .- Tj0;
            Tj0 = copy(Tj1);
            Tj1 = copy(Tj);
        end
        res = Tj1;
    end
    if length(res)==1
        return res[1];
    else
        return res;
    end
end


"""
Forms all possible combinations of elements in the grids contained in xs
	Input: 	xs: 	tuple of grids to be combines
"""
function combine_vecs(xs)
    nvars = length(xs);                     #number of vectors in question
    nxs = map(x->length(xs[x]),1:nvars);    #length of each vector
    vec_locs = zeros(prod(nxs),nvars);
    for j=1:nvars
        #num_rep = nxs[j].^(nvars .- collect(1:nvars)); #number of times to repeat block
        if j==nvars
            num_rep = 1;
        else
            num_rep = prod(nxs[j+1:end]);
        end
        block = repeat((xs[j][1:nxs[j]])',prod(nxs[1:j-1]),1);  #block to repeat
        vec_locs[:,j] = repeat(block[:],num_rep);
        #Int(prod(nxs)/length(block)));
    end
    return vec_locs;
end



"""
Computes Chebyshev coefficients (assuming collocation) for any number of state variables\n
	Input:
    -------
    `ns`: 		tuple for each dimension of the state space\n
	``Ys``: 		Function of interest evaluated at fnodes\n
	``fnodes``: 	array containing the nodes (Chebyshev nodes) at which Ys is evaluated\n
    -------
	Output:
    -------
    ``coffs``:		coefficients of Chebyshev approximation\n
	``vec_locs``:	combination of all grid point elements in 2D array\n
    -------
"""
function chebycoeffs(ns,Ys,fnodes)
    ndim = length(ns);
    nvert = prod(ns);
    k = ns[1];
    @assert prod(ns.==k) "This program only works for collocation algorithms"
    coffs = zeros(ns);

    #vectorizing coefficients and data
    coffs = coffs[:];
    Ys_vert = Ys[:];
    #identifying for entire ndim grid where 1st coordinate is for each dimension
    #e.g. [a b;c d] --> the zeros in dimension 1 are (1,1) and (1,2) while in
    # dimension 2 they are (1,1) and (2,1)
    tup = tuple();
    for j=1:ndim
        tup = tuple(tup...,collect(1:k));
    end
    vec_locs = Int.(combine_vecs(tup));

    coffs[1] = mean(Ys);
    #including non-Cheby related term
    to_rem = (vec_locs.==1);
    coffs[2:end] = ((1/k)^ndim).*(2 .^(ndim.-sum(to_rem[2:end,:],dims=2)));

    for j=2:nvert
        to_cheb = setdiff(1:ndim,findall(to_rem[j,:]));    #those dimensions for which
        #we need to evaluate Chebyshev polynomial
        cheby_prod = ones(nvert,1);
        for jj in to_cheb
            tcheb_jj = Tcheb(fnodes[jj],vec_locs[j,jj]-1);
            cheby_prod = cheby_prod.*tcheb_jj[vec_locs[:,jj]];
        end
        coffs[j] = coffs[j].*sum(Ys_vert.*cheby_prod);
    end
    coffs = reshape(coffs,ns);
    return coffs, vec_locs;
end



"""
Evaluates a Chebyshev approximation given coefficients for the polynomial
	Input:\n
    -----
    coffs: 	coefficients of polynomial (see 1st output of chebycoeffs)\n
	T_deg: 	combination of all grid point elements in 2D array (see 2nd output of chebycoeffs)\n
	xs: 	tuple of arrays at which polynomial is evaluated\n
    -----
    # Example\n
    ```
    julia> f(x,y) = sin(x) + cos(y)
    julia> ns = (5,5); fnodes = (chebynodes(5),chebynodes(5));
    julia> Ys = sin.(Z(fnodes[1],0,1)) .+ (cos.(Z(fnodes[1],0,1)))';
    julia> coffs,vec_locs = chebycoeffs(ns,Ys,fnodes);
    julia> xs = (collect(-1:0.01:1),collect(-1:0.01:1));
    julia> F_hat = Cheby_eval(coffs,vec_locs,xs)
    ```
"""
function Cheby_eval(coffs,T_deg,xs)
    nvars = size(xs,1);
    nnodes = map(x->length(xs[x]),1:nvars);
    ns = size(coffs);
    k = ns[1];
    @assert prod(ns.==k) "This program only works for collocation algorithms"
    #@assert prod(nnodes[2:end].==1) "Can only take the first dimension with multiple entries"
    nvert = prod(ns);
    T_val = zeros(nvert,nvars,prod(nnodes)); #evaluated Chebyshev

    X_val = combine_vecs(xs);

    for j=1:k
        for jj=1:nvars
            to_place = findall(T_deg[:,jj].==j);
            if length(X_val[:,jj])>1
                T_val[to_place,jj,:] .= repeat(Tcheb(X_val[:,jj],j-1)',length(to_place));
            else
                T_val[to_place,jj,:] .= repeat([Tcheb(X_val[:,jj],j-1)],length(to_place));
            end
        end
    end
    r = sum(repeat(coffs[:],1,prod(nnodes)).*dropdims(prod(T_val,dims=2),dims=2),dims=1)';
    r = reshape(r,nnodes...);
    r = dropdims(r,dims=tuple(findall(size(r) .== 1)...));
    if prod(nnodes)==1
        return r[1];
    else
        return r;
    end
end
