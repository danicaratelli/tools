function NewtonRaphson(f,fprime,xstart,Xs,tol=1e-8,maxiter=10000)
    fnew(x) = x - f(x)/fprime(x);

    dd = Inf;
    it = 0;
    xloop = copy(xstart);
    @assert ((xloop>=Xs[1]) && (xloop<=Xs[end])) "Starting value is not in function support Xs"
    while (abs(f(xloop))>tol) && (it<maxiter)
        xloop = fnew(xloop);
        if xloop<Xs[1]
            xloop = Xs[1];
        elseif xloop>Xs[end]
            xloop = Xs[end];
        end
        it+=1;
    end
    return xloop;
end


function NewtonRaphson_secant(f,x0,x1,Xs,tol=1e-8,maxiter=10000)

    dd = Inf;
    it = 0;
    @assert ((x0>=Xs[1]) && (x0<=Xs[end]) && (x1>=Xs[1]) && (x1<=Xs[end])) "Starting value is not in function support Xs"
    while (abs(f(x1))>tol) && (it<maxiter)
        xnew = x1  - f(x1)*(x1 - x0)/(f(x1)-f(x0));
        if xnew<Xs[1]
            xnew = Xs[1];
        elseif xnew>Xs[end]
            xnew = Xs[end];
        end
        x0 = copy(x1);
        x1 = copy(xnew);
        it+=1;
    end
    println(it)
    return x1;
end


function NewtonRaphson_MV(f,J,xstart,Xs,tol=1e-8,maxiter=10000)
    fnew(x) = x - inv(J(x))*f(x);

    dd = Inf;
    it = 0;
    ndim = length(xstart);
    xloop = copy(xstart);
    #making sure starting x is in support
    for n=1:ndim
        @assert ((xloop[n]>=Xs[n,1]) && (xloop[n]<=Xs[n,end])) "Starting value is not in function support Xs"
    end
    while (sum(abs.(f(xloop)))>tol) && (it<maxiter)
        xloop = fnew(xloop);
        for n=1:ndim
            if xloop[n]<Xs[n,1]
                xloop[n] = Xs[n,1];
            elseif xloop[n]>Xs[n,end]
                xloop[n] = Xs[n,end];
            end
        end
        it+=1;
    end
    println(it)
    return xloop;
end
