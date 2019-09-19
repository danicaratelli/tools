using LinearAlgebra

function NewtonRaphson(f, fprime, xstart, Xs, tol = 1e-8, maxiter = 10000)
    fnew(x) = x - f(x) / fprime(x)

    dd = Inf
    it = 0
    xloop = copy(xstart)
    @assert ((xloop >= Xs[1]) && (xloop <= Xs[end])) "Starting value is not in function support Xs"
    while (abs(f(xloop)) > tol) && (it < maxiter)
        xloop = fnew(xloop)
        if xloop < Xs[1]
            xloop = Xs[1]
        elseif xloop > Xs[end]
            xloop = Xs[end]
        end
        it += 1
    end
    return xloop
end

"""
    NewtonRaphson_secant(f,  x0, x1, Xboundry, tolerance, maxiter)

    This function computes the zeros for univariate function f given starting
        value x0 and x1 using the secant method.
        It looks for a solution only within Xboundry.
        The algorithm terminates if f(x) is smaller than indicated by tolerance
        or if the number of iterations go beyond maxiter.
        For more details see Ken Judd (p. 158) or Heer & Maussner (p. 609-610)
"""

function NewtonRaphson_secant(f, x0, x1, Xs, tol = 1e-8, maxiter = 10000)

    dd = Inf
    it = 0
    @assert ((x0 >= Xs[1]) &&
             (x0 <= Xs[end]) && (x1 >= Xboundry[1]) && (x1 <= Xboundry[end])) "Starting value is not in function support Xs"
    while (abs(f(x1)) > tol) && (it < maxiter)
        xnew = x1 - f(x1) * (x1 - x0) / (f(x1) - f(x0))
        if xnew < Xboundry[1]
            xnew = Xboundry[1]
        elseif xnew > Xboundry[end]
            xnew = Xboundry[end]
        end
        x0 = copy(x1)
        x1 = copy(xnew)
        it += 1
    end
    return x1
end

"""
    NewtonRaphson_MV(f, J, x0, Xboundry, tolerance, maxiter)

    This function computes the zeros for function f given the true analytic
        Jacobian J from a starting value x0.
        It looks for a solution only within Xboundry.
        The algorithm terminates if f(x) is smaller than indicated by tolerance
        or if the number of iterations go beyond maxiter.
        For more details see Ken Judd (p. 168) or Heer & Maussner (p. 612)
"""
function NewtonRaphson_MV(f, J, x0, Xs, tol = 1e-8, maxiter = 10000)
    fnew(x) = x - inv(J(x)) * f(x)

    dd = Inf
    it = 0
    ndim = length(x0)
    xloop = copy(x0)
    #making sure starting x is in support
    for n = 1:ndim
        @assert ((xloop[n] >= Xboundry[n, 1]) && (xloop[n] <= Xboundry[n, end])) "Starting value is not in function support Xs"
    end
    while (sum(abs.(f(xloop))) > tol) && (it < maxiter)
        xloop = fnew(xloop)
        for n = 1:ndim
            if xloop[n] < Xboundry[n, 1]
                xloop[n] = Xboundry[n, 1]
            elseif xloop[n] > Xboundry[n, end]
                xloop[n] = Xboundry[n, end]
            end
        end
        it += 1
    end
    return xloop
end


"""
    Broyden_secant(f, J, x0, Xboundry, tolerance, maxiter)

    This function computes the zeros for function f given an initial guess of the
        Jacobian J (usually initialized to the identity matrix) from a starting
        value x0. It looks for a solution only within Xboundry.
        The algorithm terminates if f(x) is smaller than indicated by tolerance
        or if the number of iterations go beyond maxiter.
        For more details see Ken Judd (p. 170) or Heer & Maussner (p. 614)
"""
function Broyden_secant(f, J, x0, Xs, tol = 1e-8, maxiter = 10000)
    Jnew = J;

    it = 0;
    ndim = length(x0)
    xold = copy(x0)
    #making sure starting x is in support
    for n = 1:ndim
        @assert ((xold[n] >= Xboundry[n, 1]) && (xold[n] <= Xboundry[n, end])) "Starting value is not in function support Xs"
    end
    while (sum(abs.(f(xold))) > tol) && (it < maxiter)
        xnew = xold - inv(Jnew)*f(xold);
        for n = 1:ndim
            if xnew[n] < Xboundry[n, 1]
                xnew[n] = Xboundry[n, 1]
            elseif xnew[n] > Xboundry[n, end]
                xnew[n] = Xboundry[n, end]
            end
        end
        w = xnew-xold;
        Jnew = Jnew + (((f(xnew)-f(xold))-Jnew*w)*w')/(w'*w);
        xold = copy(xnew);
        it += 1
    end
    return xold
end
