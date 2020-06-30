using Polynomials, LinearAlgebra

# Defining functions to map between [a,b]->[-1,1] (X) & between [-1,1]->[a,b] (Z)
X(z,a,b) = (2z)./(b-a) .- (a + b)/(b-a);
Z(x,a,b) = (x .+ 1).*(b-a)/2 .+ a;


## Chebyshev approximation
"""
    Cheby_approx(n::Int64,m::Int64,a::Number,b::Number)
    Function outputs the coefficients for the Chebyshev polynomial approximation.
    The inputs are:
   -------
   ``n``   :       Int64,       order of Chebyshev polynomial\n
   ``m``   :       Int64,       Chebyshev nodes at which to evaluate polynomial\n
   ``a``   :       Number,      lower bound of support of original function\n
   ``b``   :       Number,      upper bound of support of original function\n
   -------
    The output is:
   -------
   ``Cheby_coffs``:   Array{Float64,n+1},     (n+1) coefficients for approx.\n
   -------
   ```
"""
function Cheby_approx(n::Int64,m::Int64,a::Number,b::Number)
    #Step 0: initializing vector of coefficients
    Cheby_coffs = zeros(n+1);
    # Step 1:
    xs = sort(map(x->cos((2x-1)*(π/(2*m))),1:m));  #compute Chebyshev nodes
    zs = Z(xs,a,b);     #corresponding nodes in [a,b] interval
    # Step 2:
    Fs = f.(zs);    #evaluating function at original points
    # Step 3: computing coefficients
    Cheby_coffs[1] = (1/m)*sum(Fs);
    Ts = Cheby_T_vec(i,xs)
    Cheby_coffs[2:end] = map(i->(2/m)*sum(Fs.*Cheby_T_vec(i,xs)),1:n);
    #returning coefficients
    return Cheby_coffs;
end


#evaluation (vecotrize and sum)
"""
    Cheby_T(n::Int64,x::Number)
    Function evaluted Chebyshev basis element T_n at x:
   -------
   ``n``   :       Int64,       degree of Chebyshev polynomial\n
   ``x``   :       Number,      argument where polynomial is evaluated\n
   -------
   ```
   The output is:
  -------
  ``pCheb(X)``:   Number,       T_n(x).\n
  -------
  ```
"""
function Cheby_T(n::Int64,x::Number)
    degvec = zeros(n+1); degvec[n+1] = 1;
    pCheb = ChebyshevT(degvec);
    return pCheb(x);
end

"""
    Cheby_T_vec(n::Int64,x::Array{Float64,1})
    Function evaluted Chebyshev basis element T_n at x:
   -------
   ``n``   :       Int64,               degree of Chebyshev polynomial\n
   ``x``   :       Array{Float64,1},    vector where polynomial is evaluated\n
   -------
   ```
   The output is:
  -------
  ``pCheb.(X)``:   Number,       T_n.(x).\n
  -------
  ```
"""
function Cheby_T_vec(n::Int64,x::Array{Float64,1})
    degvec = zeros(n+1); degvec[n+1] = 1;
    pCheb = ChebyshevT(degvec);
    return pCheb.(x);
end


## Example: approximating max(0,x-1) using n=5 and m=30 over [0,2]
f(x) = max(0,x-1);
Cheby_coffs = Cheby_approx(5,30,a,b);
f_approx = ChebyshevT(Cheby_coffs);

using PyPlot; pygui(true);
plot(collect(range(a,stop=b,length=100)),f.(collect(range(a,stop=b,length=100))))
plot(collect(range(a,stop=b,length=100)),f_approx.(collect(range(-1,stop=1,length=100))))

plot(Z(collect(range(-1,stop=1,length=100))),res)
