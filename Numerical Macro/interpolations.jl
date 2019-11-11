function simple_interpolation(xs,ys,x)
    if x ∈ xs
        return ys[findfirst(xs.==x)];
    else
        x_pre = findfirst(xs.<=x);
        x_post = findfirst(xs.>x);
        y = ys[x_pre] + ((x - xs[x_pre])/(xs[x_post] - xs[x_pre]))*(ys[x_post] - ys[x_pre]);
        return y;
    end
end
