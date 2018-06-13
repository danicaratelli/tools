function forecast(A,C,Y_in,n,p,h)
    #INPUTS:
            #A: stacked matrix with coefficients
            #C: Constant from VAR
            #Y_in: input data, i.e. Yₜ
            #n: number of variables
            #p: number of lags
            #h: forecast ahead
    #OUTPUT:
            #Y_out: Yₜ₊ₕ    
    
    A_tmp = eye(size(A,1),size(A,2));

    for i=1:h-1;
        A_tmp = A_tmp+(A)^i;
    end

    Aₕ₋₁ = A_tmp;
    return Y_out = (A^h)*Y_in+(Aₕ₋₁*C);
end
