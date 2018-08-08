function irf(mod::model,h::Int64,λ,D,v_type,b=0,K=200,bw=68)
    M = copy(mod);
    T = M.T;
    n = M.n;
    p = M.lags;
    data = M.data;

    k=1;
    d=MvNormal(zeros(n),M.variance);
    if b!=0
        IRFs = zeros(h+1,n,n,K);
        YY = zeros(T,n,K);
        while k<=K
            Y_out = zeros(T,n);
            Y_out[1:p,:] = data[1:p,:];
            M_t = copy(M);
            for t=p+1:T
                E = rand(d,1);
                Y_in = Y_out[1:t-1,:];
                M_t.Y = Y_in;
                fore_t = forecast(M_t,1);
                Y_out[t,:] = fore_t.+E;
            end

            YY[:,:,k] = Y_out;
            if v_type == "bvar"
                M_tmp = bvar_m(Y_out,p,λ,D[:]);
            else
                M_tmp = var_m(Y_out,p);
            end
            IRFs[:,:,:,k] = irf_help(M_tmp,h);
            k = k+1;
        end
        #quantiles
        qnt_high = zeros(h+1,n,n);
        qnt_low = zeros(h+1,n,n);
        qnt_med = zeros(h+1,n,n);
        for ii = 1:n
            for jj = 1:n
                set_sort = sort(IRFs[:,ii,jj,:],2);
                qnt_high[:,ii,jj] = set_sort[:,Int( K*(50+(bw/2))/100)];
                qnt_low[:,ii,jj] = set_sort[:,Int( K*(50-(bw/2))/100)];
                qnt_med[:,ii,jj] = set_sort[:,Int( K*(50/100))];
            end
        end
    else
        qnt_high = zeros(h+1,n,n);
        mod_new = var_m(data[p+1:end,:],p);
        qnt_med = medians = irf_help(mod_new,h);
        qnt_low = zeros(h+1,n,n);
    end

    return qnt_high, qnt_med, qnt_low;
end

function irf_help(M::model,h::Int64)
    T = M.T;
    n = M.n;
    p = M.lags;

    B = M.coeff; #model coefficients
    A = zeros(n*p,n*p);
    A[1:n,:] = B[1:end-1,:]';
    A[n+1:end,1:n*(p-1)] = eye(n*(p-1));

    CC = zeros(n*p,n);
    Su = M.variance;
    CC[1:n,1:n] = chol(Su)';
    JJ = [eye(n) zeros(n,n*(p-1))]; #indicator matrix to select correct entries
    Aj = eye(n*p);

    IRF = zeros(h+1,n,n);

    for j=0:h
        IRF[j+1,:,:] = JJ*Aj*CC*diagm(1./diag(CC)); #last multiplicatioon is to normalize shock to 1 standard deviation
        Aj = Aj*A;
    end
    return IRF;
end
