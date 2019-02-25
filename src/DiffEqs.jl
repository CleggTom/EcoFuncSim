function lotka_volterra(dc,c,p::Parameters,t)
    for i in 1:p.n_sp
        dc[i] = c[i] * (p.Com.U[i] - p.Com.R[i])
        for j in 1:p.n_sp
            dc[i] -= c[i] * c[j] * p.Com.a[i,j]
        end
    end
end

function mean_lotka_volterra(dc,c,p::Parameters,t)
    for i in 1:p.n_sp
        dc[i] = (c[i] * (p.Com.U[i] - p.Com.R[i])) -
                (c[i] * c[i] * p.Com.a[i,i]) -
                (c[i] * (sum(p.Com.a[i,:])-p.Com.a[i,i]) * (sum(c)/p.n_sp) )
    end
end
