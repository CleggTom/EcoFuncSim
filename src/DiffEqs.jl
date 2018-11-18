function Lotka(dc,c,p::Parameters,t)

    for i in 1:p.n_sp
        dc[i] = c[i] * (p.Com.U[i] - p.Com.R[i])
        for j in 1:p.n_sp
            dc[i] -= c[i] * c[j] * p.Com.a[i,j]
        end

        # if (c[i] + dc[i]) < 0
        #     c[i] = 0.0
        # end
    end
end
