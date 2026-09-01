function C = QDM(O, M, xmpt)
    f_obs  = fitter(reshape(O, [], 1), 'verbosity', 0, 'sortby', 'AIC');
    f_mod  = fitter(reshape(M, [], 1), 'verbosity', 0, 'sortby', 'AIC');

    n      = length(xmpt);
    C      = zeros(n, 1);
    window = 120; % ±10 tahun = ±120 bulan

    for t = 1:n
        t_start = max(1, t - window);
        t_end   = min(n, t + window);
        xmpt_win = xmpt(t_start:t_end);

        try
            f_xmpt = fitter(reshape(xmpt_win, [], 1), 'verbosity', 0, 'sortby', 'AIC');
            epst = cdf(f_xmpt, xmpt(t));
            epst = max(epst, 1e-6);
            epst = min(epst, 1-1e-6);
            delt  = xmpt(t) - icdf(f_mod, epst);
            xhatt = icdf(f_obs, epst);
            C(t)  = xhatt + delt;
            if isinf(C(t)) || isnan(C(t))
                C(t) = NaN;
            end
        catch
            C(t) = NaN;
        end
    end
end