function Q = QM(O, M, xmpt)
    f_obs = fitter(reshape(O, [], 1), 'verbosity', 0, 'sortby', 'AIC');
    f_mod = fitter(reshape(M, [], 1), 'verbosity', 0, 'sortby', 'AIC');
    epst  = cdf(f_mod, xmpt(:));
    epst  = max(epst, 1e-6);
    epst  = min(epst, 1-1e-6);
    Q     = icdf(f_obs, epst);
    Q(isinf(Q) | isnan(Q)) = NaN;
end