function fval = WeighedAverageSelling(r, D, P, var_ub)
    fval = sum(r/sum(P)*P.^2);
end