function [Fstat, pval] = do_ftest(Ru, Rr, m, df)

    % Ru = R squared statistic from unrestricted model
    % Rr = R squared from restricted model
    % m = number of restrictions (e.g. if all parameters are zero it is K)
    % df = degrees of freedom (sample size used in OLS estimation)

    Fstat = ((Ru - Rr) / (m)) / ((1 - Ru) / (df));
    pval       = 1 - fcdf(Fstat,m,df);


end