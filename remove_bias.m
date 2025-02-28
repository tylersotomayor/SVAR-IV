function [Beta, corrections] = remove_bias(beta,c,p,boot_beta)

% Function to compute and correct for small sample bias in bootstrap

% Inputs:   beta        = (Np+1 x N) matrix of estimated coefficients (Np x N) if
%                          no constant is included 
%           p           = integer VAR lag order
%           c           = 1 if constant required
%           boot_beta   = (Np+1 x N x nboot) array of bootstrapped VAR
%           coefficients

% Outputs:  Beta        = (N x Np+1) matrix of bias corrected VAR
%                          parameters
%           corrections = integer number of total corrections needed 

    % Compute the bias of the parameter estimates as the bootstrap mean
    % estimates minus the point estimate
    bias = mean(boot_beta,3) - beta';

    % Check if the estimates imply an unstable run
    BigA = companionMatrix(beta,c,p);
    ev = abs(eig(BigA));

    if max(ev) >=1
        Beta = beta';
    else
        Beta = beta' - bias;
    end

    % If the original estimates are unstable, correct iteratively towards
    % stationarity

    corrections = 1;
    delta_loop = 1;
    bias_loop = bias;
    beta_loop = Beta;

    while max(ev) >= 1

        Beta = beta_loop - bias_loop .* delta_loop;
        BigA_loop = companionMatrix(Beta',c,p);
        ev = abs(eig(BigA_loop));

        delta_loop = delta_loop - 0.01;
        corrections = corrections + 1;

    end

end