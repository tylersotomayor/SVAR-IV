function [upper_corrected, lower_corrected, corrections] = bootstrapChol_corrected(beta,boot_beta,p,c,y,residuals,nboot,horizon,prc,bootscheme)

% This code uses the shortcut proposed by Kilian that the bias term is
% hardly affected by the second bootstrap so that we can save one loop 

if c==1
    N = (size(beta,1) - 1) / p;
else
    N = size(beta,1) / p;
end

% Compute the bias of the parameter estimates as the bootstrap mean
% estimates minus the point estimate
bias = mean(boot_beta,3) - beta';

% Check if the point estimates imply an unstable run
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
    BigA_loop = companionMatrix(Beta,c,p);
    ev = abs(eig(BigA_loop));

    delta_loop = delta_loop - 0.01;
    corrections = corrections + 1;
    
end

% Bootstrap again, now using the bias adjusted parameter estimates
[~, upper_corrected, lower_corrected, ~] = bootstrapChol(y,p,c,Beta',residuals,nboot,horizon,prc, bootscheme);


end
