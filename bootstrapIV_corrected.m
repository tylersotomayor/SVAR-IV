function [upper, lower, meanirf, medianirf] = bootstrapIV_corrected(y,p,c,beta,residuals,Z, nboot1, nboot2, adjustZ, adjustu, policyvar, horizon, prc)

[T, N] = size(y);
var_estimates = zeros(N, size(beta,1), nboot1);

for b1=1:nboot1
    
    [varboot] = bootstrapVAR(y,p,c,beta,residuals,'residual');
    [betaloop, ~] = VAR(varboot, p, 1);
    var_estimates(:,:,b1) = betaloop';
    
end
    
% Now remove the bias
[beta_new, ~] = remove_bias(beta,c,p,var_estimates);

% Save the bootstrapped IRFs here
ivirf_boot = zeros(N,horizon+1, nboot2);

for b2=1:nboot2
    
    % Bootstrap the VAR
    [varboot,rademacher] = bootstrapVAR(y,p,c,beta_new',residuals,'wild');
    
    % Compute residuals and VAR coefficients
    [betaloop, residuals_loop] = VAR(varboot, p, 1);

    % Adjust the sizes as indicated to have same time periods
    u_p_loop = residuals_loop(:,policyvar);
    u_p_loop_final = u_p_loop(adjustu(1):adjustu(2));
    Z_final = Z(adjustZ(1):adjustZ(2));
    u_q_final = residuals_loop(adjustu(1):adjustu(2),:);
    u_q_final(:,policyvar) = [];
    
    % Bootstrap the instrument
    Znew = Z_final.*rademacher(adjustu(1):adjustu(2));

    % Compute Wold IRFs
    wold_loop = woldirf(betaloop,c,p,horizon);
    
    % Compute first-stage LS
    [~, uhat, ~] = myols(u_p_loop_final, Znew,0);
    
    % Compute second-stage LS
    sq_sp = (uhat'*uhat)\uhat'*u_q_final;

    % We normalize sp = 1 
    s = zeros(N,1);
    looper = 1:1:N;
    looper(policyvar) = [];
    s(looper) = sq_sp;
    s(policyvar) = 1;

    % Now we can compute the structural IRFs as
    ivirf = zeros(N, horizon+1);
    for h=1:horizon+1
        ivirf(:,h) = wold_loop(:,:,h)*s;
    end

  ivirf_boot(:,:,b2) = ivirf; 

    
end

up = (50 + prc/2);
low = (50 - prc/2);

upper = prctile(ivirf_boot,up,3);
lower = prctile(ivirf_boot,low,3);
medianirf = prctile(ivirf_boot, 50,3);
meanirf = mean(ivirf_boot,3);


end