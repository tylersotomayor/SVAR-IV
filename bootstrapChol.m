function [bootchol, upper, lower, boot_beta] = bootstrapChol(y,p,c,beta,residuals,nboot,horizon,prc,bootscheme)

[T, N] = size(y);
bootchol = zeros(N,N,horizon+1,nboot);
boot_beta = zeros(N, size(beta,1), nboot);

for b=1:nboot
    
    varboot = bootstrapVAR(y,p,c,beta,residuals,bootscheme);
    [betaloop, err_loop] = VAR(varboot, p, 1);
    boot_beta(:,:,b) = betaloop';
    wold_loop = woldirf(betaloop,c,p,horizon);
    sigma_loop = (err_loop' * err_loop) ./ (T - 1 - p - N*p);
    S_loop = chol(sigma_loop, 'lower');
    cholirf_loop = choleskyIRF(wold_loop, S_loop);
    bootchol(:,:,:,b) = cholirf_loop;
    
end

up = (50 + prc/2);
low = (50 - prc/2);

upper = prctile(bootchol,up,4);
lower = prctile(bootchol,low,4);


end