function [comp, N] = companionMatrix(beta, c, p)

if c==1
    N = (size(beta,1) - 1) / p;
    beta = beta(2:end,:);
else
    N = size(beta,1) / p;
end

comp = zeros(N*p,N*p);

comp(1:N,:) = beta';
comp(N+1:end,1:N*(p-1)) = eye(N*(p-1));


end