function [chol] = choleskyIRF(wold, S)

[N,~, horizon] = size(wold);
chol = zeros(N,N,horizon);

for h=1:horizon
    
    chol(:,:,h) = wold(:,:,h) * S;
    
end

end