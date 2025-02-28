function irfwold = woldirf(beta, c, p, horizon)
%%
[BigA, N] = companionMatrix(beta,c,p);

irfwold = zeros(N,N,horizon + 1);

for h=1:horizon+1
    
    temp = BigA^(h-1);
    irfwold(:,:,h) = temp(1:N, 1:N);
    
end

%%
end
