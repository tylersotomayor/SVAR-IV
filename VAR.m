function [beta, residuals]  = VAR(y,p,c)

[T, ~] = size(y);
yfinal = y(p+1:T,:);
if c == 1
    X = [ones(T-p,1), lagmakerMatrix(y,p)];
else
    X = lagmakerMatrix(y,p);
end

beta = (X'*X)\X'*yfinal;
residuals = yfinal - X*beta;

end