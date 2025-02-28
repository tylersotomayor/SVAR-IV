function [beta, fitted, err, r2, fitted_partial] = myols(y, X, c)

[T,N] = size(X);

if c==1
    X = [ones(T,1) X];
else
    X = X;
end

beta = (X'*X)\X'*y;

fitted = X*beta;

if c == 1
    fitted_partial = X(:,2:end)*beta(2:end);
else
    fitted_partial = fitted;
end

err = y - fitted;

rss = err'*err;
tss = sum((y - mean(y)).^2);
r2 = 1 - rss/tss;

end