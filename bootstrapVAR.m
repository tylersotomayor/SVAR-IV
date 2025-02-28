function [ynext, rademacher] = bootstrapVAR(y,p,c,beta,residuals,bootscheme)


[T, N] = size(y);
yinit = y(1:p, :);

if c==1
    const = beta(1,:);
    pi = beta(2:end,:);
else
    const = 0;  
    pi = beta;
end

ynext=zeros(T,N);
ynext(1:p,:) = yinit;
yinit=reshape(flipud(yinit)',1,[]);

if isequal(bootscheme,'residual')
    for i=1:T-p 
      ynext(p+i,:)= const + yinit*pi + residuals(randi(T-p),:);  
      yinit = reshape(flipud(ynext(i+1:p+i,:))',1,[]);  
    end
    
elseif isequal(bootscheme,'wild')
    rademacher = 1-2*(rand(T-p,1)>0.5); % Draw a uniform random number between 0 and 1, check if larger than 0.5 (equal probs)
                                        % If yes, equal 1, if false equal 0. 
                                        % 1 - 2 * this indicator gives 1 or -1.
    residuals_wild = residuals .* rademacher;                                   
                                        
    for i=1:T-p 
      ynext(p+i,:)= const + yinit*pi + residuals_wild(i,:);  
      yinit = reshape(flipud(ynext(i+1:p+i,:))',1,[]);  
    end

    
end


end