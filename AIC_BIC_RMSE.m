function [AIC BIC RMSE]=AIC_BIC_RMSE(x0,xc,n,k)
% Compute the Akaike Information Criterion (AIC), Bayesian Information Criterion (BIC), and root mean square error (RMSE)
% x0 is the empirical frequency of the joint distribution; xc is the theoretical probability of the joint distribution; n is the sample size; k is the number of copula function parameters
msum=0;
for i=1:n
    msum=msum+(x0(i)-xc(i))^2;                                                                                                                                                                                                                                                                                                
end
RMSE=sqrt(msum/(n-k));
MSE=msum/(n-k);
AIC=n*log(MSE)+2*k;
BIC=n*log(MSE)+k*log(n);
return
