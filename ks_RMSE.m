function [h,p,RMSE]=ks_RMSE(data,distribution,scale)
% Using the Kolmogorov–Smirnov (K–S) test and root mean square error (RMSE) to test the goodness of fit for the marginal distributions.
% Data is the variable data that needs to be tested for goodness of fit; distribution is the type of functional distribution; scale is the time scale.
A1=[];
for is=1:scale
    A1=[A1,data(is:length(data)-scale+is)];
end
    XS=sum(A1,2);
nseas=12;
for is=1:nseas
    tind=is:nseas:length(XS);
    Xn=XS(tind);
    n=length(Xn);
    [p1,X1]=ecdf(Xn);
    if isequal(distribution,'exp') % Exponential (Exp) Function
        [a1,a2]=expPWM(Xn);   
        u=expcdf(X1,a1);  
    elseif isequal(distribution,'wbl') % Weibull (Wbl) Function
        wblp=wblfit(Xn); b1=wblp(1);b2=wblp(2);   
        u=wblcdf(X1,b1,b2);   
    elseif isequal(distribution,'gamma') %  Gamma (Gam) Function
        gamp=gamfit(Xn);A=gamp(1);B=gamp(2);   
        u=gamcdf(X1,A,B);      
    elseif isequal(distribution,'p3') % Pearson III (P-III) Function
        xa=mean(Xn);xst=std(Xn);
       cs=skewness(Xn,0);
        A=(2/cs)^2;
        B=sqrt(xst^2/A);
        a0=xa-A*B;
        X=X1-a0;
        u=gamcdf(X,A,B);     
    elseif isequal(distribution,'logn') %  Log-normal (Logn) Function
        lognp=lognfit(Xn);A=lognp(1);B=lognp(2); 
        u=logncdf(X1,A,B);      
    elseif isequal(distribution,'gp') %  Generalized Pareto (Gp) Function
        gpp=gpfit(Xn); c1=gpp(1);c2=gpp(2);  
        u=gpcdf(X1,c1,c2);    
    elseif isequal(distribution,'gev') % Generalized Extreme value (Gev) Function
        gevp=gevfit(Xn);        
        K=gevp(1);sigma=gevp(2);mu=gevp(3);
        u=gevcdf(X1,K,sigma,mu);
    elseif isequal(distribution,'Lol')  % Log-Logistic (Lol) Function
        u=lolcdf(Xn,X1);
     elseif isequal(distribution,'norm') % Normal (Norm) Function
        [A,B]=normfit(Xn);           
        u=normcdf(X1,A,B);    
   end
    [h(is),p(is)]= kstest(X1,'CDF',[X1,u]);
    msum=0;
    for i=1:n
        msum=msum+(p1(i)-u(i))^2;
    end
    RMSE(is)=sqrt(msum/(n-1));
end
