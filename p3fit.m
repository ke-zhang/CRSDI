function [a0,A,B]=p3fit(data)
% Parameter estimation of Pearson III (P-III) distribution
xa=mean(data);xst=std(data);
       cs=skewness(data,0);
        A=(2/cs)^2;
        B=sqrt(xst^2/A);% p3分布矩法估计参数
        a0=xa-A*B;
end
