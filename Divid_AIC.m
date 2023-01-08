function [AIC,BIC,RMSE]=Divid_AIC(data,distribution,scale,kind)%针对三变量
% Compute the optimal joint distribution of different vegetation types for fitting precipitation-EVI-LST (P-E-L) based on the Akaike Information Criterion (AIC), Bayesian Information Criterion (BIC), and root mean square error (RMSE)
% data is the data of three variables of shrubland, the first is precipitation, the second is enhanced vegetation index (EVI), the third is land surface temperature (LST); distribution is copula function types;scale is time scale;kind is vegetation types.
P=data(:,1);
E=data(:,2);
L=data(:,3);
P1=[]; L1=[]; E1=[];AIC=[];BIC=[];RMSE=[];
for is=1:scale
    P1=[P1,P(is:length(P)-scale+is)];
    L1=[L1,L(is:length(L)-scale+is)];
    E1=[E1,E(is:length(E)-scale+is)];
end
P=sum(P1,2);
L=sum(L1,2);
E=sum(E1,2);
nseas=12;
for i=1:nseas
    tind=i:nseas:length(P);
    P2=P(tind);
    L2=L(tind);
    E2=E(tind);
    n=length(P2);
    if isequal(kind,'gra')
        [a1,a2,a3,b1,b2,c1,c2,c3]=uni_paragra(P2,E2,L2);
        Fp=gevcdf(P2,a1,a2,a3);Fe=logncdf(E2,b1,b2);Fl=gamcdf(L2-c1,c2,c3);
    elseif isequal(kind,'for')
        [a1,a2,a3,b1,b2,b3,c1,c2]=uni_parafor(P2,E2,L2);
        Fp=gamcdf(P2-a1,a2,a3);Fe=gamcdf(E2-b1,b2,b3);Fl=logncdf(L2,c1,c2);
    elseif isequal(kind,'shr')
        [a1,a2,a3,b1,b2,b3,c1,c2,c3]=uni_parashr(P2,E2,L2);
        Fp=gevcdf(P2,a1,a2,a3);Fe=gevcdf(E2,b1,b2,b3);Fl=gevcdf(L2,c1,c2,c3);
    elseif isequal(kind,'cro')
        [a1,a2,a3,b1,b2,b3,c1,c2,c3]=uni_paracro(P2,E2,L2);
        Fp=gevcdf(P2,a1,a2,a3);Fe=gevcdf(E2,b1,b2,b3);Fl=gevcdf(L2,c1,c2,c3);
    end
    if min(Fp)<=0
        f_min=find(Fp==min(Fp));
        Fp(f_min)=min(Fp)+0.00000000000000001;
    elseif    max(Fp)>=1
        f_max=find(Fp==max(Fp));
        Fp(f_max)=max(Fp)-0.0000001;
    end
    if min(Fe)<=0
        f_min=find(Fe==min(Fe));
        Fe(f_min)=min(Fe)+0.00000000000000001;
    elseif    max(Fe)>=1
        f_max=find(Fe==max(Fe));
        Fe(f_max)=max(Fe)-0.0000001;
    end
    if min(Fl)<=0
        f_min=find(Fl==min(Fl));
        Fl(f_min)=min(Fl)+0.00000000000000001;
    elseif    max(Fl)>=1
        f_max=find(Fl==max(Fl));
        Fl(f_max)=max(Fl)-0.0000001;
    end
    EmPro3_123=ecdf3(P2,L2,E2);
    if isequal(distribution,'Gaussian')
        % Gaussian copula
        rho123=Gaussian3_parameter(Fp,Fe,Fl);
        ThPro3_123=copulacdf('Gaussian',[Fp Fe Fl],rho123);
        [AIC3_123, BIC3_123, RMSE3_123]=AIC_BIC_RMSE(EmPro3_123,ThPro3_123,n,3);
    elseif isequal(distribution,'Studentt')
        % Student t copula
        [rho123,NU]=StudentT3_parameter(Fp,Fe,Fl);
        ThPro3_123=copulacdf('t',[Fp Fe Fl],rho123,NU);
        [AIC3_123, BIC3_123, RMSE3_123]=AIC_BIC_RMSE(EmPro3_123,ThPro3_123,n,4);
    elseif isequal(distribution,'Gumbel')
        % Gumbel copula
        paramhat1 = copulafit('Gumbel',[Fp,Fe]);
        c2 = copulacdf('Gumbel',[Fp,Fe],paramhat1);
        paramhat2 = copulafit('Gumbel',[c2,Fl]);
        ThPro3_123 = copulacdf('Gumbel',[c2,Fl],paramhat2);
        [AIC3_123, BIC3_123, RMSE3_123]=AIC_BIC_RMSE(EmPro3_123,ThPro3_123,n,1);
    elseif isequal(distribution,'Frank')
        % Frank copula
        paramhat1 = copulafit('Frank',[Fp,Fe]);
        c2=copulacdf('Frank',[Fp,Fe],paramhat1);
        paramhat2 = copulafit('Frank',[c2,Fl]);
        ThPro3_123=copulacdf('Frank',[c2,Fl],paramhat2);
        [AIC3_123, BIC3_123, RMSE3_123]=AIC_BIC_RMSE(EmPro3_123,ThPro3_123,n,1);
    elseif isequal(distribution,'Clayton')
        % Clayton copula
        paramhat1 = copulafit('Clayton',[Fp,Fe]);
        c2 = copulacdf('Clayton',[Fp,Fe],paramhat1);
        paramhat2 = copulafit('Clayton',[c2,Fl]);
        ThPro3_123= copulacdf('Clayton',[c2,Fl],paramhat2);
        [AIC3_123, BIC3_123, RMSE3_123]=AIC_BIC_RMSE(EmPro3_123,ThPro3_123,n,1);
    end
    AIC(i,1)=AIC3_123;
    BIC(i,1)=BIC3_123;
    RMSE(i,1)=RMSE3_123;
end
end