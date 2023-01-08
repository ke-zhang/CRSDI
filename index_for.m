function [DI]=index_for(data,scale)
% Calculate the comprehensive remote sensing drought index (CRSDI) of forest land
% data is the data of three variables of forest land, the first is precipitation, the second is enhanced vegetation index (EVI), the third is land surface temperature (LST); scale is time scale.
P=data(:,1);
E=data(:,2);
L=data(:,3);
P1=[]; L1=[]; E1=[];
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
    [a1,a2,a3,b1,b2,b3,c1,c2]=uni_parafor(P2,E2,L2);
        Fp=gamcdf(P2-a1,a2,a3);Fe=gamcdf(E2-b1,b2,b3);Fl=logncdf(L2,c1,c2);
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
    paramhat = copulafit('Gumbel',[Fp,Fe]);
    ThPro2_13 = copulacdf('Gumbel',[Fp,Fe],paramhat);
    paramhat2 = copulafit('Gumbel',[ThPro2_13,Fl]);
    ThPro3_123 = copulacdf('Gumbel',[ThPro2_13,Fl],paramhat2);
    ThPro=ThPro2_13-ThPro3_123;
    T(tind,1)=ThPro;
end
n=length(T);
R=tiedrank(T);
PG=(R-0.44)./(n+0.12);% Gringorten plotting position formulae
DI=norminv(PG);
end
