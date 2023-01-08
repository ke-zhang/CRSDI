function [a1,a2,a3,b1,b2,b3,c1,c2]=uni_parafor(x,y,z) 
% Calculate the parameters of the Marginal distribution function that best fit the three variables of the forest land.
% x is precipitation; y is enhanced vegetation index (EVI); z is land surface temperature (LST)
[a1,a2,a3]=p3fit(x);
[b1,b2,b3]=p3fit(y);
lognp=lognfit(z);c1=lognp(1);c2=lognp(2);
end
