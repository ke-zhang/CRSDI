function [a1,a2,a3,b1,b2,c1,c2,c3]=uni_paragra(x,y,z) 
% Calculate the parameters of the Marginal distribution function that best fit the three variables of the grassland.
% x is precipitation; y is enhanced vegetation index (EVI); z is land surface temperature (LST)
gevp=gevfit(x);a1=gevp(1);a2=gevp(2);a3=gevp(3);
lognp=lognfit(y);b1=lognp(1);b2=lognp(2);
[c1,c2,c3]=p3fit(z);
end
