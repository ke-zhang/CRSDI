function [a1,a2,a3,b1,b2,b3,c1,c2,c3]=uni_paracro(x,y,z) 
% Calculate the parameters of the Marginal distribution function that best fit the three variables of the cropland.
% x is precipitation; y is enhanced vegetation index (EVI); z is land surface temperature (LST)
gevp=gevfit(x);a1=gevp(1);a2=gevp(2);a3=gevp(3);
geve=gevfit(y);b1=geve(1);b2=geve(2);b3=geve(3);
gevl=gevfit(z);c1=gevl(1);c2=gevl(2);c3=gevl(3);
end
