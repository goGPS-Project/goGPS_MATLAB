function [ztd_corr, pwv_corr] = getZenithDelayCorrection(coo1, p1, t1, coo2, p2,t2)
%
% Compute ztd_corr and iwv_corr at one receiver height from another
%
% Corresponding about the script: alessandra.mascitelli (at) polimi.it
%
% OUTPUT
%   ztd_corr     zenith total delay correction [mm]
%   pwv_corr     integrated water vapor correction [mm]
%
tc2tk = 273.15; %temperature conversion from Celsius to Kelvin
Tisa = 288.15; %temperature International Standard Atmosphere [K]
R = 8.31432;  %universal gas constant for air [J K^-1 mol^-1]
rv = 461.5; %specific gas constant - Water vapor [N m/kg K]
md = 0.0289644; %molar mass of dry air [kg mol^-1]
g = 9.80665; %gravitational acceleration [m/s^2]
k3 = 377600; %constant[K^2/mbar]
k2 = 17; %constant[K/mbar]
%
if nargin<6
%compute p2 from p1
    p2=p1*exp(-(g*md*(coo2(H)-coo1(H)))/(R*Tisa));
%compute t2 from t1
    t2=t1+0.0065(coo1(H)-coo2(H));
end
