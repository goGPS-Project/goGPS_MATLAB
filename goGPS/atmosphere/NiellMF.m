function hmf = NiellMF(DOY,lat,h,elev)

% Niell Mapping Function
%
% %%%%%%%%%% HELP %%%%%%%%%%
%
%hmf = NiellMF(DOY,lat,h,elev)
%
% Input Data
% DOY is day of year.
% lat is latitude, degrees.
% long is longitude, degrees.
% h is height, meters.
% elev is elevation angle, degrees.
%
% Transwritten by Phakphong Homniam
% September 20, 2002
% Original Mathcad source code by Boonsap Witchayangkoon, 2000
% http://gps-ppp.blogspot.com/2009/05/neill-mapping-function-source-code-for.html

% Coefficient for changing DOY to radian
day2rad = 2*pi/365.25;

% Coefficient for changing degrees to radians
deg2rad = pi/180;

% Latitude array for the hydrostatic mapping function coefficients
lat_hmf = [15,30,45,60,75]';

% Average of coefficients a, b and c corresponding to the given latitude.
% See Table 5.4
abc_avg = [1.2769934 2.9153695 62.610505
    1.2683230 2.9152299 62.837393
    1.2465397 2.9288445 63.721774
    1.2196049 2.9022565 63.824265
    1.2045996 2.9024912 64.258455]*10^-3;

% Amplitude of coefficients a, b and c corresponding to the given latitude.
abc_amp = [0 0 0
    1.2709626 2.1414979 9.0128400
    2.6523662 3.0160779 4.3497037
    3.4000452 7.2562722 84.795348
    4.1202191 11.723375 170.37206]*10^-5;

% Height correction
a_ht = 2.53*10^-5;
b_ht = 5.49*10^-3;
c_ht = 1.14*10^-3;

% Wet Mapping Function
% See Table5.5
lat_wmf = lat_hmf;
abc_w2po = [5.8021897*10^-4 1.4275268*10^-3 4.3472961*10^-2
    5.6794847*10^-4 1.5138625*10^-3 4.6729510*10^-2
    5.8118019*10^-4 1.4572752*10^-3 4.3908931*10^-2
    5.9727542*10^-4 1.5007428*10^-3 4.4626982*10^-2
    6.1641693*10^-4 1.7599082*10^-3 5.4736038*10^-2];

hs_km = h/1000;
DOY_atm = DOY-28;
if lat < 0
    DOY_atm = DOY_atm+365.25/2;
    lat = abs(lat);
end
DOYr_atm = DOY_atm*day2rad;
cost = cos(DOYr_atm);
if lat <= lat_hmf(1,1)
    a = abc_avg(1,1);
    b = abc_avg(1,2);
    c = abc_avg(1,3);
end
for i = 1:4
    if (lat >= lat_hmf(i,1)) & (lat <= lat_hmf(i+1,1))
        dlat = (lat-lat_hmf(i))/(lat_hmf(i+1)-lat_hmf(i));
        daavg = abc_avg(i+1,1)-abc_avg(i,1);
        daamp = abc_amp(i+1,1)-abc_amp(i,1);
        aavg = abc_avg(i,1)+dlat*daavg;
        aamp = abc_amp(i,1)+dlat*daamp;
        a = aavg+aamp*cost; % + or - i don't sure
        dbavg = abc_avg(i+1,2)-abc_avg(i,2);
        dbamp = abc_amp(i+1,2)-abc_amp(i,2);
        bavg = abc_avg(i,2)+dlat*dbavg;
        bamp = abc_amp(i,2)+dlat*dbamp;
        b = bavg+bamp*cost; % + or - i don't sure
        dcavg = abc_avg(i+1,3)-abc_avg(i,3);
        dcamp = abc_amp(i+1,3)-abc_amp(i,3);
        cavg = abc_avg(i,3)+dlat*dcavg;
        camp = abc_amp(i,3)+dlat*dcamp;
        c = aavg+aamp*cost; % + or - i don't sure
    end
end
if lat >= lat_hmf(5,1)
    a = abc_avg(5,1)+abc_amp(5,1)*cost; % + or - i don't sure
    b = abc_avg(5,2)+abc_amp(5,2)*cost; % + or - i don't sure
    c = abc_avg(5,3)+abc_amp(5,3)*cost; % + or - i don't sure
end
sine = sin(elev*deg2rad);
cose = cos(elev*deg2rad);

beta = b/(sine+c);
gamma = a/(sine+beta);
topcon = 1+(a/(1+b/(1+c)));
hmf(1,1) = topcon/(sine+gamma);
hmf(2,1) = -topcon*cose/((sine+gamma)^2*...
    (1-a/((sine+beta)^2*(1-b/(sine*c)^2))));

beta = b_ht/(sine+c_ht);
gamma = a_ht/(sine+beta);
topcon = 1+a_ht/(1+b_ht/(1+c_ht));
ht_corr_coef = 1/sine-topcon/(sine+gamma);
ht_corr = ht_corr_coef*hs_km;
hmf(1,1) = hmf(1,1)+ht_corr;
dhcc_del = -cose/sine^2+topcon*cose/((sine+gamma)^2*...
    (1-a_ht/((sine+beta)^2*(1-b_ht/(sine*c_ht)^2))));
dht_corr_del = dhcc_del*hs_km;
hmf(2,1) = hmf(2,1)+dht_corr_del;
dht_corr_del = dhcc_del*hs_km;
hmf(2,1) = hmf(2,1)+dht_corr_del;
if lat <= lat_wmf(1,1)
    alat = abc_w2po(1,1);
    blat = abc_w2po(2,1);
    clat = abc_w2po(3,1);
end
for i = 1:4
    if (lat >= lat_wmf(i,1)) & (lat <= lat_wmf(i+1,1))
        dll = (lat-lat_wmf(i,1))/(lat_wmf(i+1,1)-lat_wmf(i,1));
        da = abc_w2po(i+1,1)-abc_w2po(i,1);
        alat = abc_w2po(i,1)+dll*da;
        db = abc_w2po(i+1,2)-abc_w2po(i,2);
        blat = abc_w2po(i,2)+dll*db;
        dc = abc_w2po(i+1,3)-abc_w2po(i,3);
        clat = abc_w2po(i,3)+dll*dc;
    end
end
if lat >= lat_wmf(5,1)
    alat = abc_w2po(5,1);
    blat = abc_w2po(5,2);
    clat = abc_w2po(5,3);
end
sinelat = sin(elev*deg2rad);
coselat = cos(elev*deg2rad);
betalat = blat/(sinelat+clat);
gammalat = alat/(sinelat+betalat);
topconlat = 1+alat/(1+blat/(1+clat));
wmf(1,1) = topconlat/(sinelat+gammalat);
wmf(2,1) = -topconlat/((sinelat+gammalat)^2*(coselat-alat/...
    ((sinelat+betalat)^2*coselat*(1-blat/(sinelat*clat)^2))));
hmf(3,1) = wmf(1,1);
hmf(4,1) = wmf(2,1);
