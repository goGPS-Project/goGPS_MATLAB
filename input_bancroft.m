function [pos, B] = input_bancroft(pr, sv, time, Eph)

% SYNTAX:
%   [pos, B] = input_bancroft(pr, sv, time, Eph);
%
% INPUT:
%   pr = observed code pseudoranges
%   sv = visible satellites configuration
%   time = GPS time
%   Eph = ephemerides matrix
%
% OUTPUT:
%   pos = ROVER ground position (X,Y,Z)
%   B = Bancroft matrix
%
% DESCRIPTION:
%   Prepare input for Bancroft algorithm and run the algorithm.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.1 alpha
%
% Copyright (C) Kai Borre 
% Kai Borre and C.C. Goad 11-24-96
%
% Adapted by Mirko Reguzzoni, Eugenio Realini, 2009
%
%----------------------------------------------------------------------------------------------

global v_light

%number of pseudorange observations
m = size(pr,1);

%first guess Earth center
pos = zeros(4,1);

%Bancroft matrix initialization
B = [];

%loop on satellites
for jsat = 1:m

    %satellite position correction (clock and Earth rotation)
    [X, tcorr] = sat_corr(Eph, sv(jsat), time, pr(jsat));

    corrected_pseudorange = pr(jsat) + v_light * tcorr;
    B(jsat,1) = X(1);
    B(jsat,2) = X(2);
    B(jsat,3) = X(3);
    B(jsat,4) = corrected_pseudorange;
end

%approximate position by Bancroft algorithm
pos = bancroft(B);  

%clock offset computation
clock_offset = pos(4) / v_light;
pos(4) = clock_offset * 1.e+9;
