function [posR, dtR] = input_bancroft(pr, posS, dtS)

% SYNTAX:
%   [posR, dtR] = input_bancroft(pr, posS, dtS);
%
% INPUT:
%   pr = observed code pseudoranges
%   posS = satellite position (X, Y, Z)
%   dtS = satellite clock correction
%
% OUTPUT:
%   posR = ROVER ground position (X,Y,Z)
%   dtR = receiver clock error
%
% DESCRIPTION:
%   Prepare input for Bancroft algorithm and run the algorithm.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.1.3 alpha
%
% Copyright (C) Kai Borre
% Kai Borre and C.C. Goad 11-24-96
%
% Adapted by Mirko Reguzzoni, Eugenio Realini, 2009
%----------------------------------------------------------------------------------------------

global v_light

%number of pseudorange observations
m = size(pr,1);

%Bancroft matrix initialization
B = [];

%loop on satellites
for jsat = 1:m

    if (~isempty(posS))
        corrected_pseudorange = pr(jsat) + v_light * dtS(jsat);
        B(jsat,1) = posS(jsat,1);
        B(jsat,2) = posS(jsat,2);
        B(jsat,3) = posS(jsat,3);
        B(jsat,4) = corrected_pseudorange;
    end
end

%approximate position by Bancroft algorithm
pos = bancroft(B);
posR = pos(1:3);

%clock offset computation
dtR = pos(4) / v_light;
