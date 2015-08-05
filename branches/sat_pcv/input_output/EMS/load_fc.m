function [prc_E, GPS_time02] = load_fc(iodp_mask, prn_mask, MT, msg, GPS_time)

% SYNTAX:
%   [prc_E, GPS_time02] = load_fc(iodp_mask, prn_mask, MT, msg, GPS_time);
%
% INPUT:
%   iodp_mask = IODP masks [vector]
%   prn_mask  = PRN masks [matrix]
%   MT  = message types [vector]
%   msg = EGNOS message strings [matrix]
%   GPS_time = GPS time of the messages
%
% OUTPUT:
%   prc_E = pseudorange corrections (PRC) [matrix]
%   GPS_time02 = GPS time of the pseudorange corrections
%
% DESCRIPTION:
%   Load the Fast Corrections (FC).

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) 2009-2014 Mirko Reguzzoni, Eugenio Realini
%
% Code contributed by Giuliano Sironi, 2011
%----------------------------------------------------------------------------------------------
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
%----------------------------------------------------------------------------------------------
         
%keep only the MTs that contain the PRC
%WARNING: the messages are not provided at regular intervals
r_MT = find(MT == 2 | MT == 3 | MT == 4 | MT == 24 | MT == 0);

%number of message types
n02 = length(r_MT);

%initialization
IODF  = NaN(n02,1);
PRC   = NaN(n02,13);
UDREI = NaN(n02,13);
PRN   = NaN(n02,13);
n_SV  = NaN(n02,1);

%GPS time of the messages
GPS_time02 = GPS_time(r_MT,:);

for i = 1 : n02
    
    [mt, iodf, prc, udrei, sv, n_sv] = ems2prc(msg(r_MT(i),:), iodp_mask, prn_mask); %#ok<ASGLU>
    
    %MT(i)     = mt;
    IODF(i)    = iodf;
    PRC(i,:)   = prc;   %pseudorange correction vector
    UDREI(i,:) = udrei;
    PRN(i,:)   = sv;
    n_SV(i)    = n_sv;
end

nGPSsat = 32;
prc_E = zeros(n02, nGPSsat);

%initialize the PRC matrix
for k = 1 : n_SV(1)
    %i_prc = find(SV == PRN(1, k));
    if PRN(1, k) > 0
        if UDREI(1,k) < 14
            prc_E(1, PRN(1, k)) = PRC(1,k);
        end
    end
end

%store the corrections for the next epochs
for i = 2 : n02
    prc_E(i,:) = prc_E(i-1,:);
    for j = 1 : n_SV(i)
        if PRN(i, j) > 0
            if UDREI(i,j) < 14
                prc_E(i, PRN(i, j)) = PRC(i,j);
            else 
                prc_E(i, PRN(i, j)) = 0;
            end
        end
    end
end
