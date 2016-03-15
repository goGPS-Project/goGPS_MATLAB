function [corr, distSR_corr] = relativistic_range_error_correction(XR, XS)

% SYNTAX:
%   [corr, distSR_corr] = relativistic_range_error_correction(XR, XS);
%
% INPUT:
%   XR = receiver position  (X,Y,Z)
%   XS = satellite position (X,Y,Z)
%
% OUTPUT:
%   corr = relativistic range error correction term
%   distSR_corr = corrected satellite-receiver distance
%
% DESCRIPTION:
%   Computation of the relativistic range error correction term (Shapiro
%   signal propagation delay).

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) 2009-2014 Mirko Reguzzoni, Eugenio Realini
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

distR = sqrt(sum(XR.^2 ,1));
distS = sqrt(sum(XS.^2 ,2));

XR_mat = XR(:,ones(size(XS,1),1))';
distSR = sqrt(sum((XS-XR_mat).^2 ,2));

% switch char(Eph(31))
%     case 'G'
%         %GM = goGNSS.GM_GPS;
        GM = 3.986005e14;
%     case 'R'
%         %GM = goGNSS.GM_GLO;
%         GM = 3.9860044e14;
%     case 'E'
%         %GM = goGNSS.GM_GAL;
%         GM = 3.986004418e14;
%     case 'C'
%         %GM = goGNSS.GM_BDS;
%         GM = 3.986004418e14;
%     case 'J'
%         %GM = goGNSS.GM_QZS;
%         GM = 3.986005e14;
%     otherwise
%         fprintf('Something went wrong in ecc_anomaly.m\nUnrecongized Satellite system!\n');
%         %GM = goGNSS.GM_GPS;
%         GM = 3.986005e14;
% end

corr = 2*GM/(goGNSS.V_LIGHT^2)*log((distR + distS + distSR)./(distR + distS - distSR));

distSR_corr = distSR + corr;
