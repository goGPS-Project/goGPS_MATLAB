function [dt_S_SP3] = interpolate_SP3_clck(time, SP3_time, SP3_clck)

% SYNTAX:
%   [dt_S_SP3] = interpolate_SP3_clck(time, SP3_time, SP3_clck);
%
% INPUT:
%   timeb = beginning of the interpolation timespan (GPS time)
%   timee = end of the interpolation timespan (GPS time)
%   SP3_time = precise ephemeris epochs (GPS time)
%   SP3_clck = satellite clock errors [s]
%   sat = satellite PRN
%
% OUTPUT:
%   dt_S_SP3  = interpolated clock correction
%
% DESCRIPTION:
%   SP3 (precise ephemeris) clock correction 1-second interpolation by spline.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.2.0 beta
%
% Copyright (C) 2009-2011 Mirko Reguzzoni, Eugenio Realini
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

%number of seconds in a quarter of an hour
quarter_sec = 900;

%find the SP3 epoch closest to the interpolation time
[~, p] = min(abs(SP3_time - time));

b = SP3_time(p) - time;

%extract the SP3 clocks
if (b>0)
    SP3_c = [SP3_clck(p-1) SP3_clck(p)];
    u = 1 - b/quarter_sec;
else
    SP3_c = [SP3_clck(p) SP3_clck(p+1)];
    u = -b/quarter_sec;
end

dt_S_SP3  = 0;

if (isempty(find(SP3_c >= 999999, 1)))

    %linear interpolation (clock)
    dt_S_SP3 = (1-u)*SP3_c(1) + u*SP3_c(2);

%     plot([0 1],SP3_c,'o',u,dt_S_SP3,'.')
%     pause
end
