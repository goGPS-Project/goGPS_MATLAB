function [pos_S_SP3, dt_S_SP3] = interpolate_SP3(timeb, timee, SP3_time, SP3_coor, SP3_clck, wait_dlg)

% SYNTAX:
%   [pos_S_SP3, dt_S_SP3] = interpolate_SP3(timeb, timee, SP3_time, SP3_coor, SP3_clck, wait_dlg);
%
% INPUT:
%   timeb = beginning of the interpolation timespan (GPS time)
%   timee = end of the interpolation timespan (GPS time)
%   SP3_time = precise ephemeris timestamps (GPS time)
%   SP3_coor = satellite coordinates  [m]
%   SP3_clck = satellite clock errors [s]
%   wait_dlg = optional handler to waitbar figure
%
% OUTPUT:
%   pos_S_SP3 = interpolated coordinates
%   dt_S_SP3  = interpolated clock error
%
% DESCRIPTION:
%   SP3 (precise ephemeris) coordinates 1-second interpolation by Lagrange
%   polynomials. SP3 clock error 1-second interpolation by spline.

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

if (nargin > 5)
    waitbar(0.5,wait_dlg,'Interpolating SP3 (precise ephemeris) coordinates...')
end

%number of seconds in a quarter of an hour
quarter_sec = 900;

%length of the interpolation span
int_len = timee - timeb + 1;

%find the SP3 epoch closest to the beginning of the interpolation timespan
[~, min_p] = min(abs(SP3_time - timeb));

b = SP3_time(min_p) - timeb;

pos_S_SP3 = zeros(3,32, int_len);
dt_S_SP3  = zeros(32, int_len);

for sat = 1 : 32
    
    p = min_p;

    int_X = [];
    int_Y = [];
    int_Z = [];
    int_c = [];
    
    %extract the SP3 coordinates and clocks
    SP3_X = []; SP3_Y = []; SP3_Z = []; SP3_c = [];
    for i = -4 : +4
        SP3_X = [SP3_X SP3_coor(1,sat,p+i)];
        SP3_Y = [SP3_Y SP3_coor(2,sat,p+i)];
        SP3_Z = [SP3_Z SP3_coor(3,sat,p+i)];
        SP3_c = [SP3_c SP3_clck(sat,p+i)];
    end
    
    if (isempty(find(SP3_c >= 999999, 1)))
        
        x  = 1 : 9;
        u  = 4 : 1/quarter_sec : 6-1/quarter_sec;
        
        %spline interpolation (clock)
        y = SP3_c;
        int_c = [int_c interp1(x, y, u, 'spline')];

        u = u - (1/quarter_sec)*int_c;
        
        %Lagrange interpolation (coordinates)
        y1 = SP3_X; y2 = SP3_Y; y3 = SP3_Z;
        int_X = [int_X LagrangeInter(x, y1, u)];
        int_Y = [int_Y LagrangeInter(x, y2, u)];
        int_Z = [int_Z LagrangeInter(x, y3, u)];
        
        while (SP3_time(p+1) < timee)
            
            p = p + 2;
            
            %extract the SP3 coordinates of interest
            SP3_X = []; SP3_Y = []; SP3_Z = []; SP3_c = [];
            for i = -4 : +4
                SP3_X = [SP3_X SP3_coor(1,sat,p+i)];
                SP3_Y = [SP3_Y SP3_coor(2,sat,p+i)];
                SP3_Z = [SP3_Z SP3_coor(3,sat,p+i)];
                SP3_c = [SP3_c SP3_clck(sat,p+i)];
            end
            
            x  = 1 : 9;
            u  = 4 : 1/quarter_sec : 6-1/quarter_sec;
            
            %spline interpolation (clock)
            y = SP3_c;
            int_ctmp = interp1(x, y, u, 'spline');
            int_c = [int_c int_ctmp];
            
            u = u - (1/quarter_sec)*int_ctmp;
            
            %Lagrange interpolation (coordinates)
            y1 = SP3_X; y2 = SP3_Y; y3 = SP3_Z;
            int_X = [int_X LagrangeInter(x, y1, u)];
            int_Y = [int_Y LagrangeInter(x, y2, u)];
            int_Z = [int_Z LagrangeInter(x, y3, u)];
        end
        
        e = SP3_time(p+1) - timee;
        
        int_X(1:(quarter_sec-b)) = [];
        int_Y(1:(quarter_sec-b)) = [];
        int_Z(1:(quarter_sec-b)) = [];
        int_c(1:(quarter_sec-b)) = [];
        
        int_X(end-e+2:end) = [];
        int_Y(end-e+2:end) = [];
        int_Z(end-e+2:end) = [];
        int_c(end-e+2:end) = [];
        
        pos_S_SP3(1,sat,:) = int_X;
        pos_S_SP3(2,sat,:) = int_Y;
        pos_S_SP3(3,sat,:) = int_Z;
        
        dt_S_SP3(sat,:) = int_c;
    end
end

if (nargin > 5)
    waitbar(1,wait_dlg)
end