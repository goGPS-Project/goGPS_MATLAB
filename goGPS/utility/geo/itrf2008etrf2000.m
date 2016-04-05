function [X_out] = itrf2008etrf2000(X_in, date)

% SYNTAX:
%   [X_out] = itrf2008etrf2000(X_in, date);
%
% INPUT:
%   X_in = input coordinates (XYZ)
%   date = reference date
%
% OUTPUT:
%   X_out = output coordinates (XYZ)
%
% DESCRIPTION:
%   ITRF2008 to ETRF2000 coordinate converter.

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

%decimal year (YYYY.DDD)
[~, frac] = date2doy(datenum(date));
t = date(1) + frac;

%ITRF2008 - ETRF2000 parameters
% (http://etrs89.ensg.ign.fr/memo-V8.pdf)

%translation [mm]
T1 =  52.1; T1dot =  0.1;
T2 =  49.3; T2dot =  0.1;
T3 = -58.5; T3dot = -1.8;

%scale factor
D = 1.34e-9; Ddot = 0.08e-9;

%rotation [mas]
R1 =  0.891; R1dot =  0.081;
R2 =  5.390; R2dot =  0.490;
R3 = -8.712; R3dot = -0.792;

%propagate parameters to epoch t
T1 = T1 + T1dot*(t - 2000.0);
T2 = T2 + T2dot*(t - 2000.0);
T3 = T3 + T3dot*(t - 2000.0);
D  = D  + Ddot *(t - 2000.0);
R1 = R1 + R1dot*(t - 2000.0);
R2 = R2 + R2dot*(t - 2000.0);
R3 = R3 + R3dot*(t - 2000.0);

%convert translation parameters to [m]
T1 = T1 * 1e-3;
T2 = T2 * 1e-3;
T3 = T3 * 1e-3;

%convert rotation parameters to [rad]
R1 = R1 * 4.848136e-9;
R2 = R2 * 4.848136e-9;
R3 = R3 * 4.848136e-9;

%translation vector
T = [T1; T2; T3];

%rotation/scale matrix
R = [D -R3 R2; R3 D -R1; -R2 R1 D];

%7-parameters Helmert transformation
X_out = X_in + T + R * X_in;
