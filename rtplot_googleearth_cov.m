function rtplot_googleearth_cov (t, pos_R, covpos_R, link_filename, kml_filename, data)

% SYNTAX:
%   rtplot_googleearth_cov (t, pos_R, covpos_R, link_filename, kml_filename, data);
%
% INPUT:
%   t = epoch (t=1,2,...)
%   pos_R = ROVER position (X,Y,Z)
%   covpos_R = ROVER position covariance matrix
%   link_filename = name of the file used for the Google Earth link
%   kml_filename = name of the .KML file
%   data = date expressed as [year,month,day,hour,minutes,seconds)
%
% DESCRIPTION:
%   Real-time representation of the ROVER position using Google Earth.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.1 pre-alpha
%
% Copyright (C) 2009 Mirko Reguzzoni*, Eugenio Realini**
%
% * Laboratorio di Geomatica, Polo Regionale di Como, Politecnico di Milano, Italy
% ** Media Center, Osaka City University, Japan
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

global a e
global conf_sat
global GE_path GE_append
global x_circle

%rover position coordinates X Y Z
X = pos_R(1);
Y = pos_R(2);
Z = pos_R(3);

%conversion from cartesian to geodetic coordinates
[phi, lam, h] = cart2geod(X, Y, Z);

%computation of the Earth local radius
NM = a / sqrt(1 - e^2 * (sin(phi))^2);
MM = NM * (1 - e^2) / (1 - e^2 * (sin(phi))^2);
RM = sqrt(NM*MM);

%conversion into metric coordinates
[EST, NORD] = geod2plan(phi,lam);

%covariance propagation
covpos_R = global2localCov(covpos_R, pos_R);

T = chol(covpos_R(1:2,1:2));        % Cholesky decomposition
for j = 1 : size(x_circle,1)        % ellipse computation
    x_ellipse(j,:) = x_circle(j,:) * T + [EST, NORD];

    %approximate conversion to geodetic coordinates
    delta(1,:) = x_ellipse(j,:) - [EST, NORD];
    deltaPhi = delta(1,2) / RM;
    deltaLam = delta(1,1) / (RM * cos(phi));
    geod_ellipse(j,:) = [deltaLam, deltaPhi] + [lam, phi];
end

%conversion from radians to degrees
lam = lam*180/pi;
phi = phi*180/pi;
geod_ellipse(:,:) = geod_ellipse(:,:)*180/pi;

%-------------------------------------------------------------------------------
% INITIALIZATION LINK FILE, KML FILE AND GOOGLE EARTH LAUNCH
%-------------------------------------------------------------------------------

if (t == 1)

    %creation of the file to be used as a link to the kml file
    KML_link_write(link_filename,lam,phi,h);

    %initialization of the kml file
    fkml=fopen(kml_filename, 'wt');
    fprintf(fkml, '<?xml version="1.0" standalone="yes"?>\n');
    fprintf(fkml, '<kml creator="goGPS" xmlns="http://earth.google.com/kml/2.2">\n');
    fprintf(fkml, '  <Document>\n');
    fprintf(fkml, '    <name><![CDATA[%s]]></name>\n', kml_filename);
    fprintf(fkml, '    <Snippet><![CDATA[created by goGPS]]></Snippet>\n');
    fprintf(fkml, '  </Document>\n</kml>');
    fclose(fkml);

    %run google earth (in background)
    current_path = pwd;
    current_path(find(current_path == '\')) = '/';
    [command] = sprintf([GE_path ' %s/%s &'],current_path,link_filename);
    system(command);
    %pause(5);
end

%-------------------------------------------------------------------------------
% KML FILE UPDATE / RE-WRITE
%-------------------------------------------------------------------------------

if (GE_append == 0)
    %re-writing kml file shows just the current position
    KML_write_cov(kml_filename,lam,phi,h,geod_ellipse,sum(abs(conf_sat)),data);
else
   %updating keeps memory of all the plotted positions (kml file can become too big)
   KML_update(kml_filename,lam,phi,h,sum(abs(conf_sat)),data);
end
