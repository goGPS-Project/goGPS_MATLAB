function rtplot_googleearth_cov (t, pos_R, pos_M, covpos_R, date)

% SYNTAX:
%   rtplot_googleearth_cov (t, pos_R, pos_M, covpos_R, date);
%
% INPUT:
%   t = epoch (t=1,2,...)
%   pos_R = ROVER position (X,Y,Z)
%   pos_M = MASTER position (X,Y,Z)
%   covpos_R = ROVER position covariance matrix
%   date = date expressed as [year,month,day,hour,minutes,seconds)
%
% DESCRIPTION:
%   Real-time visualization on Google Earth.

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

global a_GPS e_GPS
global conf_sat
global link_filename kml_filename
global x_circle

a = a_GPS;
e = e_GPS;

%rover position coordinates X Y Z
XR = pos_R(1);
YR = pos_R(2);
ZR = pos_R(3);

%conversion from cartesian to geodetic coordinates
[phiR, lamR, hR] = cart2geod(XR, YR, ZR);

if (sum(abs(pos_M)) ~= 0)
    %rover position coordinates X Y Z
    XM = pos_M(1);
    YM = pos_M(2);
    ZM = pos_M(3);
    
    %conversion from cartesian to geodetic coordinates
    [phiM, lamM, hM] = cart2geod(XM, YM, ZM);
else
    phiM = 0;
    lamM = 0;
    hM = 0;
end

%computation of the Earth local radius
NM = a / sqrt(1 - e^2 * (sin(phiR))^2);
MM = NM * (1 - e^2) / (1 - e^2 * (sin(phiR))^2);
RM = sqrt(NM*MM);

%conversion into metric coordinates
[EAST, NORTH] = geod2plan(phiR,lamR);

%covariance propagation
covpos_R = global2localCov(covpos_R, pos_R);

T = chol(covpos_R(1:2,1:2));        % Cholesky decomposition
for j = 1 : size(x_circle,1)        % ellipse computation
    x_ellipse(j,:) = x_circle(j,:) * T + [EAST, NORTH];

    %approximate conversion to geodetic coordinates
    delta(1,:) = x_ellipse(j,:) - [EAST, NORTH];
    deltaPhi = delta(1,2) / RM;
    deltaLam = delta(1,1) / (RM * cos(phiR));
    geod_ellipse(j,:) = [deltaLam, deltaPhi] + [lamR, phiR];
end

%conversion from radians to degrees
lamR = lamR*180/pi;
phiR = phiR*180/pi;
lamM = lamM*180/pi;
phiM = phiM*180/pi;
geod_ellipse(:,:) = geod_ellipse(:,:)*180/pi;

%-------------------------------------------------------------------------------
% INITIALIZATION LINK FILE, KML FILE AND GOOGLE EARTH LAUNCH
%-------------------------------------------------------------------------------

if (t == 1)

    %creation of the file to be used as a link to the kml file
    KML_link_write(link_filename,lamR,phiR,hR);

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
    current_path(current_path == '\') = '/';
    [command] = sprintf('"%s/%s"&',current_path,link_filename);
    system(command);
    %pause(5);
end

%-------------------------------------------------------------------------------
% KML FILE WRITING
%-------------------------------------------------------------------------------

%re-writing kml file to show current position in real-time
KML_write_cov(lamR,phiR,hR,lamM,phiM,hM,geod_ellipse,sum(abs(conf_sat)),date);
