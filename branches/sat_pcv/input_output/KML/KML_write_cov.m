function KML_write_cov (lamR,phiR,hR,lamM,phiM,hM,ellipse,nsat,date)

% SYNTAX:
%   KML_write_cov (lamR,phiR,hR,lamM,phiM,hM,ellipse,nsat,date);
%
% INPUT:
%   lamR = rover longitude [degrees]
%   phiR = rover latitude [degrees]
%   hR   = rover orthometric height [m]
%   lamM = master longitude [degrees]
%   phiM = master latitude [degrees]
%   hM   = master orthometric height [m]
%   ellipse = coordinates defining the error ellipse [lam, phi] [degrees]
%   nsat = number of visible satellites
%   date = date expressed as [year,month,day,hour,minutes,seconds)
%
% DESCRIPTION:
%   Write a KML file (Goole Earth).

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

global kml_filename n_sys

%-------------------------------------------------------------------------------
% PARAMETERS FOR THE KML FILE
%-------------------------------------------------------------------------------

%"clampedToGround" plots the points attached to the ground
%"absolute" uses the height defined in the tag <coordinates>;
%N.B. Google Earth uses orthometric heights
z_pos = 'clampToGround';
%z_pos = 'absolute';

%URL to load the icon for the points
iconR = 'http://maps.google.com/mapfiles/kml/pal2/icon26.png';

%URL to load the icon for the master station position
iconM = 'http://maps.google.com/mapfiles/kml/shapes/square.png';

%string representing the ARGB color of the points
min_nsat_LS = 3 + n_sys;
if (nsat >= min_nsat_LS)
    point_colorR = 'fff5005a';
else
    point_colorR = 'ff0000ff';
end

%point size
scaleR = 0.4;

%label color
label_colorR = point_colorR;

%label size
label_scaleR = 0.7;

%string representing the ARGB color of the points
point_colorM = 'ff00ffff';

%point size
scaleM = 0.8;

%label color
label_colorM = point_colorM;

%label size
label_scaleM = 0.7;

%-------------------------------------------------------------------------------
% FORMAT DATE FOR KML
%-------------------------------------------------------------------------------

year = num2str(date(1,1));
month = num2str(date(1,2));
day = num2str(date(1,3));
hour = num2str(date(1,4));
minute = num2str(date(1,5));
second = num2str(floor(date(1,6)));

[null, ncifre] = size(year); %#ok<ASGLU>
if (ncifre == 1)
    [year] = sprintf('200%s',year);
elseif (ncifre == 2)
    [year] = sprintf('20%s',year);
end

[null, ncifre] = size(month); %#ok<ASGLU>
if (ncifre == 1)
    [month] = sprintf('0%s',month);
end

[null, ncifre] = size(day); %#ok<ASGLU>
if (ncifre == 1)
    [day] = sprintf('0%s',day);
end

[null, ncifre] = size(hour); %#ok<ASGLU>
if (ncifre == 1)
    [hour] = sprintf('0%s',hour);
end

[null, ncifre] = size(minute); %#ok<ASGLU>
if (ncifre == 1)
    [minute] = sprintf('0%s',minute);
end

[null, ncifre] = size(second); %#ok<ASGLU>
if (ncifre == 1)
    [second] = sprintf('0%s',second);
end

%-------------------------------------------------------------------------------
% FILE WRITING
%-------------------------------------------------------------------------------

fkml=fopen(kml_filename,'wt');

while (fkml == -1)
    fkml=fopen(kml_filename,'wt');
end

pos = find(kml_filename == '/');
kml_name = kml_filename(pos(end)+1:end);

fprintf(fkml, '<?xml version="1.0" encoding="UTF-8"?>\n');
fprintf(fkml, '<kml xmlns="http://www.opengis.net/kml/2.2" xmlns:gx="http://www.google.com/kml/ext/2.2" xmlns:kml="http://www.opengis.net/kml/2.2" xmlns:atom="http://www.w3.org/2005/Atom">\n');
fprintf(fkml, '<Document>\n');
fprintf(fkml, '\t<name>%s</name>\n', kml_name);
fprintf(fkml, '\t<snippet>created by goGPS</snippet>\n');
fprintf(fkml, '\t\t<Style id="go1">\n');
fprintf(fkml, '\t\t\t<IconStyle>\n');
fprintf(fkml, '\t\t\t\t<color>%s</color>\n',point_colorR);
fprintf(fkml, '\t\t\t\t<colorMode>normal</colorMode>\n');
fprintf(fkml, '\t\t\t\t<scale>%.2f</scale>\n',scaleR);
fprintf(fkml, '\t\t\t\t<Icon>\n');
fprintf(fkml, '\t\t\t\t\t<href>%s</href>\n',iconR);
fprintf(fkml, '\t\t\t\t</Icon>\n');
fprintf(fkml, '\t\t\t</IconStyle>\n');
fprintf(fkml, '\t\t\t<LabelStyle>\n');
fprintf(fkml, '\t\t\t\t<color>%s</color>\n',label_colorR);
fprintf(fkml, '\t\t\t\t<scale>%s</scale>\n',label_scaleR);
fprintf(fkml, '\t\t\t</LabelStyle>\n');
fprintf(fkml, '\t\t</Style>\n');
fprintf(fkml, '\t\t<Style id="master">\n');
fprintf(fkml, '\t\t\t<IconStyle>\n');
fprintf(fkml, '\t\t\t\t<color>%s</color>\n',point_colorM);
fprintf(fkml, '\t\t\t\t<colorMode>normal</colorMode>\n');
fprintf(fkml, '\t\t\t\t<scale>%.2f</scale>\n',scaleM);
fprintf(fkml, '\t\t\t\t<Icon>\n');
fprintf(fkml, '\t\t\t\t\t<href>%s</href>\n',iconM);
fprintf(fkml, '\t\t\t\t</Icon>\n');
fprintf(fkml, '\t\t\t</IconStyle>\n');
fprintf(fkml, '\t\t\t<LabelStyle>\n');
fprintf(fkml, '\t\t\t\t<color>%s</color>\n',label_colorM);
fprintf(fkml, '\t\t\t\t<scale>%s</scale>\n',label_scaleM);
fprintf(fkml, '\t\t\t</LabelStyle>\n');
fprintf(fkml, '\t\t</Style>\n');
fprintf(fkml, '\t\t<Style id="ellipse">\n');
fprintf(fkml, '\t\t\t<LineStyle>\n');
fprintf(fkml, '\t\t\t<color>ff0000ff</color>\n');
fprintf(fkml, '\t\t\t<width>2</width>\n');
fprintf(fkml, '\t\t\t</LineStyle>\n');
fprintf(fkml, '\t\t\t<PolyStyle>\n');
fprintf(fkml, '\t\t\t<color>330000ff</color>\n');
fprintf(fkml, '\t\t\t</PolyStyle>\n');
fprintf(fkml, '\t\t</Style>\n');
if (lamM ~= 0 | phiM ~= 0 | hM ~= 0)
    fprintf(fkml, '\t\t<Placemark>\n');
    fprintf(fkml, '\t\t\t<name>Master station</name>\n');
    fprintf(fkml, '\t\t\t<description><![CDATA[ <i>Latitude:</i> %.8f &#176;<br/> <i>Longitude:</i> %.8f &#176;<br/> <i>Elevation:</i> %.1f m<br/> <i>Time:</i> %s-%s-%s %s:%s:%s]]></description>\n',phiM,lamM,hM,year,month,day,hour,minute,second);
    fprintf(fkml, '\t\t\t<styleUrl>#master</styleUrl>\n');
    fprintf(fkml, '\t\t\t<Point>\n');
    fprintf(fkml, '\t\t\t\t<altitudeMode>%s</altitudeMode>\n',z_pos);
    fprintf(fkml, '\t\t\t\t<coordinates>%.8f,%.8f,%.3f</coordinates>\n',lamM,phiM,hM);
    fprintf(fkml, '\t\t\t</Point>\n');
    fprintf(fkml, '\t\t</Placemark>\n');
end
fprintf(fkml, '\t\t<Placemark>\n');
fprintf(fkml, '\t\t\t<name>%d</name>\n', nsat);
fprintf(fkml, '\t\t\t<description><![CDATA[ <i>Latitude:</i> %.8f &#176;<br/> <i>Longitude:</i> %.8f &#176;<br/> <i>Elevation:</i> %.1f m<br/> <i>Satellites:</i> %d <br/> <i>Time:</i> %s-%s-%s %s:%s:%s]]></description>\n',phiR,lamR,hR,nsat,year,month,day,hour,minute,second);
fprintf(fkml, '\t\t\t<styleUrl>#go1</styleUrl>\n');
fprintf(fkml, '\t\t\t<Point>\n');
fprintf(fkml, '\t\t\t\t<altitudeMode>%s</altitudeMode>\n',z_pos);
fprintf(fkml, '\t\t\t\t<coordinates>%.8f,%.8f,%.3f</coordinates>\n',lamR,phiR,hR);
fprintf(fkml, '\t\t\t</Point>\n');
fprintf(fkml, '\t\t</Placemark>\n');
fprintf(fkml, '\t\t<Placemark>\n');
fprintf(fkml, '\t\t\t<name>Error ellipse</name>\n');
fprintf(fkml, '\t\t\t<styleUrl>#ellipse</styleUrl>\n');
fprintf(fkml, '\t\t\t<Polygon>\n');
fprintf(fkml, '\t\t\t\t<extrude>0</extrude>\n');
fprintf(fkml, '\t\t\t\t<tessellate>0</tessellate>\n');
fprintf(fkml, '\t\t\t\t<altitudeMode>%s</altitudeMode>\n',z_pos);
fprintf(fkml, '\t\t\t\t<outerBoundaryIs>\n');
fprintf(fkml, '\t\t\t\t\t<LinearRing>\n');
fprintf(fkml, '\t\t\t\t\t\t<coordinates>\n');
for i = 1 : size(ellipse,1)
    fprintf(fkml, '\t\t\t\t\t\t%.15f,%.15f,0\n',ellipse(i,1),ellipse(i,2));
end
fprintf(fkml, '\t\t\t\t\t\t</coordinates>\n');
fprintf(fkml, '\t\t\t\t\t</LinearRing>\n');
fprintf(fkml, '\t\t\t\t</outerBoundaryIs>\n');
fprintf(fkml, '\t\t\t</Polygon>\n');
fprintf(fkml, '\t\t</Placemark>\n');
fprintf(fkml, '</Document>\n</kml>');
fclose(fkml);
