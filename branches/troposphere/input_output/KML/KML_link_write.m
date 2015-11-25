function KML_link_write (filename, lam, phi, h)

% SYNTAX:
%   KML_write (filename,lam,phi,h,nsat,date);
%
% INPUT:
%   filename = name of the file with extension
%   lam = longitude
%   phi = latitude
%   h = orthometric height
%
% DESCRIPTION:
%   Write a link to a KML file (Goole Earth).

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

%-------------------------------------------------------------------------------
% PARAMETERS FOR THE KML FILE
%-------------------------------------------------------------------------------

%"clampedToGround" plots the points attached to the ground
%"absolute" uses the height defined in the tag <coordinates>;
%N.B. Google Earth uses orthometric heights
z_pos = 'clampToGround';
%z_pos = 'absolute';

%URL to load the icon for the points
icon = 'http://maps.google.com/mapfiles/kml/pal2/icon10.png';

%point size
scale = 0.5;

%-------------------------------------------------------------------------------
% SCRITTURA DEL FILE
%-------------------------------------------------------------------------------

flink=fopen(filename,'wt');
fprintf(flink, '<?xml version="1.0" encoding="UTF-8"?>\n');
fprintf(flink, '<kml xmlns="http://www.opengis.net/kml/2.2" xmlns:gx="http://www.google.com/kml/ext/2.2" xmlns:kml="http://www.opengis.net/kml/2.2" xmlns:atom="http://www.w3.org/2005/Atom">\n');
fprintf(flink, '<Document>\n');
fprintf(flink, '\t<open>1</open>\n');
fprintf(flink, '\t<ScreenOverlay>\n');
fprintf(flink, '\t\t<name>goGPS logo</name>\n');
fprintf(flink, '\t\t<visibility>1</visibility>\n');
fprintf(flink, '\t\t<Icon>\n');
fprintf(flink, '\t\t\t<href>./ge-goGPS.png</href>\n');
fprintf(flink, '\t\t</Icon>\n');
fprintf(flink, '\t\t<overlayXY x="0" y="1" xunits="fraction" yunits="fraction"/>\n');
fprintf(flink, '\t\t<screenXY x="0" y="1" xunits="fraction" yunits="fraction"/>\n');
fprintf(flink, '\t\t<rotationXY x="0" y="0" xunits="fraction" yunits="fraction"/>\n');
fprintf(flink, '\t\t<size x="0" y="0" xunits="fraction" yunits="fraction"/>\n');
fprintf(flink, '\t</ScreenOverlay>\n');
fprintf(flink, '\t<NetworkLink>\n');
fprintf(flink, '\t\t<name>goGPS</name>\n');
fprintf(flink, '\t\t<Link id="ID">\n');
fprintf(flink, '\t\t\t<href>./goGPS.kml</href>\n');
fprintf(flink, '\t\t\t<refreshMode>onInterval</refreshMode>\n');
fprintf(flink, '\t\t\t<refreshInterval>0.5</refreshInterval>\n');
fprintf(flink, '\t\t</Link>\n');
fprintf(flink, '\t</NetworkLink>\n');
fprintf(flink, '\t<Placemark>\n');
fprintf(flink, '\t\t<name>Starting location</name>\n');
fprintf(flink, '\t\t<visibility>0</visibility>\n');
fprintf(flink, '\t\t<description><![CDATA[ <i>Latitude:</i> %.8f &#176;<br/> <i>Longitude:</i> %.8f &#176;<br/> <i>Elevation:</i> %.3f m]]></description>\n',phi,lam,h);
fprintf(flink, '\t\t<LookAt>\n');
fprintf(flink, '\t\t\t<longitude>%.8f</longitude>\n',lam);
fprintf(flink, '\t\t\t<latitude>%.8f</latitude>\n',phi);
fprintf(flink, '\t\t\t<altitude>0</altitude>\n');
fprintf(flink, '\t\t\t<range>120</range>\n');
fprintf(flink, '\t\t\t<tilt>30</tilt>\n');
fprintf(flink, '\t\t\t<heading>0</heading>\n');
fprintf(flink, '\t\t</LookAt>\n');
fprintf(flink, '\t\t<Style>\n');
fprintf(flink, '\t\t\t<IconStyle>\n');
fprintf(flink, '\t\t\t\t<colorMode>normal</colorMode>\n');
fprintf(flink, '\t\t\t\t<scale>%.2f</scale>\n',scale);
fprintf(flink, '\t\t\t\t<Icon>\n');
fprintf(flink, '\t\t\t\t\t<href>%s</href>\n',icon);
fprintf(flink, '\t\t\t\t</Icon>\n');
fprintf(flink, '\t\t\t</IconStyle>\n');
fprintf(flink, '\t\t</Style>\n');
fprintf(flink, '\t\t<Point>\n');
fprintf(flink, '\t\t\t<altitudeMode>%s</altitudeMode>\n',z_pos);
fprintf(flink, '\t\t\t<coordinates>%.8f,%.8f,%.3f</coordinates>\n',lam,phi,h);
fprintf(flink, '\t\t</Point>\n');
fprintf(flink, '\t</Placemark>\n');
fprintf(flink, '</Document>\n</kml>');
fclose(flink);
