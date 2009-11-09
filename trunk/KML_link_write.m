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

%-------------------------------------------------------------------------------
% PARAMETERS FOR THE KML FILE
%-------------------------------------------------------------------------------

%"clampedToGround" plots the points attached to the ground
%"absolute" uses the height defined in the tag <coordinates>;
%N.B. Google Earth uses orthometric heights
z_pos = 'clampedToGround';
%z_pos = 'absolute';

%URL to load the icon for the points
icon = 'http://maps.google.com/mapfiles/kml/pal2/icon10.png';

%point size
scale = 0.5;

%-------------------------------------------------------------------------------
% SCRITTURA DEL FILE
%-------------------------------------------------------------------------------

flink=fopen(filename,'wt');
fprintf(flink, '<?xml version="1.0" standalone="yes"?>\n');
fprintf(flink, '<kml creator="goGPS" xmlns="http://earth.google.com/kml/2.2">\n');
fprintf(flink, '  <Document>\n');
fprintf(flink, '    <open>1</open>\n');
fprintf(flink, '    <ScreenOverlay>\n');
fprintf(flink, '      <name>goGPS logo</name>\n');
fprintf(flink, '      <visibility>1</visibility>\n');
fprintf(flink, '      <Icon>\n');
fprintf(flink, '        <href>ge-goGPS.png</href>\n');
fprintf(flink, '      </Icon>\n');
fprintf(flink, '      <overlayXY x="0" y="1" xunits="fraction" yunits="fraction"/>\n');
fprintf(flink, '      <screenXY x="0" y="1" xunits="fraction" yunits="fraction"/>\n');
fprintf(flink, '      <rotationXY x="0" y="0" xunits="fraction" yunits="fraction"/>\n');
fprintf(flink, '      <size x="0" y="0" xunits="fraction" yunits="fraction"/>\n');
fprintf(flink, '    </ScreenOverlay>\n');
fprintf(flink, '    <NetworkLink>\n');
fprintf(flink, '      <name>goGPS</name>\n');
fprintf(flink, '      <Link id="ID">\n');
fprintf(flink, '        <href>goGPS.kml</href>\n');
fprintf(flink, '        <refreshMode>onInterval</refreshMode>\n');
fprintf(flink, '        <refreshInterval>0.2</refreshInterval>\n');
fprintf(flink, '      </Link>\n');
fprintf(flink, '    </NetworkLink>\n');
fprintf(flink, '    <Placemark>\n');
fprintf(flink, '      <name>Starting location</name>\n');
fprintf(flink, '      <visibility>0</visibility>\n');
fprintf(flink, '      <Point>\n');
fprintf(flink, '        <altitudeMode>%s</altitudeMode>\n',z_pos);
fprintf(flink, '        <coordinates>%.6f,%.6f,%.6f</coordinates>\n',lam,phi,h);
fprintf(flink, '      </Point>\n');
fprintf(flink, '      <Snippet></Snippet>\n');
fprintf(flink, '      <Style>\n');
fprintf(flink, '        <IconStyle>\n');
fprintf(flink, '          <Icon>\n');
fprintf(flink, '            <href>%s</href>\n',icon);
fprintf(flink, '          </Icon>\n');
fprintf(flink, '          <colorMode>normal</colorMode>\n');
fprintf(flink, '          <scale>%.2f</scale>\n',scale);
fprintf(flink, '        </IconStyle>\n');
fprintf(flink, '      </Style>\n');
fprintf(flink, '      <description><![CDATA[ <i>Latitude:</i> %.6f &#176;<br/> <i>Longitude:</i> %.6f &#176;<br/> <i>Elevation:</i> %.1f m]]></description>\n',phi,lam,h);
fprintf(flink, '      <LookAt>\n');
fprintf(flink, '        <longitude>%.6f</longitude>\n',lam);
fprintf(flink, '        <latitude>%.6f</latitude>\n',phi);
fprintf(flink, '        <altitude>0</altitude>\n');
fprintf(flink, '        <range>120</range>\n');
fprintf(flink, '        <tilt>30</tilt>\n');
fprintf(flink, '        <heading>0</heading>\n');
fprintf(flink, '      </LookAt>\n');
fprintf(flink, '    </Placemark>\n');
fprintf(flink, '  </Document>\n</kml>');
fclose(flink);
