%display information
fprintf('Writing KML file...\n');
%"clampToGround" plots the points attached to the ground
%"absolute" uses the height defined in the tag <coordinates>;
%N.B. Google Earth uses orthometric heights
z_pos = 'clampToGround';
%z_pos = 'absolute';
%URL to load the icon for the points
iconR = 'http://maps.google.com/mapfiles/kml/pal2/icon26.png';
iconM = 'http://maps.google.com/mapfiles/kml/shapes/square.png';
iconP = 'http://maps.google.com/mapfiles/kml/shapes/square.png';
good_point_colorR = 'fff5005a';
bad_point_colorR = 'ff0000ff';
dyn_point_colorR = 'ff00ffff';
point_colorM = 'ff00ffff';
point_colorP = 'ff32cd32';
%point size
scaleR = 0.2;
scaleM = 0.8;
scaleP = 0.8;
%line color and thickness (rover track)
line_colorR = 'fff5005a';
line_widthR = 1;
%line color and thickness (stop-go-stop direction)
line_colorG = 'ff0000ff';
line_widthG = 4;
%label color
label_colorM = point_colorM;
label_colorP = point_colorP;
%label size
label_scaleM = 0.7;
label_scaleP = 0.7;
%initialization

phiM = zeros(1, nSol / (1 + state.isForwardBackwardKF()));
lamM = zeros(1, nSol / (1 + state.isForwardBackwardKF()));
hM = zeros(1, nSol / (1 + state.isForwardBackwardKF()));

%threshold on KHDOP
if (o1 == 1)
    KHDOP_thres = median(KHDOP);
else
    KHDOP_thres = 2;
end

%if relative positioning (i.e. with master station)
if (goGNSS.isDD(mode) || mode == goGNSS.MODE_RT_NAV)
    %master station coordinates
    for i = 1 : nSol / (1 + state.isForwardBackwardKF())
        if (sum(abs(pos_M(:,i))) ~= 0)
            XM = pos_M(1,i);
            YM = pos_M(2,i);
            ZM = pos_M(3,i);
            
            %conversion from cartesian to geodetic coordinates
            [phiM(i), lamM(i), hM(i)] = cart2geod(XM, YM, ZM);
            
            %conversion from radians to degrees
            lamM(i) = lamM(i)*180/pi;
            phiM(i) = phiM(i)*180/pi;
        else
            lamM(i) = 0;
            phiM(i) = 0;
            hM(i) = 0;
        end
    end
end

pos = find(filerootOUT == '/');
if (isempty(pos))
    pos = find(filerootOUT == '\');
end
kml_name = checkPath(filerootOUT(pos(end)+1:end));

%file saving (Google Earth KML)
fid_kml = fopen([filerootOUT '.kml'], 'wt');
fprintf(fid_kml, '<?xml version="1.0" encoding="UTF-8"?>\n');
fprintf(fid_kml, '<kml xmlns="http://www.opengis.net/kml/2.2" xmlns:gx="http://www.google.com/kml/ext/2.2" xmlns:kml="http://www.opengis.net/kml/2.2" xmlns:atom="http://www.w3.org/2005/Atom">\n');
fprintf(fid_kml, '<Document>\n');
fprintf(fid_kml, '\t<name>%s</name>\n', [kml_name '.kml']);
fprintf(fid_kml, '\t<snippet>created by goGPS</snippet>\n');
fprintf(fid_kml, '\t\t<Style id="go1">\n');
fprintf(fid_kml, '\t\t\t<IconStyle>\n');
fprintf(fid_kml, '\t\t\t\t<color>%s</color>\n',good_point_colorR);
fprintf(fid_kml, '\t\t\t\t<scale>%.2f</scale>\n',scaleR);
fprintf(fid_kml, '\t\t\t\t<Icon>\n');
fprintf(fid_kml, '\t\t\t\t\t<href>%s</href>\n',iconR);
fprintf(fid_kml, '\t\t\t\t</Icon>\n');
fprintf(fid_kml, '\t\t\t</IconStyle>\n');
fprintf(fid_kml, '\t\t</Style>\n');
fprintf(fid_kml, '\t\t<Style id="go2">\n');
fprintf(fid_kml, '\t\t\t<IconStyle>\n');
fprintf(fid_kml, '\t\t\t\t<color>%s</color>\n',bad_point_colorR);
fprintf(fid_kml, '\t\t\t\t<scale>%.2f</scale>\n',scaleR);
fprintf(fid_kml, '\t\t\t\t<Icon>\n');
fprintf(fid_kml, '\t\t\t\t\t<href>%s</href>\n',iconR);
fprintf(fid_kml, '\t\t\t\t</Icon>\n');
fprintf(fid_kml, '\t\t\t</IconStyle>\n');
fprintf(fid_kml, '\t\t</Style>\n');
fprintf(fid_kml, '\t\t<Style id="go3">\n');
fprintf(fid_kml, '\t\t\t<IconStyle>\n');
fprintf(fid_kml, '\t\t\t\t<color>%s</color>\n',dyn_point_colorR);
fprintf(fid_kml, '\t\t\t\t<scale>%.2f</scale>\n',scaleR);
fprintf(fid_kml, '\t\t\t\t<Icon>\n');
fprintf(fid_kml, '\t\t\t\t\t<href>%s</href>\n',iconR);
fprintf(fid_kml, '\t\t\t\t</Icon>\n');
fprintf(fid_kml, '\t\t\t</IconStyle>\n');
fprintf(fid_kml, '\t\t</Style>\n');
fprintf(fid_kml, '\t\t<Style id="master">\n');
fprintf(fid_kml, '\t\t\t<IconStyle>\n');
fprintf(fid_kml, '\t\t\t\t<color>%s</color>\n',point_colorM);
fprintf(fid_kml, '\t\t\t\t<colorMode>normal</colorMode>\n');
fprintf(fid_kml, '\t\t\t\t<scale>%.2f</scale>\n',scaleM);
fprintf(fid_kml, '\t\t\t\t<Icon>\n');
fprintf(fid_kml, '\t\t\t\t\t<href>%s</href>\n',iconM);
fprintf(fid_kml, '\t\t\t\t</Icon>\n');
fprintf(fid_kml, '\t\t\t</IconStyle>\n');
fprintf(fid_kml, '\t\t\t<LabelStyle>\n');
fprintf(fid_kml, '\t\t\t\t<color>%s</color>\n',label_colorM);
fprintf(fid_kml, '\t\t\t\t<scale>%s</scale>\n',label_scaleM);
fprintf(fid_kml, '\t\t\t</LabelStyle>\n');
fprintf(fid_kml, '\t\t</Style>\n');
fprintf(fid_kml, '\t\t<Style id="ppos">\n');
fprintf(fid_kml, '\t\t\t<IconStyle>\n');
fprintf(fid_kml, '\t\t\t\t<color>%s</color>\n',point_colorP);
fprintf(fid_kml, '\t\t\t\t<colorMode>normal</colorMode>\n');
fprintf(fid_kml, '\t\t\t\t<scale>%.2f</scale>\n',scaleP);
fprintf(fid_kml, '\t\t\t\t<Icon>\n');
fprintf(fid_kml, '\t\t\t\t\t<href>%s</href>\n',iconP);
fprintf(fid_kml, '\t\t\t\t</Icon>\n');
fprintf(fid_kml, '\t\t\t</IconStyle>\n');
fprintf(fid_kml, '\t\t\t<LabelStyle>\n');
fprintf(fid_kml, '\t\t\t\t<color>%s</color>\n',label_colorP);
fprintf(fid_kml, '\t\t\t\t<scale>%s</scale>\n',label_scaleP);
fprintf(fid_kml, '\t\t\t</LabelStyle>\n');
fprintf(fid_kml, '\t\t</Style>\n');
fprintf(fid_kml, '\t\t<Style id="goLine1">\n');
fprintf(fid_kml, '\t\t\t<LineStyle>\n');
fprintf(fid_kml, '\t\t\t\t<color>%s</color>\n',line_colorR);
fprintf(fid_kml, '\t\t\t\t<width>%d</width>\n',line_widthR);
fprintf(fid_kml, '\t\t\t</LineStyle>\n');
fprintf(fid_kml, '\t\t</Style>\n');
if (flag_stopGOstop && goGNSS.isPP(mode)) %stop-go-stop and post-processing
    fprintf(fid_kml, '\t\t<Style id="goLine2">\n');
    fprintf(fid_kml, '\t\t\t<LineStyle>\n');
    fprintf(fid_kml, '\t\t\t\t<color>%s</color>\n',line_colorG);
    fprintf(fid_kml, '\t\t\t\t<width>%d</width>\n',line_widthG);
    fprintf(fid_kml, '\t\t\t</LineStyle>\n');
    fprintf(fid_kml, '\t\t</Style>\n');
end
if (goGNSS.isDD(mode) || mode == goGNSS.MODE_RT_NAV) %relative positioning
    for i = 1 : length(phiM)
        if (lamM(i) ~= 0 || phiM(i) ~= 0 || hM(i) ~= 0)
            if (i == 1) || (lamM(i)~=lamM(i-1) || phiM(i)~=phiM(i-1) || hM(i)~=hM(i-1))
                fprintf(fid_kml, '\t\t<Placemark>\n');
                fprintf(fid_kml, '\t\t\t<name>Master station</name>\n');
                fprintf(fid_kml, '\t\t\t<description><![CDATA[ <i>Latitude:</i> %.8f &#176;<br/> <i>Longitude:</i> %.8f &#176;<br/> <i>Elevation (ellips.):</i> %.1f m<br/>]]></description>\n',phiM(i),lamM(i),hM(i));
                fprintf(fid_kml, '\t\t\t<styleUrl>#master</styleUrl>\n');
                fprintf(fid_kml, '\t\t\t<Point>\n');
                fprintf(fid_kml, '\t\t\t\t<altitudeMode>%s</altitudeMode>\n',z_pos);
                fprintf(fid_kml, '\t\t\t\t<coordinates>%.8f,%.8f,%.3f</coordinates>\n',lamM(i),phiM(i),h_ortho(i));
                fprintf(fid_kml, '\t\t\t</Point>\n');
                fprintf(fid_kml, '\t\t</Placemark>\n');
            end
        end
    end
end
fprintf(fid_kml, '\t\t<Placemark>\n');
fprintf(fid_kml, '\t\t<name>Rover track</name>\n');
fprintf(fid_kml, '\t\t\t<styleUrl>#goLine1</styleUrl>\n');
fprintf(fid_kml, '\t\t\t<LineString>\n');
fprintf(fid_kml, '\t\t\t\t<coordinates>\n\t\t\t\t\t');
for i = 1 : nSol / (1 + state.isForwardBackwardKF())
    id = (state.getForwardBackwardKF() < 0) * (nSol/2) + (state.getForwardBackwardKF() > 0) * (nSol + 1 - 2 * i) + i;
    fprintf(fid_kml, '%.8f,%.8f,0 ',lam_KAL(id),phi_KAL(id));
end
fprintf(fid_kml, '\n\t\t\t\t</coordinates>\n');
fprintf(fid_kml, '\t\t\t</LineString>\n');
fprintf(fid_kml, '\t\t</Placemark>\n');
if (flag_stopGOstop && flag_var_dyn_model && mode == goGNSS.MODE_PP_KF_CP_DD)
    
    [P1Lat, P1Lon] = cart2geod(P1_GLB(1), P1_GLB(2), P1_GLB(3));
    [P2Lat, P2Lon] = cart2geod(P2_GLB(1), P2_GLB(2), P2_GLB(3));
    
    fprintf(fid_kml, '\t\t<Placemark>\n');
    fprintf(fid_kml, '\t\t<name>Estimated direction</name>\n');
    fprintf(fid_kml, '\t\t\t<styleUrl>#goLine2</styleUrl>\n');
    fprintf(fid_kml, '\t\t\t<LineString>\n');
    fprintf(fid_kml, '\t\t\t\t<coordinates>\n\t\t\t\t\t');
    fprintf(fid_kml, '%.8f,%.8f,0 ',P1Lon*180/pi,P1Lat*180/pi);
    fprintf(fid_kml, '%.8f,%.8f,0 ',P2Lon*180/pi,P2Lat*180/pi);
    fprintf(fid_kml, '\n\t\t\t\t</coordinates>\n');
    fprintf(fid_kml, '\t\t\t</LineString>\n');
    fprintf(fid_kml, '\t\t</Placemark>\n');
end
fprintf(fid_kml, '\t\t<Folder>\n');
fprintf(fid_kml, '\t\t<name>Rover positioning</name>\n');
for i = 1 : nSol / (1 + state.isForwardBackwardKF())
    id = (state.getForwardBackwardKF() < 0) * (nSol/2) + (state.getForwardBackwardKF() > 0) * (nSol + 1 - 2 * i) + i;
    
    fprintf(fid_kml, '\t\t<Placemark>\n');
    if (pivot_OUT(id) == 0)
        fprintf(fid_kml, '\t\t\t<styleUrl>#go3</styleUrl>\n');
    elseif ~(mode == goGNSS.MODE_PP_BLK_CP_DD_STATIC) && (KHDOP(id)>KHDOP_thres)
        fprintf(fid_kml, '\t\t\t<styleUrl>#go2</styleUrl>\n');
    else
        fprintf(fid_kml, '\t\t\t<styleUrl>#go1</styleUrl>\n');
    end
    fprintf(fid_kml, '\t\t\t<Point>\n');
    fprintf(fid_kml, '\t\t\t\t<altitudeMode>%s</altitudeMode>\n',z_pos);
    fprintf(fid_kml, '\t\t\t\t<coordinates>%.8f,%.8f,%.3f</coordinates>\n',lam_KAL(i),phi_KAL(i),h_KAL(i));
    fprintf(fid_kml, '\t\t\t</Point>\n');
    fprintf(fid_kml, '\t\t</Placemark>\n');
end
fprintf(fid_kml, '\t\t</Folder>\n');

if (mode_vinc == 0) && ((mode == goGNSS.MODE_PP_KF_C_SA) || (mode == goGNSS.MODE_PP_KF_CP_SA) || (mode == goGNSS.MODE_PP_KF_C_DD) || (mode == goGNSS.MODE_PP_KF_CP_DD) || (mode == goGNSS.MODE_RT_NAV) || (mode == goGNSS.MODE_PP_BLK_CP_DD_STATIC))
    if (o1 == 1) && (nSol ~= 0)
        %static positioning coordinates
        phiP = phi_KAL(end);
        lamP = lam_KAL(end);
        hP   = h_KAL(end);
        
        fprintf(fid_kml, '\t\t<Placemark>\n');
        fprintf(fid_kml, '\t\t\t<name>Static positioning</name>\n');
        fprintf(fid_kml, '\t\t\t<description><![CDATA[ <i>Latitude:</i> %.8f &#176;<br/> <i>Longitude:</i> %.8f &#176;<br/> <i>Elevation (ellips.):</i> %.1f m<br/>]]></description>\n',phiP,lamP,hP);
        fprintf(fid_kml, '\t\t\t<styleUrl>#ppos</styleUrl>\n');
        fprintf(fid_kml, '\t\t\t<Point>\n');
        fprintf(fid_kml, '\t\t\t\t<altitudeMode>%s</altitudeMode>\n',z_pos);
        fprintf(fid_kml, '\t\t\t\t<coordinates>%.8f,%.8f,%.3f</coordinates>\n',lamP,phiP,hP);
        fprintf(fid_kml, '\t\t\t</Point>\n');
        fprintf(fid_kml, '\t\t</Placemark>\n');
    end
end
fprintf(fid_kml, '</Document>\n</kml>');
fclose(fid_kml);
